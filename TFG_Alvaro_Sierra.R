# TRABAJO DE FIN DE GRADO

# Título: Econometría espacial aplicada al estudio de las disparidades regionales de la renta per cápita en España

# Autor: Álvaro Sierra García
# Directora: Beatriz Manotas Hidalgo

# Julio de 2025

# Versión de R utilizada: 4.4.2

##################################################################################################################

##### Carga de librerías
library(eurostat)
library(giscoR)
library(dplyr)
library(stringr)
library(units)
library(sf)
library(spdep)
library(spatialreg)
library(ggplot2)
library(tibble)

##### Descarga de datos

## Configuraciones previas para descarga de Eurostat
options(eurostat_cache = TRUE, eurostat_update = TRUE, eurostat_httpfix = TRUE)
Sys.setenv(LANG = "en_US.UTF-8")

## Función para estandarizar columna año
fix_time <- function(df) { if ("TIME_PERIOD" %in% names(df)) dplyr::rename(df, año = TIME_PERIOD) else df }

# PIB total
pib <- get_eurostat("nama_10r_3gdp", time_format = "num") %>%
  fix_time() %>%
  filter(str_detect(geo, "^ES"),
         unit == "MIO_EUR",
         nchar(geo) == 5) %>% # NUTS-3
  transmute(geo, año, pib = values)

# Población
poblacion <- get_eurostat("demo_r_pjangrp3", time_format = "num") %>%
  fix_time() %>%
  filter(str_detect(geo, "^ES"), sex == "T") %>%
  group_by(geo, año) %>%
  summarise(poblacion = sum(values, na.rm = TRUE), .groups = "drop") %>%
  filter(nchar(geo) == 5) # NUTS-3

# Paro
paro_nuts2 <- get_eurostat("lfst_r_lfu3rt", time_format = "num") %>%
  fix_time() %>%
  filter(str_detect(geo, "^ES"),
         sex == "T", age == "Y15-74", isced11 == "TOTAL",
         nchar(geo) == 4) %>% # NUTS-2 (no hay datos NUTS-3, proyección a NUTS-3 más adelante)
  transmute(nuts2 = geo, año, tasa_paro = values)

# GISCO (cartografía)
gisco_nuts3_esp <- gisco_get_nuts(year = 2021, resolution = "1", nuts_level = 3) %>%
  st_make_valid() %>%
  filter(CNTR_CODE == "ES") %>%
  select(geo = NUTS_ID, NAME_LATN, geometry)

###############################################################
# 4.2: Análisis de autocorrelación espacial
###############################################################

## 1. Construcción DataFrame

# Selección del año a utilizar
años_pib_poblacion <- intersect(unique(pib$año), unique(poblacion$año)) # NUTS-3
años_paro <- unique(paro_nuts2$año)  # NUTS-2
años_comun  <- intersect(años_pib_poblacion, años_paro)

año_uso <- max(años_comun, na.rm = TRUE)
año_uso

# Proyección paro a NUTS-3
gisco_nuts3_a_nuts2 <- gisco_nuts3_esp %>%
  st_drop_geometry() %>%
  transmute(geo, nuts2 = substr(geo, 1, 4))

paro_nuts3 <- gisco_nuts3_a_nuts2 %>%
  left_join(paro_nuts2 %>% filter(año == año_uso), by = "nuts2") %>%
  transmute(geo, año = año_uso, tasa_paro)

# DataFrame final para el año 2022
df <- pib %>% filter(año == año_uso) %>%
  inner_join(poblacion  %>% filter(año == año_uso), by = c("geo","año")) %>%
  mutate(pib_pc = (pib * 1e6) / poblacion) %>%
  inner_join(paro_nuts3, by = c("geo","año")) %>%
  filter(geo %in% gisco_nuts3_esp$geo)


## 2. Datos espaciales, vecindades y pesos

# Proyección métrica y área (km2)
gisco_nuts3_esp <- st_transform(gisco_nuts3_esp, 3035)
gisco_nuts3_esp$area_km2 <- set_units(st_area(gisco_nuts3_esp), km^2) |> drop_units()

# Datos espaciales
dat_sp <- gisco_nuts3_esp %>%
  left_join(df, by = "geo") %>%
  mutate(
    lpibpc = log(pib_pc),
    densidad = poblacion / area_km2,
    año = año_uso
  ) %>%
  filter(!is.na(lpibpc), !is.na(tasa_paro), !is.na(densidad)) %>%
  arrange(geo)

# Vecindades y listas de pesos
nb_queen <- poly2nb(dat_sp, queen = TRUE,  row.names = dat_sp$geo)
nb_rook  <- poly2nb(dat_sp, queen = FALSE, row.names = dat_sp$geo)

lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
lw_rook  <- nb2listw(nb_rook,  style = "W", zero.policy = TRUE)

## 3. I de Moran y C de Geary (global)
set.seed(12345)

tests_globales <- function(var, listw, listw_name, nperm = 999) {
  x <- dat_sp[[var]]
  mt <- moran.test(x, listw, zero.policy = TRUE)
  mc <- moran.mc(x, listw, nsim = nperm, zero.policy = TRUE)
  gt <- geary.test(x, listw, zero.policy = TRUE)
  data.frame(
    variable   = var,
    matriz_w   = listw_name,
    moran_I    = unname(mt$estimate[1]),
    moran_EI   = unname(mt$estimate[2]),
    moran_var  = unname(mt$estimate[3]),
    moran_p    = mt$p.value,
    moran_mc_I = unname(mc$statistic),
    moran_mc_p = mc$p.value,
    geary_C    = unname(gt$estimate[1]),
    geary_EC   = unname(gt$estimate[2]),
    geary_var  = unname(gt$estimate[3]),
    geary_p    = gt$p.value
  )
}

resumen_global <- bind_rows(
  tests_globales("lpibpc",    lw_queen, "queen"),
  tests_globales("tasa_paro",lw_queen, "queen"),
  tests_globales("lpibpc",    lw_rook,  "rook"),
  tests_globales("tasa_paro",lw_rook,  "rook")
)
resumen_global


## 3. LISA (Moran local) para lpibpc con queen
y   <- dat_sp$lpibpc
y_z <- as.numeric(scale(y))
Wy  <- lag.listw(lw_queen, y_z, zero.policy = TRUE)

li <- localmoran(y, lw_queen, zero.policy = TRUE)  # Ii, E(Ii), Var, Z, p
lisa_df <- dat_sp %>% st_drop_geometry() %>%
  mutate(
    y_z = y_z,
    Wy_z = Wy,
    Ii   = li[,1],
    Ii_z = li[,4],
    Ii_p = li[,5]
  )

clasificacion <- ifelse(lisa_df$y_z >= 0 & lisa_df$Wy_z >= 0, "High–High",
               ifelse(lisa_df$y_z <  0 & lisa_df$Wy_z <  0, "Low–Low",
                      ifelse(lisa_df$y_z >= 0 & lisa_df$Wy_z <  0, "High–Low", "Low–High")))
significativo  <- ifelse(lisa_df$Ii_p <= 0.05, "Significativo", "No significativo")

lisa_df <- lisa_df %>% mutate(clasificacion = clasificacion, significativo = significativo)
lisa_df$significativo[is.na(lisa_df$significativo)] <- "No significativo"

lisa_sp <- dat_sp %>%
  select(geo, NAME_LATN, geometry) %>%
  left_join(lisa_df %>% select(geo, y_z, Wy_z, Ii, Ii_p, clasificacion, significativo), by = "geo")

## 4. Gráficos
p_scatter <- ggplot(lisa_df, aes(x = y_z, y = Wy_z)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(title = paste0("Gráfico de dispersión de Moran: log(PIB per cápita) - ", as.character(año_uso)),
       x = "z(log PIB per cápita)", y = "W z(log PIB per cápita)") +
  theme_minimal()
ggsave(filename = paste0("moran_", as.character(año_uso), ".png"),
       plot = p_scatter, width = 7, height = 5, dpi = 300)

lisa_sp$cluster <- factor(ifelse(lisa_sp$significativo=="Significativo", lisa_sp$clasificacion, "No significativo"),
                          levels = c("High–High","Low–Low","High–Low","Low–High","No significativo"))
p_lisa <- ggplot(lisa_sp) +
  geom_sf(aes(fill = cluster), linewidth = 0.15, color = "grey40") +
  labs(title = paste0("Mapa LISA: log(PIB per cápita) - ", as.character(año_uso)), fill = "Cluster") +
  theme_minimal()
ggsave(filename = paste0("LISA_", as.character(año_uso), ".png"),
       plot = p_lisa, width = 7, height = 6, dpi = 300)


###############################################################
# 4.3: Estimación de modelos espaciales: SAR y SEM
###############################################################

## 1. OLS
ols <- lm(lpibpc ~ tasa_paro + densidad, data = dat_sp)
summary(ols)

# LM test clásicos y robustos (diagnóstico de especificación espacial)
lm_tests <- lm.LMtests(ols, lw_queen, test = c("LMlag","LMerr","RLMlag","RLMerr","SARMA"),
                       zero.policy = TRUE)
lm_tests

# Moran en residuos OLS
moran_ols <- moran.test(residuals(ols), lw_queen, zero.policy = TRUE)
moran_ols


## 2. SAR y SEM con W = queen
sar_q <- lagsarlm(lpibpc ~ tasa_paro + densidad, data = dat_sp,
                  listw = lw_queen, method = "eigen", zero.policy = TRUE)
summary(sar_q)


sem_q <- errorsarlm(lpibpc ~ tasa_paro + densidad, data = dat_sp,
                    listw = lw_queen, method = "eigen", zero.policy = TRUE)
summary(sem_q)

# Moran en residuos de SAR y SEM
moran_sar <- moran.test(residuals(sar_q), lw_queen, zero.policy = TRUE)
moran_sar
moran_sem <- moran.test(residuals(sem_q), lw_queen, zero.policy = TRUE)
moran_sem

## 3. AIC comparado
AIC = c(AIC(ols), AIC(sar_q), AIC(sem_q))
AIC

## 4. Impactos del SAR
imp_q   <- impacts(sar_q, listw = lw_queen, R = 1000)
imp_q

# Monte Carlo
imp_sum <- summary(imp_q, zstats = TRUE, short = FALSE)
imp_MC <- as.data.frame(imp_sum$tot$statistics) %>%
  rownames_to_column("variable")
imp_MC


###############################################################
# 4.4: Comprobación de robustez y contraste de especificaciones
###############################################################

##### 4.4.1 Grado medio de los grafos queen y rook

grado_q <- card(nb_queen)
grado_r <- card(nb_rook)

resumen_grados <- tibble(
  esquema = c("queen","rook"),
  n_nodos = c(length(grado_q), length(grado_r)),
  aislados = c(sum(grado_q == 0), sum(grado_r == 0)),
  grado_medio = c(mean(grado_q), mean(grado_r)),
  mediana = c(median(grado_q), median(grado_r)),
  p25 = c(quantile(grado_q, .25), quantile(grado_r, .25)),
  p75 = c(quantile(grado_q, .75), quantile(grado_r, .75)),
  min = c(min(grado_q), min(grado_r)),
  max = c(max(grado_q), max(grado_r)),
  sd  = c(sd(grado_q), sd(grado_r))
)
resumen_grados

##### 4.4.2 Estabilidad de residuos estandarizados por provincia
# SEM con rook
sem_r <- errorsarlm(lpibpc ~ tasa_paro + densidad, data = dat_sp,
                    listw = lw_rook, method = "eigen", zero.policy = TRUE)

# Estandarizar residuos
res_q <- scale(residuals(sem_q))
res_r <- scale(residuals(sem_r))

# Estadísticos de estabilidad
corr_qr <- cor(res_q, res_r, use = "complete.obs")
corr_qr
rmse_qr <- sqrt(mean((res_q - res_r)^2, na.rm = TRUE))
rmse_qr

# Gráfico de residuos estandarizados (Q vs R)
res_df <- data.frame(
  geo = dat_sp$geo,
  res_q = as.numeric(res_q),
  res_r  = as.numeric(res_r)
)

p_res <- ggplot(res_df, aes(x = res_q, y = res_r)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(title = paste0("Residuos estandarizados (queen vs rook) - ", as.character(año_uso)),
       x = "Residuos estandarizados (queen)",
       y = "Residuos estandarizados (rook)") +
  theme_minimal()
ggsave(filename = paste0("residuos_", as.character(año_uso), ".png"),
       plot = p_res, width = 7, height = 5, dpi = 300)

# Top 10 provincias
top_10_absoluto <- res_df %>%
  mutate(abs_q = abs(res_q), abs_r = abs(res_r)) %>%
  arrange(desc(abs_q)) %>%
  select(geo, res_q, res_r) %>%
  head(10)
top_10_absoluto

##### 4.4.3 Impactos directos e indirectos de la tasa de paro bajo distintos esquemas de vecindad

# Estimar SAR con rook
sar_r <- lagsarlm(lpibpc ~ tasa_paro + densidad, data = dat_sp,
                  listw = lw_rook, method = "eigen", zero.policy = TRUE)

# Impactos puntuales (directo, indirecto, total) en cada esquema
imp_q
imp_r <- impacts(sar_r, listw = lw_rook,  R = 1000)
imp_r