# ==============================================================================
# Importación de paquetes y preparación ---
# ==============================================================================
library(fs)
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(ggtime)
library(tsibble)
library(feasts)
library(fable)
library(tseries)
library(FinTS)
library(lmtest)
library(urca) 

ARIMA <- fable::ARIMA

# Fijar ruta de referencia
here::i_am("BoxJenkins_20261/BoxJenkins.R")
BASE_DIR <- here::here("BoxJenkins_20261")
DATA_DIR <- fs::path(BASE_DIR, "datos")
ruta_datos <- fs::path(DATA_DIR, "serie.csv")

# --- Funciones auxiliares ---
grilla <- function(..., nrow, ncol) {
  graficos <- list(...)
  if (length(graficos) > nrow * ncol) stop("La cantidad de graficos supera el tamano de la grilla.")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow, ncol = ncol)))
  for (i in seq_along(graficos)) {
    fila <- ceiling(i / ncol)
    columna <- ((i - 1) %% ncol) + 1
    print(graficos[[i]], vp = grid::viewport(layout.pos.row = fila, layout.pos.col = columna))
  }
  grid::popViewport()
}

# ==============================================================================
# Carga y Transformación de Datos ---
# ==============================================================================
datos_base <- read_delim(
  ruta_datos,
  delim = ",",
  col_names = c("fecha", "precio"),  
  locale = locale(decimal_mark = "."),
  show_col_types = FALSE,
  skip = 1 
)

fechas_datos_base <- yearmonth(seq(
  from = as.Date("2010-05-01"),
  by = "month",
  length.out = nrow(datos_base)
))

datos_tbl <- datos_base |>
  mutate(fecha = fechas_datos_base) |>
  as_tsibble(index = fecha)

datos_serie <- datos_tbl |>
  mutate(precio = as.numeric(precio)) |>
  filter(!is.na(precio))

datos_vector <- datos_serie$precio

# ==============================================================================
# PASO 1: Identificación ====
# ==============================================================================
print(
  ggtime::autoplot(datos_tbl, precio) +
    ggtitle("Precio internacional, 2010-2026") +
    xlab("Fecha") + ylab("Precio(USD)") + theme_light()
)

# Análisis de Estacionariedad en la Serie Original
grafico_fac_datos <- datos_tbl |> ACF(precio, lag_max = 15) |> ggtime::autoplot() + ggtitle("FAC del precio") + ylim(-1, 1) + theme_light()
grafico_facp_datos <- datos_tbl |> PACF(precio, lag_max = 15) |> ggtime::autoplot() + ggtitle("FACP del precio") + ylim(-1, 1) + theme_light()
grilla(grafico_fac_datos, grafico_facp_datos, nrow = 1, ncol = 2)

# --- Tests de Raíz Unitaria ---
adf_result <- ur.df(datos_serie$precio, type = "drift", selectlags = "AIC")
cat("=== Test ADF ===\n")
print(summary(adf_result))

kpss_result <- ur.kpss(datos_serie$precio, type = "mu", lags = "short")
cat("\n=== Test KPSS ===\n")
print(summary(kpss_result))

# --- Transformaciones Estacionarias ---
datos_serie <- datos_serie |>
  mutate(
    diff_datos = difference(precio),
    log_datos = log(precio),
    log_diff = difference(log_datos)
  )

# --- Gráfica de la serie de tiempo en primera diferencia (Precios diferenciados)
print(
  datos_serie |>
    filter(!is.na(diff_datos)) |>
    ggtime::autoplot(diff_datos) +
    ggtitle("Serie Diferenciada (Primera Diferencia del Precio)") +
    xlab("Fecha") +
    ylab("D(Precio)") + 
    theme_light()
)

# --- Gráfica de la serie de tiempo en diferencia logarítmica (Retornos/Tasas de cambio)
print(
  datos_serie |>
    filter(!is.na(log_diff)) |>
    ggtime::autoplot(log_diff) +
    ggtitle("Diferencia del Logaritmo de la Serie de Precios") +
    xlab("Fecha") +
    ylab("D(Log(Precio))") + 
    theme_light()
)

# --- Gráficas de Estacionariedad (Logaritmo Diferenciado) ---
grafico_fac_log_diff <- datos_serie |>
  filter(!is.na(log_diff)) |>
  ACF(log_diff, lag_max = 24) |>
  ggtime::autoplot() + ggtitle("FAC de la diferencia del logaritmo") + ylim(-1, 1) + theme_light()

grafico_facp_log_diff <- datos_serie |>
  filter(!is.na(log_diff)) |>
  PACF(log_diff, lag_max = 24) |>
  ggtime::autoplot() + ggtitle("FACP de la diferencia del logaritmo") + ylim(-1, 1) + theme_light()

grilla(grafico_fac_log_diff, grafico_facp_log_diff, nrow = 1, ncol = 2)

# ==============================================================================
# PASO 2: Estimación (Selección e Inferencia) =========================
# ==============================================================================

# 1. Definición de la función para armar el data frame de criterios
arma_seleccion_df = function(ts_object, AR.m, MA.m, d, bool_trend, metodo){
  index = 1
  df = data.frame(p = double(), d = double(), q = double(), AIC = double(), BIC = double())
  for (p in 0:AR.m) {
    for (q in 0:MA.m)  {
      posible_modelo <- tryCatch({
        arima(ts_object, order = c(p, d, q), include.mean = bool_trend, method = metodo)
      }, error = function(e) { NULL })
      if (!is.null(posible_modelo)) {
        df[index,] = c(p, d, q, AIC(posible_modelo), BIC(posible_modelo))
      } else {
        df[index,] = c(p, d, q, NA, NA)
      }
      index = index + 1
    }
  }  
  return(df)
}

# 2. Funciones auxiliares para extraer los modelos con menor criterio
arma_min_AIC = function(df){
  df2 = df %>% filter(AIC == min(AIC, na.rm = TRUE))
  return(df2)
}

arma_min_BIC = function(df){
  df2 = df %>% filter(BIC == min(BIC, na.rm = TRUE))
  return(df2)
}

# 3. Preparación de datos y ejecución de la búsqueda (hasta orden 6)
datos_para_seleccion <- datos_serie$log_diff[!is.na(datos_serie$log_diff)]
tabla_criterios <- arma_seleccion_df(datos_para_seleccion, AR.m = 6, MA.m = 6, d = 0, TRUE, "ML")

# Abre la tabla completa en una pestaña interactiva de RStudio
View(tabla_criterios) 

# También puedes imprimir las primeras filas en la consola para verificar
cat("\n=== Primeras filas de la tabla de criterios ===\n")
print(head(tabla_criterios, 15))

# Seleccionar de forma automática los mejores modelos de la tabla
mejor_modelo_AIC <- arma_min_AIC(tabla_criterios)
mejor_modelo_BIC <- arma_min_BIC(tabla_criterios)

# Imprimir los resultados analíticos exactos en la consola
cat("\n==============================================\n")
cat("  RESULTADOS DE SELECCIÓN DE MODELO")
cat("\n==============================================\n")
cat("\n=== MEJOR MODELO SEGÚN CRITERIO AIC ===\n")
print(mejor_modelo_AIC)

cat("\n=== MEJOR MODELO SEGÚN CRITERIO BIC (PARSIMONIA) ===\n")
print(mejor_modelo_BIC)
cat("\n==============================================\n")

# ==============================================================================

# Estimación final en Fable para un ARIMA(0,1,1) sin constante
# (Recuerda cambiar el pdq(0,1,1) si tras ver la tabla prefieres usar el orden del AIC o BIC)
fit_modelo <- datos_serie |>
  select(fecha, precio) |>
  model(modelo_optimo = fable::ARIMA(log(precio) ~ 0 + pdq(0, 1, 1) + PDQ(0, 0, 0)))

cat("\n=== Reporte de Estimación Fable ===\n")
report(fit_modelo)
# ==============================================================================
# PASO 3: Validación de supuestos =========================
# ==============================================================================
residuales_tbl <- fit_modelo |> residuals() |> filter(!is.na(.resid))
residuales <- residuales_tbl$.resid

# --- Gráficas e Inferencia de Ruido Blanco ---
grafico_fac_residuales <- residuales_tbl |> ACF(.resid, lag_max = 24) |> ggtime::autoplot() + ggtitle("FAC de los residuales") + ylim(-1, 1) + theme_light()
grafico_facp_residuales <- residuales_tbl |> PACF(.resid, lag_max = 24) |> ggtime::autoplot() + ggtitle("FACP de los residuales") + ylim(-1, 1) + theme_light()
grilla(grafico_fac_residuales, grafico_facp_residuales, nrow = 1, ncol = 2)

cat("\n--- Prueba Ljung-Box (H0: Errores No Autocorrelacionados) ---\n")
rezagos_ljung_box <- c(6, 12, 18, 24)
prueba_ljung_box <- lapply(rezagos_ljung_box, function(lag) {
  prueba <- Box.test(residuales, lag = lag, type = "Ljung-Box", fitdf = 1) # fitdf = 1 dado el parámetro ma1
  data.frame(lb_stat = as.numeric(prueba$statistic), lb_pvalue = prueba$p.value)
}) |> bind_rows()
rownames(prueba_ljung_box) <- rezagos_ljung_box
print(round(prueba_ljung_box, 6))

# --- Diagnóstico de Homocedasticidad (Efectos ARCH) ---
residuales_tbl <- residuales_tbl |> mutate(residuo_cuadrado = .resid^2)
arch_test <- FinTS::ArchTest(residuales, lags = 12)
cat("\n--- Prueba ARCH de heterocedasticidad (H0: Errores Homocedásticos) ---\n")
print(arch_test)

# --- Diagnóstico de Normalidad ---
car::qqPlot(residuales, main = "Q-Q plot de los residuos", xlab = "Cuantiles teóricos", ylab = "Residuales")

jb_test <- tseries::jarque.bera.test(residuales)
cat("\n--- Prueba Jarque-Bera (H0: Errores Siguen una Distribución Normal) ---\n")
print(jb_test)

# ==============================================================================
# PASO 4: Pronóstico =========================
# ==============================================================================
pronostico <- fit_modelo |> forecast(h = 12)

tabla_pronostico <- pronostico |> as_tibble() |> select(fecha, .mean) |> rename(pronostico = .mean)
intervalos <- pronostico |> hilo(level = 95) |> fabletools::unpack_hilo(`95%`) |> as_tibble() |> 
  select(fecha, `95%_lower`, `95%_upper`) |> rename(limite_inferior = `95%_lower`, limite_superior = `95%_upper`)

tabla_pronostico <- left_join(tabla_pronostico, intervalos, by = "fecha")
cat("\n=== Tabla de Pronóstico (Nivel original de precios) ===\n")
print(tabla_pronostico)

print(
  ggtime::autoplot(datos_serie, precio) +
    ggtime::autolayer(pronostico, level = 95) +
    ggtitle("Pronóstico Ex-Ante 12 Meses Adelante (ARIMA(0,1,1))") +
    xlab("Fecha") + ylab("Precio") + theme_light()
)
 
