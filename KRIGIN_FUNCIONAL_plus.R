# ===============================
# 1. Cargar librerías necesarias
# ===============================
library(sf)          # Para manejar datos geoespaciales
library(sp)          # Funciones espaciales básicas
library(ggplot2)     # Visualización gráfica
library(dplyr)       # Manipulación de datos
library(fda)         # Análisis funcional de datos
library(gstat)       # Geoestadística y kriging
library(plotly)      # Visualización 3D interactiva
library(viridis)     # Escalas de color
library(mvtnorm)     # Para manejo de distribuciones multivariadas y covarianzas

# ===============================
# 2. Cargar y preparar los datos 
# ===============================
# Leer datos geoespaciales desde archivo JSON
gdf <- st_read("G:/Univalle/puntos_proyect.json") 

# Verificación de estructura de datos
if(!"DATE" %in% names(gdf)) stop("La columna DATE no existe en los datos")
if(nrow(gdf) == 0) stop("Datos geoespaciales vacíos")

# Uso de coordenadas originales
gdf_df <- gdf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(X = X, Y = Y) %>%
  mutate(DATE = as.Date(gdf$DATE))  # Conversión a tipo Date

# Visualización de datos transformados
ggplot(gdf_df, aes(X, Y, color = DATE)) +
  geom_point() +
  scale_color_viridis(trans = "date") +  # Escala de color por fecha
  labs(title = "Datos de línea costera - Coordenadas originales") +
  theme_minimal()

# ===============================
# 3. Análisis Funcional 
# ===============================
# Definición de puntos comunes para interpolación
X_common <- seq(min(gdf_df$X), max(gdf_df$X), length.out = 100)

# Función de interpolación segura con manejo de errores
safe_interp <- function(x, y, xout) {
  tryCatch({
    approx(x, y, xout = xout, rule = 2)$y
  }, error = function(e) {
    warning("Error en interpolación: ", e$message)
    rep(NA, length(xout))
  })
}

# Interpolación 
Y_matrix <- gdf_df %>%
  group_by(DATE) %>%
  group_map(~ safe_interp(.x$X, .x$Y, X_common)) %>%
  purrr::discard(~any(is.na(.x))) %>%
  do.call(cbind, .)

# Función para selección óptima de bases B-spline usando validación cruzada
cv_basis <- function(data, max_basis = 15) {
  cv_errors <- sapply(5:max_basis, function(nb) {
    basis <- create.bspline.basis(range(X_common), nbasis = nb)
    fdPar <- fdPar(basis)
    mean(smooth.basis(X_common, data, fdPar)$SSE)
  })
  which.min(cv_errors) + 4
}

# Creación de base óptima y objeto funcional
optimal_basis <- cv_basis(Y_matrix)
fd_basis <- create.bspline.basis(range(X_common), nbasis = optimal_basis)
fd_obj <- Data2fd(X_common, Y_matrix, fd_basis)

# Visualización de curvas funcionales
plot(fd_obj, main = "Representación Funcional de la Línea de Costa")

# ===============================
# 4. FPCA y Selección de Componentes 
# ===============================
# Análisis de componentes principales funcionales
fpca <- pca.fd(fd_obj, nharm = 10)

# Asignación de nombres a los scores
colnames(fpca$scores) <- paste0("score.", 1:ncol(fpca$scores))

# Selección dinámica de componentes principales
cumvar <- cumsum(fpca$values)/sum(fpca$values)
n_components <- which(cumvar >= 0.95)[1]
fpca$harmonics <- fpca$harmonics[1:n_components]
fpca$scores <- fpca$scores[, 1:n_components, drop = FALSE]  # Mantener estructura de matriz

# ===============================
# 5. MODELADO TEMPORAL CON KRIGING
# ===============================
# 5.1 Preparación de datos espaciales 
# ===============================
# Conversión de fechas a numérico (usando años decimales)
dates_numeric <- as.numeric(format(gdf_df$DATE, "%Y")) + 
  as.numeric(format(gdf_df$DATE, "%j"))/365

data_vario <- data.frame(
  time = dates_numeric,
  y = jitter(0, amount = 1e-6),
  fpca$scores
) %>% 
  distinct(time, .keep_all = TRUE)
coordinates(data_vario) <- ~ time + y

# Validación de estructura espacial
if(!inherits(data_vario, "SpatialPointsDataFrame")) {
  stop("Error en conversión a objeto espacial")
}
# ===============================
# 5.2 Función de ajuste de variogramas
# ===============================
fit_variogram_robust <- function(formula, sp_data) {
  environment(formula) <- environment()
  nombre_variable <- all.vars(formula)[1]
  
  if(!nombre_variable %in% names(sp_data)) {
    stop("Variable ", nombre_variable, " no encontrada")
  }
  
  valores <- sp_data[[nombre_variable]]
  var_inicial <- var(valores, na.rm = TRUE)
  rango_temporal <- diff(range(coordinates(sp_data)[,1], na.rm = TRUE))
  
  message("Variable: ", nombre_variable, " - Varianza inicial: ", var_inicial, 
          " - Rango temporal: ", rango_temporal)
  
  if(var_inicial == 0) {
    message("Varianza cero en ", nombre_variable, ". Usando modelo Nugget.")
    return(vgm(0, "Nug", 0))
  }
  
  # Modelos candidatos extendidos
  modelos <- list(
    vgm(0.6 * var_inicial, "Exp", rango_temporal * 2, 0.4 * var_inicial),
    vgm(0.6 * var_inicial, "Ste", rango_temporal * 2, 0.4 * var_inicial, kappa = 0.5),
    vgm(0.5 * var_inicial, "Gau", rango_temporal / 2, 0.5 * var_inicial),
    vgm(0.7 * var_inicial, "Mat", rango_temporal / 3, 0.3 * var_inicial, kappa = 1.5),
    vgm(0.7 * var_inicial, "Mat", rango_temporal / 3, 0.3 * var_inicial, kappa = 2.5),
    vgm(var_inicial, "Nug", 0)
  )
  
  # Cálculo del variograma empírico (sin cambios importantes)
  ancho <- rango_temporal / 15
  vgm_emp <- tryCatch({
    variogram(formula, data = sp_data, width = ancho, cutoff = rango_temporal * 0.8)
  }, error = function(e) {
    message("Error en cálculo del variograma empírico: ", e$message)
    return(NULL)
  })
  
  if(is.null(vgm_emp)) return(NULL)
  
  ajustes <- lapply(modelos, function(m) {
    tryCatch(
      fit.variogram(vgm_emp, model = m, fit.method = 7),
      error = function(e) {
        message("Error ajustando modelo ", m$model[2], ": ", e$message)
        return(NULL)
      }
    )
  })
  
  sse <- sapply(ajustes, function(a) {
    if(!is.null(a)) {
      fitted_vals <- tryCatch({
        variogramLine(a, maxdist = max(vgm_emp$dist), n = nrow(vgm_emp))$gamma
      }, error = function(e) {
        message("Error en cálculo del SSE para modelo ", a$model[2], ": ", e$message)
        return(rep(Inf, nrow(vgm_emp)))
      })
      emp_vals <- vgm_emp$gamma
      sum((fitted_vals - emp_vals)^2)
    } else {
      Inf
    }
  })
  
  indices_validos <- which(sapply(ajustes, function(a) {
    !is.null(a) && a$range[2] > 0
  }))
  
  if(length(indices_validos) > 0) {
    sse_validos <- sse[indices_validos]
    mejor_index <- indices_validos[which.min(sse_validos)]
    mejor_ajuste <- ajustes[[mejor_index]]
  } else {
    warning("No se encontró modelo con estructura positiva para ", nombre_variable, 
            ". Se usará modelo Nugget.")
    mejor_ajuste <- vgm(var_inicial, "Nug", 0)
  }
  
  return(mejor_ajuste)
}

# ===============================
# 5.3 Ajuste de variogramas por componente
# ===============================
# Selección de componentes con varianza significativa
varianzas <- apply(fpca$scores, 2, var)
umbral_varianza <- 0.01 * max(varianzas)  # Umbral del 1% de varianza máxima
componentes_activos <- which(varianzas > umbral_varianza)

# Ajuste de variogramas para cada componente activo
fit_list <- lapply(componentes_activos, function(i) {
  componente <- paste0("score.", i)
  formula <- reformulate("1", response = componente)
  
  message("Ajustando componente: ", componente)
  
  tryCatch({
    fit_variogram_robust(formula, data_vario)
  }, error = function(e) {
    message("Error grave en ", componente, ": ", e$message)
    return(NULL)
  })
})

# Validación de ajustes
if(any(sapply(fit_list, is.null))) {
  stop("Error en ajuste de componentes: ", 
       paste(componentes_activos[sapply(fit_list, is.null)], collapse = ", "))
}

# Información de diagnóstico
print(paste("Número de componentes activos:", length(componentes_activos)))
print(head(data_vario))
print(summary(varianzas))

# 5.4 Validación y visualización
# ---------------------------------------------------------------
x11()
if(length(componentes_activos) > 0) {
  # Configuración de ventana gráfica múltiple
  dev.new()  # Nueva ventana gráfica
  par(mfrow = c(ceiling(length(componentes_activos)/2), 2), mar = c(4,4,2,1))
  
  # Graficar variogramas para cada componente
  for(i in seq_along(componentes_activos)) {
    idx <- componentes_activos[i]
    vgm_emp <- variogram(reformulate("1", paste0("score.", idx)), data_vario)
    
    if(nrow(vgm_emp) > 0) {  # Verificar datos
      plot(vgm_emp, fit_list[[i]], 
           main = paste("Componente", idx, "\nVar:", round(varianzas[idx], 2)),
           col = "darkred", pch = 19)
    } else {
      message("Variograma empírico vacío para componente ", idx)
    }
  }
} else {
  warning("No hay componentes activos para graficar")
}
# ===============================
# 6. Predicción y Incertidumbre (Con covarianza)
# ===============================
# Función para predecir scores en una fecha objetivo
predict_scores <- function(target_date) {
  # Si target_date es de clase Date, se convierte a años decimales.
  # De lo contrario, se asume que target_date ya está en la escala correcta.
  if (inherits(target_date, "Date")) {
    target_year <- as.numeric(format(target_date, "%Y")) + 
      as.numeric(format(target_date, "%j"))/365
  } else {
    target_year <- target_date
  }
  
  newdata <- data.frame(
    time = target_year,
    y = 0
  )
  coordinates(newdata) <- ~ time + y
  
  # Permitir extrapolación con maxdist = Inf
  scores_pred <- sapply(1:n_components, function(i) {
    krige(as.formula(paste0("score.", i, " ~ 1")),
          data_vario, newdata, fit_list[[i]], maxdist = Inf)$var1.pred
  })
  
  list(pred = scores_pred, covariance = cov(fpca$scores))
}

# Función para reconstruir la curva con intervalos de confianza
reconstruct_curve <- function(scores_pred, cov_matrix) {
  harmonics <- fpca$harmonics
  mean_curve <- eval.fd(X_common, mean.fd(fd_obj))
  
  # Reconstrucción de la curva media
  reconstructed <- mean_curve + rowSums(sapply(1:n_components, function(i) {
    scores_pred[i] * eval.fd(X_common, harmonics[i])
  }))
  
  # Cálculo de varianza en cada punto
  variance <- sapply(1:length(X_common), function(j) {
    phi <- sapply(1:n_components, function(i) eval.fd(X_common[j], harmonics[i]))
    as.numeric(t(phi) %*% cov_matrix %*% phi)
  })
  
  list(
    mean = reconstructed,
    lower = reconstructed - 1.96 * sqrt(variance),
    upper = reconstructed + 1.96 * sqrt(variance)
  )
}

# ===============================
# 7. Validación y Backtesting 
# ===============================
# Función para validación cruzada temporal
temporal_cv <- function() {
  rmse <- numeric(nrow(data_vario))
  coverage <- numeric(nrow(data_vario))
  
  for(i in 1:nrow(data_vario)) {
    train_data <- data_vario[-i,]
    # Aquí se fuerza la fecha de prueba a ser 1 año más que el máximo observado
    test_date <- max(data_vario$time) + 1  
    
    pred <- predict_scores(test_date)
    rec <- reconstruct_curve(pred$pred, pred$covariance)
    
    observed <- Y_matrix[, i]
    rmse[i] <- sqrt(mean((rec$mean - observed)^2))
    coverage[i] <- mean(observed >= rec$lower & observed <= rec$upper)
  }
  
  list(rmse = mean(rmse), coverage = mean(coverage))
}


# Ejecución de validación cruzada
cv_results <- temporal_cv()
message("\nResultados Validación Cruzada:")
message("RMSE Promedio: ", round(cv_results$rmse, 2))
message("Cobertura IC 95%: ", round(cv_results$coverage * 100, 1), "%")

# ===============================
# 8.Predicción (Gráfico 2D)
# ===============================
# Predicción para una fecha específica
target_date <- as.Date("2025-01-01")
prediction <- predict_scores(target_date)
recon <- reconstruct_curve(prediction$pred, prediction$covariance)

# Preparación de datos para gráfico
df_pred <- data.frame(
  X = X_common,
  Mean = recon$mean,
  Lower = recon$lower,
  Upper = recon$upper
)

colnames(df_pred) <- c("X", "Mean", "Lower", "Upper")

# Visualización interactiva con plotly
plot_ly(data = df_pred, x = ~X) %>%
  add_ribbons(ymin = ~Lower, ymax = ~Upper,
              fillcolor = "rgba(0,100,80,0.2)",
              line = list(color = "rgba(255,255,255,0)"),
              name = "IC 95%") %>%
  add_lines(y = ~Mean, line = list(color = "blue"), name = "Predicción") %>%
  layout(title = paste("Predicción para", target_date),
         xaxis = list(title = "X Relativa"),
         yaxis = list(title = "Y Relativa"))

# ===========================================
# 9.EXPORTACIÓN DE LÍNEA COSTERA
# ===========================================

# 1. Extraer y preparar los valores de predicción CORRECTAMENTE
y_pred <- as.vector(recon$mean)  # Convertir matriz a vector

# Verificación crítica
message("\n=== VERIFICACIÓN DE DATOS ===")
message("Longitud X_common: ", length(X_common))
message("Longitud y_pred: ", length(y_pred))
message("Rango X original: ", paste(range(gdf_df$X), collapse = " - "))
message("Rango Y original: ", paste(range(gdf_df$Y), collapse = " - "))
message("Rango Y predicho: ", paste(range(y_pred), collapse = " - "))

# 2. Crear data frame con estructura garantizada
coast_pred <- data.frame(
  coord_X = X_common,  # Valores X en escala original
  coord_Y = y_pred     # Valores Y en escala original
)

# 3. Verificación visual inmediata
plot(coast_pred$coord_X, coast_pred$coord_Y, type = 'l', col = 'red', lwd = 2,
     main = "Línea Costera Predicha",
     xlab = "Coordenada X", ylab = "Coordenada Y")
points(gdf_df$X, gdf_df$Y, pch = 20, cex = 0.5, col = 'blue')
legend("topright", legend = c("Predicción", "Datos originales"),
       col = c("red", "blue"), lty = c(1, NA), pch = c(NA, 20))

# 4. Creación del LINESTRING con verificación
coords <- as.matrix(coast_pred[, c("coord_X", "coord_Y")])
if(any(is.na(coords))) stop("Existen valores NA en las coordenadas")
coast_line <- st_linestring(coords)

# 5. Crear objeto sf con metadatos completos
coast_sf <- st_sf(
  geometry = st_sfc(coast_line, crs = st_crs(gdf)),
  fecha_prediccion = as.character(target_date),
  modelo = "FPCA-Kriging",
  resolucion = "100 puntos",
  fuente = "Datos Univalle",
  stringsAsFactors = FALSE
)

# 6. Guardar el archivo GeoJSON
output_final <- "G:/Univalle/linea_costera_2025_predicha.geojson"
st_write(coast_sf, output_final, delete_dsn = TRUE)

# 7. Mensaje de confirmación detallado
message("\n=== ARCHIVO GUARDADO CORRECTAMENTE ===")
message("Ubicación: ", output_final)
message("Sistema de referencia: ", st_crs(coast_sf)$input)
message("Número de puntos: ", nrow(coast_pred))
message("Extensión geográfica:")
print(st_bbox(coast_sf))
# ==========================================
# VALIDACIÓN 1: Curvas funcionales
# ==========================================
plot(fd_obj, col = viridis::viridis(ncol(fd_obj$coefs)), lwd = 2,
     main = "Curvas Funcionales Estimadas",
     xlab = "X común", ylab = "Y")
grid()

# ==========================================
# VALIDACIÓN 2: Scores FPCA
# ==========================================
# Matriz de scores
print("Matriz de Scores FPCA (primeros 5):")
print(head(fpca$scores))

# Boxplots de scores
scores_df <- as.data.frame(fpca$scores)
scores_df$fecha <- unique(gdf_df$DATE)[1:nrow(scores_df)]

scores_long <- scores_df %>%
  tidyr::pivot_longer(cols = starts_with("score."), names_to = "componente", values_to = "valor")

ggplot(scores_long, aes(x = componente, y = valor)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "Distribución de Scores por Componente",
       y = "Score", x = "Componente")

# ==========================================
# VALIDACIÓN 3: Varianza explicada FPCA
# ==========================================
var_exp <- fpca$values / sum(fpca$values)
cumvar_exp <- cumsum(var_exp)

barplot(var_exp, names.arg = paste0("PC", 1:length(var_exp)),
        main = "Varianza Explicada por Componente",
        ylab = "Proporción", col = "steelblue")
lines(cumvar_exp, type = "b", pch = 19, col = "darkred")
abline(h = 0.95, col = "gray", lty = 2)
legend("topright", legend = c("Acumulado"), col = c("darkred"), lty = 1)

# ==========================================
# VALIDACIÓN 4: Semivariogramas
# ==========================================
for(i in componentes_activos) {
  componente <- paste0("score.", i)
  cat("\nComponente:", componente, "\n")
  
  form <- reformulate("1", response = componente)
  vgm_emp <- variogram(form, data = data_vario)
  vgm_fit <- fit_list[[which(componentes_activos == i)]]
  
  print(vgm_fit)
  
  plot(vgm_emp, vgm_fit, main = paste("Semivariograma para", componente))
}

# Función para graficar el variograma empírico y el modelo ajustado
plot_variograma <- function(vgm_emp, ajuste, componente) {
  plot(vgm_emp, main = paste("Semivariograma para", componente),
       col = "darkred", pch = 19)
  lines(variogramLine(ajuste, maxdist = max(vgm_emp$dist)), col = "blue", lwd = 2)
}

# En el bucle de validación:
for(i in componentes_activos) {
  componente <- paste0("score.", i)
  form <- reformulate("1", response = componente)
  vgm_emp <- variogram(form, data = data_vario)
  ajuste <- fit_list[[which(componentes_activos == i)]]
  
  cat("\nComponente:", componente, "\n")
  print(ajuste)
  
  plot_variograma(vgm_emp, ajuste, componente)
}

# Después de ajustar los variogramas, verifica los rangos
sapply(fit_list, function(m) m$range[2]) %>% 
  as.data.frame() %>% 
  setNames("Rango") %>% 
  mutate(Componente = componentes_activos) %>% 
  mutate(Rango_Años = Rango/365) %>% 
  knitr::kable(caption = "Rangos de Variograma en Años")



# ==========================================
# VALIDACIÓN Análisis de Kriging 
# ==========================================
# Cargar librerías necesarias
library(gstat)
library(sp)
library(ggplot2)
library(dplyr)

# --- 1. Función de Validación Cruzada Mejorada ---
cv_analysis_corrected <- function(componente, sp_data, vario_model) {
  # Crear fórmula a partir del nombre del componente
  formula <- as.formula(paste0(componente, " ~ 1"))
  
  # Verificar que sp_data sea un SpatialPointsDataFrame
  if(!inherits(sp_data, "SpatialPointsDataFrame")) {
    stop("Los datos deben ser de clase SpatialPointsDataFrame")
  }
  
  # Realizar validación cruzada con manejo de errores
  cv_result <- tryCatch({
    krige.cv(formula, locations = sp_data, model = vario_model, 
             nfold = nrow(sp_data), verbose = FALSE)
  }, error = function(e) {
    message("Error en CV para ", componente, ": ", e$message)
    return(NULL)
  })
  
  if(is.null(cv_result)) return(NULL)
  
  # Calcular métricas de error (Error medio, RMSE y R²)
  residuals <- cv_result$residual
  observed <- cv_result$observed
  mean_observed <- mean(observed, na.rm = TRUE)
  
  metrics <- data.frame(
    Componente = componente,
    Error_Medio = mean(residuals, na.rm = TRUE),
    RMSE = sqrt(mean(residuals^2, na.rm = TRUE)),
    R2 = ifelse(var(observed, na.rm = TRUE) > 0,
                1 - sum(residuals^2, na.rm = TRUE) / 
                  sum((observed - mean_observed)^2, na.rm = TRUE),
                NA)
  )
  
  list(
    cv_result = cv_result,
    metrics = metrics,
    residuals = residuals,
    observed = observed
  )
}

# --- 2. Aplicar Validación Cruzada a Todos los Componentes Activos ---
# Se asume que 'componentes_activos' es un vector de números (o strings) y que en los datos la
# variable se llama "score.X" para cada componente
cv_results_corrected <- lapply(paste0("score.", componentes_activos), function(comp) {
  idx <- which(componentes_activos == as.numeric(gsub("score.", "", comp)))
  cv_analysis_corrected(comp, data_vario, fit_list[[idx]])
})
# Filtrar resultados válidos (no nulos)
valid_cv <- cv_results_corrected[!sapply(cv_results_corrected, is.null)]
cv_metrics <- if (length(valid_cv) > 0) do.call(rbind, lapply(valid_cv, function(x) x$metrics)) else NULL

if(!is.null(cv_metrics)) {
  print("Métricas de Validación Cruzada (Corregida):")
  print(cv_metrics)
} else {
  warning("No se pudieron calcular métricas de validación cruzada")
}

# --- 3. Análisis de Residuos (Gráficos) ---
if(length(valid_cv) > 0) {
  # Configurar área de gráficos (histogramas)
  op <- par(mfrow = c(ceiling(length(valid_cv) / 2), 2), mar = c(4,4,2,1))
  for(i in seq_along(valid_cv)) {
    hist(valid_cv[[i]]$residuals, breaks = 20, col = "skyblue",
         main = paste("Residuos", valid_cv[[i]]$metrics$Componente),
         xlab = "Residuos", ylab = "Frecuencia")
  }
  par(op)  # Restaurar parámetros
  
  # Configurar área para Q-Q Plots
  op2 <- par(mfrow = c(ceiling(length(valid_cv) / 2), 2), mar = c(4,4,2,1))
  for(i in seq_along(valid_cv)) {
    qqnorm(valid_cv[[i]]$residuals, main = paste("Q-Q Plot", valid_cv[[i]]$metrics$Componente))
    qqline(valid_cv[[i]]$residuals, col = "red")
  }
  par(op2)
} else {
  warning("No hay datos válidos para análisis de residuos")
}

# --- 4. Función para Comparar Modelos de Variograma con Robustez ---
compare_vario_models_robust <- function(componente, sp_data) {
  form <- as.formula(paste0(componente, " ~ 1"))
  
  # Asegurar que 'time' sea numérica (si no, convertirla)
  if(!is.numeric(sp_data$time)) {
    sp_data$time <- as.numeric(as.character(sp_data$time))
  }
  
  # Verificar que la variable 'time' tenga un rango válido
  time_range <- diff(range(sp_data$time, na.rm = TRUE))
  if(time_range <= 0) {
    message("El rango de 'time' no es válido para ", componente)
    return(NULL)
  }
  
  vg_emp <- tryCatch({
    variogram(form, data = sp_data, 
              cutoff = time_range * 0.8,
              width = time_range / 15)
  }, error = function(e) {
    message("Error calculando variograma empírico para ", componente, ": ", e$message)
    return(NULL)
  })
  
  if(is.null(vg_emp) || nrow(vg_emp) < 3) {
    message("Variograma empírico no válido para ", componente)
    return(NULL)
  }
  
  # Parámetros iniciales para los modelos candidatos
  var_total <- var(sp_data[[componente]], na.rm = TRUE)
  rango <- max(time_range * 0.3, .Machine$double.eps)
  
  models <- list(
    Exponencial = vgm(psill = 0.7 * var_total, model = "Exp", range = rango, nugget = 0.3 * var_total),
    Esférico   = vgm(psill = 0.7 * var_total, model = "Sph", range = rango, nugget = 0.3 * var_total),
    Gaussiano  = vgm(psill = 0.7 * var_total, model = "Gau", range = rango * 0.5, nugget = 0.3 * var_total),
    Matern     = vgm(psill = 0.7 * var_total, model = "Mat", range = rango, nugget = 0.3 * var_total, kappa = 1.5)
  )
  
  # Ajuste de cada uno de los modelos con manejo de errores
  fits <- lapply(models, function(m) {
    tryCatch(
      fit.variogram(vg_emp, model = m, fit.method = 7),
      error = function(e) {
        message("Error ajustando modelo ", m$model[2], " para ", componente, ": ", e$message)
        return(NULL)
      }
    )
  })
  
  # Graficar la comparación entre el variograma empírico y las líneas de los modelos ajustados
  plot(vg_emp, main = paste("Comparación de Modelos -", componente), 
       col = "darkgray", pch = 19, xlab = "Distancia", ylab = "Semivarianza")
  cols <- c("blue", "green", "red", "purple")
  for(j in seq_along(fits)) {
    if(!is.null(fits[[j]])) {
      lines(variogramLine(fits[[j]], maxdist = max(vg_emp$dist)), 
            col = cols[j], lwd = 2)
    }
  }
  legend("bottomright", legend = names(fits), col = cols, lwd = 2, cex = 0.8)
  
  # Calcular el SSE (suma de errores cuadrados) para cada modelo
  sse <- sapply(seq_along(fits), function(j) {
    fit <- fits[[j]]
    if(!is.null(fit)) {
      pred <- variogramLine(fit, maxdist = max(vg_emp$dist), n = nrow(vg_emp))$gamma
      sum((pred - vg_emp$gamma)^2)
    } else {
      NA
    }
  })
  
  # Retornar los resultados comparativos
  data.frame(
    Componente = componente,
    Modelo = names(fits),
    SSE = sse,
    stringsAsFactors = FALSE
  )
}

# --- 5. Aplicar Comparación de Modelos a los Componentes Activos ---
model_comparisons_list <- lapply(paste0("score.", componentes_activos), function(comp) {
  tryCatch(
    compare_vario_models_robust(comp, data_vario),
    error = function(e) {
      message("Error en la comparación para ", comp, ": ", e$message)
      return(NULL)
    }
  )
})
# Solo se combinan aquellos elementos no nulos
model_comparisons_list <- model_comparisons_list[!sapply(model_comparisons_list, is.null)]
model_comparisons <- if(length(model_comparisons_list) > 0) do.call(rbind, model_comparisons_list) else NULL

if(is.null(model_comparisons) || nrow(model_comparisons) == 0) {
  warning("No se generaron comparaciones de modelos válidas.")
} else {
  print("Comparación de Modelos:")
  print(model_comparisons)
  # Guardar resultados
  write.csv(model_comparisons, "comparacion_variogramas.csv", row.names = FALSE)
}

# --- 6. Gráfico de Comparación de Modelos con ggplot2 ---
if(!is.null(model_comparisons) && nrow(model_comparisons) > 0) {
  ggplot(model_comparisons, aes(x = reorder(Modelo, SSE), y = log10(SSE), fill = Modelo)) +
    geom_col() +
    facet_wrap(~ Componente, nrow = 1) +
    coord_flip() +
    labs(title = "Comparación de Modelos por Componente",
         subtitle = "Escala logarítmica de SSE (menor es mejor)",
         y = "log10(Suma de Errores Cuadrados)",
         x = "Modelo") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom")
} else {
  warning("No hay datos para el gráfico de comparación de modelos.")
}

# --- 7. Selección de los Mejores Modelos ---
if(!is.null(model_comparisons) && nrow(model_comparisons) > 0) {
  best_models <- model_comparisons %>%
    group_by(Componente) %>%
    filter(!is.na(SSE)) %>%
    slice_min(SSE)
  print("Mejores modelos para cada componente:")
  print(best_models)
} else {
  best_models <- NULL
  warning("No se pudieron identificar los mejores modelos.")
}

# --- 8. Mapa de Incertidumbre (Predicción y Varianza) ---
if(length(componentes_activos) >= 1) {
  message("\nGenerando gráficos de predicción y varianza con ggplot2...")
  
  # Crear grilla segura usando la variable 'time'
  time_seq <- seq(min(data_vario$time, na.rm = TRUE), max(data_vario$time, na.rm = TRUE), length.out = 100)
  grilla_df <- data.frame(time = time_seq, y = 0)
  coordinates(grilla_df) <- ~ time + y
  
  # Función para generar gráficos usando ggplot2
  create_ggplot <- function(data, y_var, title) {
    ggplot(data, aes(x = time, y = .data[[y_var]])) +
      geom_line(color = ifelse(y_var == "pred", "blue", "red"), linewidth = 1) +
      labs(title = title, x = "Tiempo", y = y_var) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  for(i in 1:min(4, length(componentes_activos))) {
    comp <- paste0("score.", componentes_activos[i])
    message("Procesando: ", comp)
    
    pred_data <- tryCatch({
      krige(as.formula(paste0(comp, " ~ 1")),
            locations = data_vario,
            newdata = grilla_df,
            model = fit_list[[i]])
    }, error = function(e) {
      message("Error en kriging para ", comp, ": ", e$message)
      return(NULL)
    })
    
    if(!is.null(pred_data)) {
      # Convertir predicciones a data.frame
      plot_data <- data.frame(
        time = coordinates(grilla_df)[,1],
        pred = pred_data$var1.pred,
        var = pred_data$var1.var
      )
      
      # Graficar predicción
      print(create_ggplot(plot_data, "pred", paste("Predicción -", comp)))
      # Graficar varianza
      print(create_ggplot(plot_data, "var", paste("Varianza -", comp)))
    }
  }
}

# --- 9. Ajuste de Modelos Anidados y Simples para 'score.1' ---
if(exists("data_vario") && "score.1" %in% names(data_vario)) {
  vgm_emp <- variogram(score.1 ~ 1, data_vario)
  
  # Definir un modelo anidado para score.1
  modelo_anidado_score1 <- vgm(psill = 1e10, model = "Exp", range = 5000,
                               nugget = 5e9, add.to = vgm(psill = 1e9, model = "Sph", range = 1000))
  fit_anidado <- tryCatch({
    fit.variogram(vgm_emp, modelo_anidado_score1)
  }, warning = function(w) {
    message("Aviso en ajuste de modelo anidado: ", w$message)
    fit.variogram(vgm_emp, modelo_anidado_score1)
  }, error = function(e) {
    message("Error en ajuste de modelo anidado: ", e$message)
    return(NULL)
  })
  
  if(!is.null(fit_anidado)) {
    plot(vgm_emp, fit_anidado, main = "Modelo Anidado para score.1")
    if(length(fit_list) >= 1 && !is.null(fit_list[[1]])) {
      plot(vgm_emp, fit_list[[1]], main = "Modelo Simple para score.1")
    }
  }
  
  # Ajustes de modelos adicionales para 'score.1' y 'score.4'
  modelo_score1 <- vgm(psill = 0.7 * var(data_vario$score.1, na.rm = TRUE), model = "Ste", 
                       range = diff(range(data_vario$time, na.rm = TRUE)) / 2,
                       nugget = 0.3 * var(data_vario$score.1, na.rm = TRUE), kappa = 0.5)
  modelo_score4 <- vgm(psill = 0.8 * var(data_vario$score.4, na.rm = TRUE), model = "Sph",
                       range = diff(range(data_vario$time, na.rm = TRUE)) * 0.8,
                       nugget = 0.2 * var(data_vario$score.4, na.rm = TRUE))
  
  # Graficar variogramas empíricos vs ajustados para los componentes score.1 a score.4
  par(mfrow = c(2, 2))
  for(i in 1:4) {
    comp <- paste0("score.", i)
    vgm_emp_comp <- variogram(as.formula(paste0(comp, " ~ 1")), data_vario)
    if(i <= length(fit_list) && !is.null(fit_list[[i]])) {
      plot(vgm_emp_comp, fit_list[[i]], 
           main = paste("Ajuste para", comp),
           col = "darkred", pch = 19,
           xlab = "Distancia", ylab = "Semivarianza")
      grid()
    }
  }
}

# --- 10. Validación Cruzada Adicional Usando los Mejores Modelos ---
if(!is.null(best_models) && nrow(best_models) > 0){
  cv_results <- lapply(1:4, function(i) {
    comp <- paste0("score.", i)
    best_model <- best_models$Modelo[best_models$Componente == comp]
    
    # Mapeo flexible de nombres para comparar con fit_list
    model_map <- c("Esférico" = "Sph", "Exponencial" = "Exp", 
                   "Gaussiano" = "Gau", "Matern" = "Mat")
    
    tryCatch({
      model_index <- which(sapply(fit_list, function(x) x$model[2]) == model_map[best_model])
      if(length(model_index) == 0) stop("Modelo ", best_model, " no encontrado en fit_list")
      krige.cv(formula = as.formula(paste0(comp, " ~ 1")),
               locations = data_vario,
               model = fit_list[[model_index[1]]],
               nfold = 5,
               verbose = FALSE)
    }, error = function(e) {
      message("Error en CV para ", comp, ": ", e$message)
      return(NULL)
    })
  })
  
  if(all(!sapply(cv_results, is.null))) {
    cv_metrics_cv <- sapply(cv_results, function(x) {
      c(ME = mean(x$residual, na.rm = TRUE),
        RMSE = sqrt(mean(x$residual^2, na.rm = TRUE)),
        R2 = ifelse(var(x$observed, na.rm = TRUE) > 0,
                    1 - sum(x$residual^2, na.rm = TRUE) /
                      sum((x$observed - mean(x$observed, na.rm = TRUE))^2, na.rm = TRUE),
                    NA)
      )
    })
    colnames(cv_metrics_cv) <- paste0("score.", 1:4)
    print("Métricas de Validación Cruzada (segunda etapa):")
    print(round(cv_metrics_cv, 3))
  } else {
    warning("Algunos modelos no permitieron la validación cruzada en la segunda etapa.")
  }
} else {
  warning("No se han definido best_models para validación cruzada adicional.")
}

# --- 11. Resumen de Diagnóstico Completo ---
if(!is.null(cv_metrics)) {
  diagnostic_summary <- list(
    Componentes_Activos = componentes_activos,
    Varianza_Explicada = data.frame(
      Componente = paste0("score.", componentes_activos),
      Varianza = varianzas[componentes_activos],
      Proporcion = varianzas[componentes_activos] / sum(varianzas)
    ),
    Metricas_CV = cv_metrics,
    Comparacion_Modelos = if(!is.null(model_comparisons) && nrow(model_comparisons) > 0) model_comparisons else NULL
  )
  
  print("Resumen Completo de Diagnóstico:")
  print(diagnostic_summary)
  
  # Guardar resultados en archivo (verificar las rutas)
  write.csv(cv_metrics, "G:/Univalle/metricas_validacion_corregidas.csv", row.names = FALSE)
  saveRDS(diagnostic_summary, "G:/Univalle/diagnostic_summary_corregido.rds")
} else {
  warning("No se pudo generar el resumen de diagnóstico por falta de datos")
}

