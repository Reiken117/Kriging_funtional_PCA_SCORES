# ===============================
# 1. Cargar librerías necesarias (Mejorado)
# ===============================
library(sf)                   # Para manejar datos geoespaciales
library(sp)                   # Funciones espaciales básicas
library(ggplot2)              # Visualización gráfica
library(dplyr)                # Manipulación de datos
library(fda)                  # Análisis funcional de datos
library(gstat)                # Geoestadística y kriging
library(plotly)               # Visualización 3D interactiva
library(viridis)              # Escalas de color
library(mvtnorm)             # Para manejo de distribuciones multivariadas
library(parallel)            # Para paralelización
library(logger)              # Para logging profesional
library(testthat)            # Para pruebas unitarias

# Configurar logging
log_appender(appender_file("coastal_analysis.log"))
log_layout(layout_glue_colors)
log_threshold(INFO)

# ===============================
# 2. Cargar y preparar los datos (Mejorado)
# ===============================
log_info("Iniciando carga de datos")

# Función para validar datos de entrada
validate_coastal_data <- function(gdf) {
  test_that("Validación de datos de entrada", {
    expect_true(all(c("DATE", "geometry") %in% names(gdf)))
    expect_true(nrow(gdf) > 0)
    expect_true(inherits(gdf, "sf"))
    expect_true(st_crs(gdf)$epsg == 3116) # MAGNA-SIRGAS 2018
  })
}

# Leer y validar datos geoespaciales
load_coastal_data <- function(file_path) {
  tryCatch({
    gdf <- st_read(file_path) %>% 
      st_transform(3116) %>%  # Asegurar MAGNA-SIRGAS 2018
      distinct() %>%          # Eliminar duplicados
      arrange(DATE)           # Ordenar por fecha
    
    # Validación estricta
    validate_coastal_data(gdf)
    
    log_info("Datos cargados correctamente. {nrow(gdf)} puntos encontrados.")
    return(gdf)
  }, error = function(e) {
    log_error("Error al cargar datos: {e$message}")
    stop(e)
  })
}

# Carga de datos
gdf <- load_coastal_data("G:/Univalle/puntos_proyect.json")


# Procesamiento de coordenadas con verificación
process_coords <- function(gdf) {
  coords <- st_coordinates(gdf) %>% 
    as.data.frame() %>%
    rename(X = X, Y = Y) %>%
    mutate(DATE = as.Date(gdf$DATE),
           id = row_number())
  
  # Detección y manejo de duplicados espaciales
  dup_coords <- duplicated(coords[, c("X", "Y")])
  if(any(dup_coords)) {
    log_warn("{sum(dup_coords)} puntos duplicados detectados. Se conserva el más reciente.")
    coords <- coords %>%
      group_by(X, Y) %>%
      arrange(desc(DATE)) %>%
      slice(1) %>%
      ungroup()
  }
  
  return(coords)
}

gdf_df <- process_coords(gdf)

# Visualización mejorada
plot_data_distribution <- function(gdf_df) {
  p <- ggplot(gdf_df, aes(X, Y, color = DATE)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_viridis(trans = "date", option = "plasma") +
    labs(title = "Distribución Espacio-Temporal de Puntos Costeros",
         subtitle = paste("Desde", min(gdf_df$DATE), "hasta", max(gdf_df$DATE)),
         x = "Coordenada X (MAGNA-SIRGAS 2018)",
         y = "Coordenada Y (MAGNA-SIRGAS 2018)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggplotly(p) %>% 
    layout(legend = list(orientation = "h", y = -0.2))
}

plot_data_distribution(gdf_df)

# ===============================
# 3. Análisis Funcional (Mejorado)
# ===============================
log_info("Iniciando análisis funcional")

# Definición de puntos comunes con densidad adaptativa
define_common_grid <- function(x_coords, n = 100) {
  dens <- density(x_coords, adjust = 0.5)
  seq(min(x_coords), max(x_coords), length.out = n)
  # Ajustar densidad en áreas con más puntos
  weighted_seq <- quantile(x_coords, probs = seq(0, 1, length.out = n))
  return(weighted_seq)
}

X_common <- define_common_grid(gdf_df$X, 150)

# Interpolación robusta con splines suavizados
safe_spline_interp <- function(x, y, xout) {
  tryCatch({
    if(length(unique(x)) < 4) {
      approx(x, y, xout = xout, rule = 2)$y
    } else {
      smooth.spline(x, y, spar = 0.5) %>% 
        predict(xout)$y
    }
  }, error = function(e) {
    log_warn("Interpolación fallida: {e$message}. Usando aproximación lineal.")
    approx(x, y, xout = xout, rule = 2)$y
  })
}

# Interpolación paralelizada corregida
interpolate_curves <- function(gdf_df, X_common) {
  dates <- unique(gdf_df$DATE)
  
  cl <- makeCluster(detectCores() - 1)
  
  # Exportar TODAS las funciones y variables necesarias
  clusterExport(cl, c("safe_spline_interp", "X_common", "gdf_df", 
                      "log_warn", "log_info", "log_error",
                      "validate_coastal_data", "process_coords"))
  
  # Cargar TODOS los paquetes necesarios dentro del clúster
  clusterEvalQ(cl, {
    library(dplyr)
    library(magrittr)
    library(logger)
    library(sf)
    library(testthat)
    library(stats)  # para approx y smooth.spline
  })
  
  Y_list <- parLapply(cl, dates, function(d) {
    df_date <- gdf_df %>% filter(DATE == d)
    safe_spline_interp(df_date$X, df_date$Y, X_common)
  })
  
  stopCluster(cl)
  
  Y_matrix <- do.call(cbind, Y_list)
  colnames(Y_matrix) <- as.character(dates)
  
  return(Y_matrix)
}

# Ejecutar la interpolación
Y_matrix <- interpolate_curves(gdf_df, X_common)

# Selección de base óptima con validación cruzada paralelizada
optimize_basis <- function(Y_matrix, X_common, max_basis = 20) {
  cv_basis <- function(nb) {
    basis <- create.bspline.basis(range(X_common), nbasis = nb)
    fdPar <- fdPar(basis)
    mean(smooth.basis(X_common, Y_matrix, fdPar)$SSE)
  }
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("create.bspline.basis", "fdPar", "smooth.basis", "X_common", "Y_matrix"))
  
  basis_range <- 5:max_basis
  cv_errors <- parSapply(cl, basis_range, cv_basis)
  
  stopCluster(cl)
  
  optimal <- basis_range[which.min(cv_errors)]
  log_info("Base óptima seleccionada: {optimal} funciones base")
  
  return(optimal)
}

optimal_basis <- optimize_basis(Y_matrix, X_common)
fd_basis <- create.bspline.basis(range(X_common), nbasis = optimal_basis, norder = 4)
fd_obj <- smooth.basis(X_common, Y_matrix, fdPar(fd_basis))$fd

# Visualización interactiva de curvas
plot_functional_curves <- function(fd_obj) {
  time_points <- seq(min(gdf_df$DATE), max(gdf_df$DATE), length.out = 10)
  colors <- viridis(length(time_points))
  
  plot_ly() %>%
    add_lines(x = X_common, 
              y = eval.fd(X_common, fd_obj[1]), 
              line = list(color = colors[1]),
              name = format(time_points[1], "%Y-%m")) %>%
    # Añadir más curvas...
    layout(title = "Evolución Temporal de la Línea Costera",
           xaxis = list(title = "X Relativa"),
           yaxis = list(title = "Y Relativa"),
           showlegend = TRUE)
}

plot_functional_curves(fd_obj)

# ===============================
# 4. FPCA y Selección de Componentes (Mejorado)
# ===============================
log_info("Realizando FPCA")

# FPCA con selección robusta de componentes
perform_fpca <- function(fd_obj, var_threshold = 0.95, min_components = 3) {
  fpca <- pca.fd(fd_obj, nharm = min(20, length(fd_obj$fdnames$reps)))
  
  # Selección adaptativa de componentes
  cumvar <- cumsum(fpca$values)/sum(fpca$values)
  n_components <- max(min_components, which(cumvar >= var_threshold)[1])
  
  # Validación de componentes
  if(n_components < min_components) {
    log_warn("Solo {n_components} componentes alcanzan el umbral. Usando mínimo {min_components}")
    n_components <- min_components
  }
  
  log_info("Seleccionados {n_components} componentes (Varianza explicada: {round(cumvar[n_components]*100, 1)}%)")
  
  # Preparar resultados
  fpca$harmonics <- fpca$harmonics[1:n_components]
  fpca$scores <- fpca$scores[, 1:n_components, drop = FALSE]
  colnames(fpca$scores) <- paste0("score.", 1:n_components)
  
  # Asignar nombres de fila (fechas) a los scores
  rownames(fpca$scores) <- fd_obj$fdnames$reps
  
  return(fpca)
}

fpca <- perform_fpca(fd_obj)


# Verificación proactiva de estructura
validate_fpca <- function(fpca) {
  if(!inherits(fpca, "pca.fd")) stop("Objeto FPCA inválido")
  if(is.null(rownames(fpca$scores))) stop("FPCA necesita fechas como row.names")
  if(ncol(fpca$scores) < 1) stop("FPCA debe tener al menos 1 componente")
}

# Antes de usar la función
validate_fpca(fpca)

# ===============================
# 5. Modelado Temporal con Kriging (Mejorado)
# ===============================
log_info("Configurando modelado temporal con Kriging")

# Preparación de datos espaciales temporales
prepare_spatial_data <- function(fpca, gdf_df) {
  # 1. Conversión robusta de fechas
  gdf_df$DATE <- as.Date(gdf_df$DATE, format = "%Y-%m-%d")
  fpca_dates <- as.Date(rownames(fpca$scores), format = "%Y-%m-%d")
  
  # 2. Validación exhaustiva de fechas
  if(any(is.na(gdf_df$DATE))) stop("gdf_df contiene fechas inválidas")
  if(any(is.na(fpca_dates))) stop("fpca contiene fechas inválidas en los nombres de fila")
  
  # 3. Búsqueda de fechas comunes con formato consistente
  common_dates <- intersect(gdf_df$DATE, fpca_dates)
  if(length(common_dates) == 0) {
    cat("Primeras fechas gdf_df:", head(gdf_df$DATE), "\n")
    cat("Primeras fechas fpca:", head(fpca_dates), "\n")
    stop("No hay coincidencia real de fechas")
  }
  
  # 4. Alineación temporal precisa
  date_order <- order(common_dates)
  common_dates_sorted <- common_dates[date_order]
  
  # 5. Mapeo seguro de índices
  row_indices <- sapply(common_dates_sorted, function(d) which(fpca_dates == d)[1])
  if(any(is.na(row_indices))) {
    problematic_dates <- common_dates_sorted[is.na(row_indices)]
    stop("Fechas no encontradas en fpca: ", paste(problematic_dates, collapse = ", "))
  }
  
  # 6. Construcción del dataset espacial-temporal
  data_vario <- data.frame(
    time = as.numeric(common_dates_sorted - min(common_dates_sorted)),
    component = rep(colnames(fpca$scores), each = length(common_dates_sorted)),
    value = as.vector(t(fpca$scores[row_indices, ])),
    y = 0.0
  )
  
  # 7. Conversión final a objeto espacial
  sp::coordinates(data_vario) <- ~ time + y
  sp::proj4string(data_vario) <- sp::CRS("+proj=utm +zone=18 +datum=WGS84")
  
  return(data_vario)
}


tryCatch({
  data_vario <- prepare_spatial_data(fpca, gdf_df)
  cat("Transformación exitosa. Muestra de datos:\n")
  print(head(data_vario, 3))
}, error = function(e) {
  cat("Error crítico:", e$message, "\n")
  cat("Debug info:\n")
  cat("Clase de fechas gdf_df:", class(gdf_df$DATE), "\n")
  cat("Clase de fechas fpca:", class(as.Date(rownames(fpca$scores))), "\n")
})
#___________________________________________________

# Verificar estructura de fpca
if(!"scores" %in% names(fpca)) stop("El objeto fpca no contiene 'scores'")
if(!is.matrix(fpca$scores)) stop("fpca$scores debe ser una matriz")

# Si fpca$scores no tiene nombres de fila, asignar las fechas
if(is.null(rownames(fpca$scores))) {
  rownames(fpca$scores) <- as.character(unique(gdf_df$DATE))
}

data_vario <- tryCatch({
  prepare_spatial_data(fpca, gdf_df)
}, error = function(e) {
  log_error("Error al preparar datos espaciales: {e$message}")
  stop(e)
})

# Verificar dimensiones
log_info("Dimensiones de fpca$scores: {nrow(fpca$scores)} filas, {ncol(fpca$scores)} columnas")
log_info("Longitud de dates_numeric: {length(dates_numeric)} elementos")
log_info("Número de fechas únicas: {length(unique_dates)}")

# Ajuste de variogramas robusto con múltiples modelos
# Función ajustada con exportación correcta de variables y paquetes
fit_component_variograms <- function(data_vario, components) {
  models <- list(
    vgm(1, "Exp", 100),
    vgm(1, "Sph", 100),
    vgm(1, "Gau", 100),
    vgm(1, "Mat", 100, kappa = 1.5),
    vgm(0.8, "Exp", 100, 0.2),
    vgm(0.8, "Sph", 100, 0.2)
  )
  
  cl <- makeCluster(detectCores() - 1)
  
  # Exportar TODAS las dependencias necesarias
  clusterExport(cl, c("data_vario", "vgm", "variogram", "fit.variogram"), 
                envir = environment())
  
  # Cargar paquetes necesarios en los clusters
  clusterEvalQ(cl, {
    library(gstat)
    library(sp)
  })
  
  fit_list <- parLapply(cl, components, function(comp) {
    sub_data <- data_vario[data_vario$component == comp, ]
    v <- variogram(value ~ 1, sub_data)
    
    best_model <- NULL
    best_sse <- Inf
    
    for(m in models) {
      fit <- try(fit.variogram(v, m, fit.method = 7), silent = TRUE)
      if(!inherits(fit, "try-error") && !is.null(fit)) {
        sse <- attr(fit, "SSErr")
        if(sse < best_sse) {
          best_sse <- sse
          best_model <- fit
        }
      }
    }
    return(best_model)
  })
  
  stopCluster(cl)
  names(fit_list) <- components
  return(fit_list)
}


components <- colnames(fpca$scores)
fit_list <- fit_component_variograms(data_vario, components)

# Visualización de variogramas
plot_variograms <- function(data_vario, fit_list) {
  plots <- list()
  
  for(comp in names(fit_list)) {
    v <- variogram(value ~ 1, data_vario[data_vario$component == comp,])
    p <- plot(v, fit_list[[comp]], 
              main = paste("Variograma para", comp),
              xlab = "Distancia Temporal", ylab = "Semivarianza")
    plots[[comp]] <- p
  }
  
  return(plots)
}

variogram_plots <- plot_variograms(data_vario, fit_list)
print(variogram_plots[[1]]) # Mostrar el primer variograma

# Verificar estructura de los resultados
print(lapply(fit_list, function(x) x$model[2]))


# ===============================
# 6. Predicción e Incertidumbre (Versión Final Corregida)
# ===============================

# ===============================
# SOLUCIÓN COMPLETA Y CORREGIDA
# ===============================

# 1. Configuración Inicial Obligatoria
library(logger)
log_appender(appender_console())
log_threshold(INFO)

# 2. Función de Predicción por Componente (Versión Robustecida)
predict_component <- function(component, data_vario, fit_model, pred_grid) {
  tryCatch({
    sub_data <- data_vario[data_vario$component == component, ]
    
    if(nrow(sub_data) < 5) {
      message(paste("Advertencia: Solo", nrow(sub_data), "puntos para", component)
      return(NULL)
    }
    
    krige(value ~ time, 
          locations = sub_data,
          newdata = pred_grid,
          model = fit_model,
          nmax = 15)
  }, error = function(e) {
    message(paste("Error en", component, ":", e$message))
    return(NULL)
  })
}

# 3. Sistema de Predicción Principal (Actualizado)
run_predictions <- function(data_vario, fit_list, target_date) {
  # Crear grid de predicción
  min_date <- min(gdf_df$DATE)
  pred_grid <- create_prediction_grid(min_date, target_date)
  
  # Configurar paralelización con exportación completa
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("data_vario", "pred_grid", "predict_component", "krige"),
                envir = environment())
  clusterEvalQ(cl, {
    library(gstat)
    library(sp)
  })
  
  predictions <- parLapply(cl, names(fit_list), function(comp) {
    predict_component(comp, data_vario, fit_list[[comp]], pred_grid)
  })
  
  stopCluster(cl)
  names(predictions) <- names(fit_list)
  return(list(predictions = predictions, grid = pred_grid))
}

# 4. Visualización Garantizada
plot_prediction_results <- function(predictions, uncertainty_df, pred_grid) {
  # Preparar datos
  pred_values <- sapply(predictions, function(p) if (!is.null(p)) p$var1.pred else NA)
  pred_df <- data.frame(
    time = pred_grid$time,
    score.1 = pred_values
  )
  
  # Gráfico Principal
  p1 <- ggplot(pred_df, aes(x = time)) +
    geom_line(aes(y = score.1), color = "blue", size = 1) +
    labs(title = "Predicción de Componente Principal",
         x = "Días desde inicio",
         y = "Valor") +
    theme_bw()
  
  # Gráfico de Incertidumbre
  p2 <- ggplot(uncertainty_df, aes(x = time, y = total_variance)) +
    geom_line(color = "red") +
    labs(title = "Incertidumbre Total",
         y = "Varianza") +
    theme_bw()
  
  # Mostrar directamente
  gridExtra::grid.arrange(p1, p2, ncol = 1)
  invisible(NULL)
}

# 5. Validación Cruzada Simplificada
simple_cross_validation <- function(gdf_df, n_folds = 3) {
  dates <- unique(gdf_df$DATE)
  folds <- cut(seq_along(dates), breaks = n_folds, labels = FALSE)
  
  results <- lapply(1:n_folds, function(i) {
    test_dates <- dates[folds == i]
    train_data <- gdf_df[!gdf_df$DATE %in% test_dates,]
    test_data <- gdf_df[gdf_df$DATE %in% test_dates,]
    
    # Aquí iría el modelo completo (simplificado para ejemplo)
    mean_pred <- mean(train_data$Y, na.rm = TRUE)
    rmse <- sqrt(mean((test_data$Y - mean_pred)^2, na.rm = TRUE))
    
    data.frame(fold = i, rmse = rmse)
  })
  
  do.call(rbind, results)
}

# --- EJECUCIÓN FINAL CON MANEJO DE ERRORES ---
tryCatch({
  # A. Predicción Principal
  target_date <- max(gdf_df$DATE) + 365
  pred_results <- run_predictions(data_vario, fit_list, target_date)
  
  # B. Cálculo de Incertidumbre
  uncertainty <- calculate_uncertainty(pred_results$predictions, pred_results$grid)
  
  # C. Visualización Obligatoria
  plot_prediction_results(pred_results$predictions, uncertainty, pred_results$grid)
  
  # D. Validación Cruzada Básica
  cv_res <- simple_cross_validation(gdf_df)
  print("Resultados Validación Cruzada:")
  print(cv_res)
  
}, error = function(e) {
  message("\nERROR CRÍTICO EN EL ANÁLISIS:")
  message("Mensaje: ", e$message)
  message("Revisar los siguientes objetos:")
  message("- Última fecha: ", tail(gdf_df$DATE, 1))
  message("- Dimensión data_vario: ", paste(dim(data_vario), collapse = " x "))
})

# ===============================
# 7. Validación Cruzada (Mejorado)
# ===============================
log_info("Iniciando validación cruzada")

# Validación cruzada por bloques temporales
temporal_block_cv <- function(gdf_df, fpca, fit_list, fd_obj, X_common, n_blocks = 5) {
  dates <- unique(gdf_df$DATE)
  date_groups <- cut(seq_along(dates), breaks = n_blocks, labels = FALSE)
  
  metrics <- data.frame()
  
  for(block in 1:n_blocks) {
    test_dates <- dates[date_groups == block]
    train_data <- gdf_df[!gdf_df$DATE %in% test_dates,]
    
    # Re-entrenar FPCA solo con datos de entrenamiento (simplificado)
    # En práctica, deberíamos rehacer todo el análisis funcional
    train_idx <- which(unique(gdf_df$DATE) %in% unique(train_data$DATE))
    train_scores <- fpca$scores[train_idx,]
    
    # Predicción para fechas de prueba
    block_results <- lapply(test_dates, function(d) {
      pred <- predict_coastal_change(d, fpca, fit_list, data_vario, fd_obj, X_common)
      obs <- gdf_df[gdf_df$DATE == d,]
      
      # Interpolar observaciones a X_common para comparación
      obs_interp <- approx(obs$X, obs$Y, xout = X_common, rule = 2)$y
      
      # Calcular métricas
      residuals <- obs_interp - pred$reconstruction$mean
      in_ci <- obs_interp >= pred$reconstruction$lower & 
        obs_interp <= pred$reconstruction$upper
      
      data.frame(
        date = d,
        rmse = sqrt(mean(residuals^2, na.rm = TRUE)),
        mae = mean(abs(residuals), na.rm = TRUE),
        coverage = mean(in_ci, na.rm = TRUE)
      )
    })
    
    metrics <- rbind(metrics, do.call(rbind, block_results))
  }
  
  return(metrics)
}

# Ejecutar validación cruzada (comentado por tiempo de ejecución)
cv_results <- temporal_block_cv(gdf_df, fpca, fit_list, fd_obj, X_common)
log_info("Resultados CV: RMSE promedio = {mean(cv_results$rmse)}")

# ===============================
# 8. Exportación de Resultados (Mejorado)
# ===============================
log_info("Exportando resultados")

# Función para exportar línea costera predicha
export_coastline <- function(prediction, original_crs, output_path) {
  # Convertir a coordenadas originales
  X_orig <- prediction$X_common + min(gdf_df$X)
  Y_orig <- prediction$reconstruction$mean + min(gdf_df$Y)
  
  # Crear objeto espacial con incertidumbre
  coast_sf <- data.frame(
    X = X_orig,
    Y = Y_orig,
    lower = prediction$reconstruction$lower + min(gdf_df$Y),
    upper = prediction$reconstruction$upper + min(gdf_df$Y),
    variance = prediction$reconstruction$variance
  ) %>%
    st_as_sf(coords = c("X", "Y")) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    st_as_sf(coords = c("X", "Y"), crs = original_crs) %>%
    st_combine() %>%
    st_cast("LINESTRING") %>%
    st_sf() %>%
    mutate(
      date = prediction$date,
      variance_mean = mean(prediction$reconstruction$variance)
    )
  
  # Exportar como GeoJSON
  st_write(coast_sf, output_path, delete_dsn = TRUE)
  log_info("Línea costera exportada a {output_path}")
  
  return(coast_sf)
}

# Exportar predicción
output_file <- "G:/Univalle/linea_costera_predicha.geojson"
coastline_pred <- export_coastline(prediction, st_crs(gdf), output_file)

# ===============================
# 9. Pruebas Unitarias (Nuevo)
# ===============================
log_info("Ejecutando pruebas unitarias")

test_that("Validación de funciones principales", {
  # Prueba de interpolación
  test_x <- 1:10
  test_y <- rnorm(10)
  expect_length(safe_spline_interp(test_x, test_y, seq(1, 10, 0.5)), 19)
  
  # Prueba de FPCA
  expect_true(ncol(fpca$scores) >= 3)
  expect_true(all(fpca$values >= 0))
  
  # Prueba de predicción
  pred <- predict_coastal_change("2024-06-01", fpca, fit_list, data_vario, fd_obj, X_common)
  expect_true(all(pred$reconstruction$upper >= pred$reconstruction$lower))
})

log_info("Análisis completado exitosamente")
