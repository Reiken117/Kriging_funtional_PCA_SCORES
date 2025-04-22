# ===================================================
# 1. Cargar librerías necesarias
# ===================================================
library(sf)          # Para manejo de datos geoespaciales
library(sp)          # Funciones espaciales básicas
library(ggplot2)     # Visualización gráfica
library(dplyr)       # Manipulación de datos
library(fda)         # Análisis funcional de datos
library(gstat)       # Geoestadística y kriging
library(plotly)      # Visualización 3D interactiva
library(viridis)     # Escalas de color
library(mvtnorm)     # Para manejo de distribuciones multivariadas y covarianzas
library(purrr)       # Para funciones de mapeo

# ===================================================
# 2. Cargar y preparar los datos 
# ===================================================
# Leer datos geoespaciales desde archivo JSON
gdf <- st_read("G:/Univalle/puntos_generados.geojson") 

# Verificación de la estructura de datos
if(!"DATE" %in% names(gdf)) stop("La columna DATE no existe en los datos")
if(nrow(gdf) == 0) stop("Datos geoespaciales vacíos")

# Uso de coordenadas originales y conversión de fecha
gdf_df <- gdf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(X = X, Y = Y) %>%
  mutate(DATE = as.Date(gdf$DATE))  # Conversión a tipo Date

# Visualización de los datos transformados
ggplot(gdf_df, aes(X, Y, color = DATE)) +
  geom_point() +
  scale_color_viridis(trans = "date") +
  labs(title = "Datos de línea costera - Coordenadas originales") +
  theme_minimal()

# ===================================================
# 3. Análisis Funcional y Suavizado
# ===================================================
# Definición de puntos comunes para interpolación
X_common <- seq(min(gdf_df$X), max(gdf_df$X), length.out = 500)

# Función de interpolación segura con manejo de NA
safe_interp <- function(x, y, xout) {
  tryCatch({
    approx(x, y, xout = xout, rule = 2)$y
  }, error = function(e) {
    warning("Error en interpolación para la fecha: ", unique(x))
    rep(NA, length(xout))
  })
}

# Generar Y_matrix con verificación de integridad
Y_matrix <- gdf_df %>%
  group_by(DATE) %>%
  group_map(~ {
    if (nrow(.x) >= 2) {
      safe_interp(.x$X, .x$Y, X_common)
    } else {
      warning("Fecha ", .x$DATE[1], " tiene menos de 2 puntos. Omitiendo.")
      rep(NA, length(X_common))
    }
  }) %>%
  discard(~ all(is.na(.x))) %>%  # Eliminar columnas completamente NA
  do.call(cbind, .)

if (ncol(Y_matrix) == 0) 
  stop("Y_matrix está vacío. Revise los datos de entrada.")
if (anyNA(Y_matrix)) 
  warning("Y_matrix contiene NA. Considerar imputación.")


# ===================================================
# 3.1 Optimización automática con optim
# ---------------------------------------
optimize_fdPar <- function(data, X_common) {
  error_fun <- function(params) {
    nbasis <- round(params[1])
    lambda <- 10^(params[2])
    basis <- create.bspline.basis(range(X_common), nbasis = nbasis, norder = 4)
    fdPar_obj <- fdPar(basis, lambda = lambda)
    sse_val <- smooth.basis(X_common, data, fdPar_obj)$SSE
    mean(sse_val)
  }
  init_params <- c(nbasis = 70, lambda = -8)
  res <- optim(init_params, error_fun, method = "Nelder-Mead", control = list(maxit = 200))
  list(nbasis = round(res$par[1]), lambda = 10^(res$par[2]))
}

# ===================================================
# 3.2 Análisis de sensibilidad (grid search)
# ---------------------------------------
sensitivity_analysis <- function(data, X_common, nbasis_vec = seq(60, 80, 2), 
                                 lambda_vec = 10^seq(-10, -5, 1)) {
  results <- expand.grid(nbasis = nbasis_vec, lambda = lambda_vec)
  results$error <- NA
  for (i in 1:nrow(results)) {
    nb <- results$nbasis[i]
    lam <- results$lambda[i]
    basis <- create.bspline.basis(range(X_common), nbasis = nb, norder = 4)
    fdPar_obj <- fdPar(basis, lambda = lam)
    sse_val <- smooth.basis(X_common, data, fdPar_obj)$SSE
    results$error[i] <- mean(sse_val)
  }
  return(results)
}

# Ejemplo de optimización (para la primera curva, por ejemplo)
optimal_params <- optimize_fdPar(Y_matrix[, 1], X_common)
cat("Parámetros óptimos (optim):\n- nbasis =", optimal_params$nbasis,
    "\n- lambda =", optimal_params$lambda, "\n")

# (Opcional: visualizar la sensibilidad)
sens_results <- sensitivity_analysis(Y_matrix[, 1], X_common)
ggplot(sens_results, aes(x = nbasis, y = log10(lambda), fill = error)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Superficie de Error: Sensibilidad en nbasis y lambda",
       x = "nbasis", y = "log10(lambda)") +
  theme_minimal()

# Ajuste final: Se utiliza la función cv_escala modificada o se elige la optimización
# Aquí empleamos el método original para mantener la funcionalidad
cv_escala <- function(data, n_basis_vec = seq(60, 80, 5), lambda_vec = 10^seq(-10, -5, 1)) {
  best_error <- Inf
  best_params <- list(nbasis = NA, lambda = NA)
  for (nb in n_basis_vec) {
    basis <- create.bspline.basis(range(X_common), nbasis = nb, norder = 4)
    for (lam in lambda_vec) {
      fdPar_obj <- fdPar(basis, lambda = lam)
      error <- mean(smooth.basis(X_common, data, fdPar_obj)$SSE)
      if (error < best_error) {
        best_error <- error
        best_params <- list(nbasis = nb, lambda = lam)
      }
    }
  }
  return(best_params)
}

# Optimizar parámetros
optimal_params <- cv_escala(Y_matrix)
final_basis <- create.bspline.basis(range(X_common), nbasis = optimal_params$nbasis, norder = 4)
final_fdPar <- fdPar(final_basis, lambda = optimal_params$lambda)
fd_obj <- smooth.basis(X_common, Y_matrix, final_fdPar)

if (any(is.na(fd_obj$fd$coefs))) {
  warning("Ajuste inestable. Incrementar lambda o reducir nbasis.")
} else {
  cat("Parámetros finales:\n- nbasis =", optimal_params$nbasis,
      "\n- lambda =", optimal_params$lambda, "\n")
}

# Visualización de funciones base y curvas ajustadas
plot(final_basis, main = "Funciones Base B-spline")
plot(fd_obj$fd, main = "Representación Funcional de la Línea de Costa",
     xlab = "Coordenada X", ylab = "Coordenada Y")

# ===================================================
# 4. FPCA y Selección de Componentes 
# ===================================================
fpca <- pca.fd(fd_obj$fd, nharm = 10)
colnames(fpca$scores) <- paste0("score.", 1:ncol(fpca$scores))
cumvar <- cumsum(fpca$values) / sum(fpca$values)
n_components <- which(cumvar >= 0.95)[1]
fpca$harmonics <- fpca$harmonics[1:n_components]
fpca$scores <- fpca$scores[, 1:n_components, drop = FALSE]

# ===================================================
# 5. Modelado Temporal con Kriging
# ===================================================
# ===================================================
# 5.1 Preparación de datos espaciales
# ===================================================
dates_numeric <- as.numeric(format(gdf_df$DATE, "%Y")) +
  as.numeric(format(gdf_df$DATE, "%j")) / 365

data_vario <- data.frame(
  time = dates_numeric,
  y = jitter(0, amount = 1e-6),
  fpca$scores
) %>%
  distinct(time, .keep_all = TRUE)
coordinates(data_vario) <- ~ time + y

if(!inherits(data_vario, "SpatialPointsDataFrame")) {
  stop("Error en conversión a objeto espacial")
}

# ===================================================
# 5.2 Función de ajuste de variogramas
# ===================================================
fit_variogram_robust_ext <- function(formula, sp_data) {
  environment(formula) <- environment()
  nombre_variable <- all.vars(formula)[1]
  if (!nombre_variable %in% names(sp_data)) {
    stop("Variable ", nombre_variable, " no encontrada")
  }
  valores <- sp_data[[nombre_variable]]
  var_inicial <- var(valores, na.rm = TRUE)
  rango_temporal <- diff(range(coordinates(sp_data)[,1], na.rm = TRUE))
  
  message("Variable: ", nombre_variable, " - Varianza inicial: ", var_inicial, 
          " - Rango temporal: ", rango_temporal)
  
  if (var_inicial == 0) {
    message("Varianza cero en ", nombre_variable, ". Usando modelo Nugget.")
    return(vgm(0, "Nug", 0))
  }
  
  # Modelos candidatos extendidos: Exp, Gau, dos Matérn y Nugget
  modelos <- list(
    vgm(0.6 * var_inicial, "Exp", rango_temporal * 2, 0.4 * var_inicial),
    vgm(0.6 * var_inicial, "Gau", rango_temporal / 2, 0.4 * var_inicial),
    vgm(0.7 * var_inicial, "Mat", rango_temporal / 3, 0.3 * var_inicial, kappa = 1.5),
    vgm(0.7 * var_inicial, "Mat", rango_temporal / 3, 0.3 * var_inicial, kappa = 2.5),
    vgm(var_inicial, "Nug", 0)
  )
  
  ancho <- rango_temporal / 15
  vgm_emp <- tryCatch({
    variogram(formula, data = sp_data, width = ancho, cutoff = rango_temporal * 0.8)
  }, error = function(e) {
    message("Error en cálculo del variograma empírico: ", e$message)
    return(NULL)
  })
  
  if (is.null(vgm_emp)) return(NULL)
  
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
    if (!is.null(a)) {
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
  
  indices_validos <- which(sapply(ajustes, function(a) !is.null(a) && a$range[2] > 0))
  if (length(indices_validos) > 0) {
    sse_validos <- sse[indices_validos]
    mejor_index <- indices_validos[which.min(sse_validos)]
    mejor_ajuste <- ajustes[[mejor_index]]
  } else {
    warning("No se encontró modelo adecuado para ", nombre_variable, ". Usando modelo Nugget.")
    mejor_ajuste <- vgm(var_inicial, "Nug", 0)
  }
  
  return(mejor_ajuste)
}

# ===================================================
# 5.3 Ajuste de variogramas por componente
# ===================================================
varianzas <- apply(fpca$scores, 2, var)
umbral_varianza <- 0.01 * max(varianzas)
componentes_activos <- which(varianzas > umbral_varianza)

fit_list <- lapply(componentes_activos, function(i) {
  componente <- paste0("score.", i)
  formula <- reformulate("1", response = componente)
  message("Ajustando componente: ", componente)
  tryCatch({
    fit_variogram_robust_ext(formula, data_vario)
  }, error = function(e) {
    message("Error en ", componente, ": ", e$message)
    return(NULL)
  })
})
if(any(sapply(fit_list, is.null))) {
  stop("Error en ajuste de componentes: ", paste(componentes_activos[sapply(fit_list, is.null)], collapse = ", "))
}

print(paste("Número de componentes activos:", length(componentes_activos)))
print(head(data_vario))
print(summary(varianzas))

# ===================================================
# 5.4 Validación y visualización de variogramas por componente
# ===================================================
x11()
if (length(componentes_activos) > 0) {
  dev.new()
  par(mfrow = c(ceiling(length(componentes_activos)/2), 2), mar = c(4,4,2,1))
  for(i in seq_along(componentes_activos)) {
    idx <- componentes_activos[i]
    vgm_emp <- variogram(reformulate("1", paste0("score.", idx)), data_vario)
    if(nrow(vgm_emp) > 0) {
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

# ===================================================
# 5.5 Ejemplo de Kriging Universal (con tendencia)
# ===================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# Función de predicción robusta con regularización y validación
predict_universal_curve_df <- function(target_dates, trend = TRUE) {
  # Conversión de fechas o tiempos
  if (inherits(target_dates, "Date")) {
    raw_target_times <- sapply(target_dates, function(d) {
      as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
    })
  } else {
    raw_target_times <- target_dates
    target_dates <- as.Date(paste0(floor(raw_target_times), "-01-01"))
  }
  
  # Extraer y reescalar los tiempos observados
  unique_dates <- sort(unique(gdf_df$DATE))
  raw_times <- sapply(unique_dates, function(d) {
    as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
  })
  time_range <- range(raw_times)
  observed_times <- (raw_times - time_range[1]) / diff(time_range)
  target_times <- (raw_target_times - time_range[1]) / diff(time_range)
  
  n_obs <- length(observed_times)
  sigma2 <- var(as.vector(Y_matrix), na.rm = TRUE)
  alpha <- diff(range(observed_times)) / 2
  cov_fun <- function(t1, t2) sigma2 * exp(-abs(t1 - t2) / alpha)
  
  # Calcular la matriz de covarianza y aplicar regularización (nugget)
  C <- outer(observed_times, observed_times, cov_fun)
  diag(C) <- diag(C) + 1e-6  # Agregar nugget para estabilizar la inversa
  
  # Construir la matriz extendida K_base según se incluya tendencia o no
  if (trend) {
    T_mat <- cbind(1, observed_times)
    K_base <- rbind(cbind(C, T_mat),
                    cbind(t(T_mat), matrix(0, ncol = ncol(T_mat), nrow = ncol(T_mat))))
  } else {
    one_vec <- rep(1, n_obs)
    K_base <- rbind(cbind(C, one_vec),
                    c(one_vec, 0))
  }
  
  # Verificar que K_base no tenga valores anómalos y, de ser necesario, regularizar
  if(any(!is.finite(K_base))) {
    stop("K_base contiene valores NA, NaN o Inf.")
  }
  if (kappa(K_base) > 1e12) {
    diag(K_base)[1:n_obs] <- diag(K_base)[1:n_obs] + 1e-6
  }
  
  # Ordenar según los tiempos observados
  order_idx <- order(observed_times)
  Y_matrix_ordered <- Y_matrix[, order_idx, drop = FALSE]
  
  # Predicción para cada tiempo objetivo
  res_list <- lapply(seq_along(target_times), function(i) {
    t_star <- target_times[i]
    c_vec <- sapply(observed_times, function(t) cov_fun(t, t_star))
    
    # Construcción del vector RHS
    rhs <- if (trend) {
      c(c_vec, c(1, t_star))
    } else {
      c(c_vec, 1)
    }
    
    # Resolver el sistema; si falla, se emplea la pseudoinversa
    sol <- tryCatch(solve(K_base, rhs),
                    error = function(e) MASS::ginv(K_base) %*% rhs)
    
    weights <- sol[1:n_obs]
    pred_curve <- as.vector(Y_matrix_ordered %*% weights)
    # Se utiliza tryCatch en la inversión de C para capturar potenciales errores
    variance <- sigma2 - t(c_vec) %*% tryCatch(solve(C, c_vec),
                                               error = function(e) {
                                                 stop("Error al invertir C: ", e$message)
                                               }) + sol[n_obs + 1]
    
    tibble(
      date = target_dates[i],
      coord = seq_along(pred_curve),
      value = pred_curve,
      variance = as.numeric(variance)
    )
  })
  
  bind_rows(res_list)
}

# ===================================================
# 5.5.1 Función de graficación
# -----------------------

plot_predicted_curves <- function(pred_df, ci = TRUE, alpha = 0.2) {
  # Definir el gráfico base con las curvas predichas
  p <- ggplot(pred_df, aes(x = coord, y = value, color = as.factor(date), group = date)) +
    geom_line(linewidth = 1) +
    labs(color = "Fecha", x = "Coordenada espacial", y = "Valor predicho") +
    theme_minimal()
  
  # Si se quiere mostrar la banda de confianza, se añaden las columnas upper y lower
  if (ci) {
    pred_df_ci <- pred_df %>%
      mutate(
        upper = value + 2 * sqrt(variance),
        lower = value - 2 * sqrt(variance)
      )
    # Añadir geom_ribbon con el data frame actualizado que contiene las columnas upper y lower
    p <- p + geom_ribbon(data = pred_df_ci,
                         aes(x = coord, ymin = lower, ymax = upper, fill = as.factor(date)),
                         alpha = alpha, color = NA)
  }
  
  return(p)
}

fechas <- as.Date(c("2025-01-01", "2026-01-01", "2027-01-01"))
pred_df <- predict_universal_curve_df(fechas, trend = TRUE)
plot_predicted_curves(pred_df, ci = TRUE)


# ===================================================
# 5.5.2  validación cruzada temporal
# -----------------------
# Función de validación cruzada para Kriging Universal con mejoras
cross_validate_kriging_universal <- function(years_ahead = 1:10, trend = TRUE) {
  unique_dates <- sort(unique(gdf_df$DATE))
  raw_times <- sapply(unique_dates, function(d) {
    as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
  })
  time_range <- range(raw_times)
  observed_times <- (raw_times - time_range[1]) / diff(time_range)
  
  order_idx <- order(observed_times)
  Y_ordered <- Y_matrix[, order_idx, drop = FALSE]
  times_ordered <- observed_times[order_idx]
  dates_ordered <- unique_dates[order_idx]
  
  results <- lapply(years_ahead, function(h) {
    errors <- c()
    n <- length(times_ordered)
    
    # Se inicia en 2 para asegurar al menos dos puntos en el entrenamiento
    for (i in 2:(n - h)) {
      training_idx <- seq_len(i)
      if (length(training_idx) < 2) next  # Omitir si no hay suficientes datos
      test_idx <- i + h
      if (test_idx > n) next
      
      Y_train <- Y_ordered[, training_idx, drop = FALSE]
      gdf_df_train <- gdf_df[gdf_df$DATE %in% dates_ordered[training_idx], ]
      
      # Redefinición local de la función de predicción usando los datos de entrenamiento reducidos
      pred_func <- function(target_date) {
        original_gdf <- gdf_df
        original_Y <- Y_matrix
        gdf_df <<- gdf_df_train
        Y_matrix <<- Y_train
        
        on.exit({
          gdf_df <<- original_gdf
          Y_matrix <<- original_Y
        }, add = TRUE)
        
        predict_universal_curve_df(target_date, trend = trend)
      }
      
      # Captura de errores en la predicción
      pred_df <- tryCatch(pred_func(dates_ordered[test_idx]),
                          error = function(e) {
                            message("Error en predicción para ", dates_ordered[test_idx], ": ", e$message)
                            return(NULL)
                          })
      if (is.null(pred_df)) next
      
      real_curve <- Y_ordered[, test_idx]
      rmse <- sqrt(mean((pred_df$value - real_curve)^2, na.rm = TRUE))
      errors <- c(errors, rmse)
    }
    
    tibble(horizon = h, rmse = mean(errors, na.rm = TRUE), sd = sd(errors, na.rm = TRUE))
  })
  
  bind_rows(results)
}

# Ejemplo de ejecución
cv_result <- cross_validate_kriging_universal(years_ahead = 1:8, trend = TRUE)
print(cv_result)


plot_cv_results <- function(cv_df) {
  ggplot(cv_df, aes(x = horizon, y = rmse)) +
    geom_line(linewidth = 1.2, color = "steelblue") +
    geom_point(size = 2, color = "steelblue") +
    geom_ribbon(aes(ymin = rmse - sd, ymax = rmse + sd), alpha = 0.2, fill = "lightblue") +
    labs(x = "Años hacia adelante", y = "RMSE", title = "Error de predicción vs horizonte temporal") +
    theme_minimal()
}

cv_result <- cross_validate_kriging_universal(years_ahead = 1:8)
plot_cv_results(cv_result)


# ===================================================
# 6. Predicción Funcional y Visualización (Kriging Funcional Puro)
# ===================================================

target_date <- as.Date("2025-01-01") #fecha que se busca estimar mediante kriging funcional puro

predict_functional_curve <- function(target_date) {
  if (inherits(target_date, "Date")) {
    target_time <- as.numeric(format(target_date, "%Y")) +
      as.numeric(format(target_date, "%j")) / 365
  } else {
    target_time <- target_date
  }
  unique_dates <- sort(unique(gdf_df$DATE))
  observed_times <- sapply(unique_dates, function(d) {
    as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
  })
  n_obs <- length(observed_times)
  sigma2 <- var(as.vector(Y_matrix), na.rm = TRUE)
  alpha <- diff(range(observed_times)) / 2
  cov_fun <- function(t1, t2) {
    sigma2 * exp(-abs(t1 - t2) / alpha)
  }
  C <- outer(observed_times, observed_times, cov_fun)
  c_vec <- sapply(observed_times, function(t) cov_fun(t, target_time))
  one_vec <- rep(1, n_obs)
  K <- rbind(cbind(C, one_vec),
             c(one_vec, 0))
  rhs <- c(c_vec, 1)
  sol <- solve(K, rhs)
  weights <- sol[1:n_obs]
  order_idx <- order(observed_times)
  Y_matrix_ordered <- Y_matrix[, order_idx]
  pred_curve <- as.vector(Y_matrix_ordered %*% weights)
  kriging_variance <- sigma2 - t(c_vec) %*% solve(C, c_vec) + sol[n_obs + 1]
  return(list(pred = pred_curve, variance = as.numeric(kriging_variance)))
}

func_kriging_result <- predict_functional_curve(target_date)
predicted_curve <- func_kriging_result$pred
predicted_variance <- func_kriging_result$variance
df_pred <- data.frame(
  X = X_common,
  Mean = predicted_curve,
  Lower = predicted_curve - 1.96 * sqrt(predicted_variance),
  Upper = predicted_curve + 1.96 * sqrt(predicted_variance)
)

plot_ly(data = df_pred, x = ~X) %>%
  add_ribbons(ymin = ~Lower, ymax = ~Upper,
              fillcolor = "rgba(0,100,80,0.2)",
              line = list(color = "rgba(255,255,255,0)"),
              name = "IC 95%") %>%
  add_lines(y = ~Mean, line = list(color = "blue"), name = "Predicción") %>%
  layout(title = paste("Predicción Funcional para", target_date),
         xaxis = list(title = "X Relativa"),
         yaxis = list(title = "Y Relativa"))

# ===================================================
# 7. Validación Cruzada Temporal para Kriging Funcional
# ===================================================
temporal_cv_functional <- function() {
  n_dates <- length(unique(gdf_df$DATE))
  rmse <- numeric(n_dates)
  coverage <- numeric(n_dates)
  
  unique_dates <- sort(unique(gdf_df$DATE))
  observed_times <- sapply(unique_dates, function(d) {
    as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
  })
  
  for (i in 1:n_dates) {
    idx <- setdiff(1:n_dates, i)
    order_idx <- order(observed_times)
    Y_matrix_ordered <- Y_matrix[, order_idx, drop = FALSE]
    Y_matrix_reducida <- Y_matrix_ordered[, idx, drop = FALSE]
    dates_reducidas <- observed_times[order_idx][idx]
    sigma2 <- var(as.vector(Y_matrix_reducida), na.rm = TRUE)
    alpha <- diff(range(dates_reducidas)) / 2
    cov_fun <- function(t1, t2) {
      sigma2 * exp(-abs(t1 - t2) / alpha)
    }
    n_obs <- length(dates_reducidas)
    C <- outer(dates_reducidas, dates_reducidas, cov_fun)
    target_time <- observed_times[order_idx][i]
    c_vec <- sapply(dates_reducidas, function(t) cov_fun(t, target_time))
    one_vec <- rep(1, n_obs)
    K <- rbind(cbind(C, one_vec), c(one_vec, 0))
    rhs <- c(c_vec, 1)
    sol <- solve(K, rhs)
    weights <- sol[1:n_obs]
    pred_curve <- as.vector(Y_matrix_reducida %*% weights)
    obs_curve <- Y_matrix_ordered[, i]
    rmse[i] <- sqrt(mean((pred_curve - obs_curve)^2, na.rm = TRUE))
    kriging_variance <- sigma2 - t(c_vec) %*% solve(C, c_vec) + sol[n_obs + 1]
    lower_bound <- pred_curve - 1.96 * sqrt(as.numeric(kriging_variance))
    upper_bound <- pred_curve + 1.96 * sqrt(as.numeric(kriging_variance))
    coverage[i] <- mean((obs_curve >= lower_bound) & (obs_curve <= upper_bound), na.rm = TRUE)
  }
  
  list(rmse = mean(rmse, na.rm = TRUE), coverage = mean(coverage, na.rm = TRUE))
}

cv_results_functional <- temporal_cv_functional()
message("\nResultados Validación Cruzada (Kriging Funcional):")
message("RMSE Promedio: ", round(cv_results_functional$rmse, 6))
message("Cobertura IC 95%: ", round(cv_results_functional$coverage * 100, 1), "%")

# ===================================================
# 8. Bootstrap para Intervalos de Confianza
# ===================================================
bootstrap_prediction <- function(target_date, n_boot = 100) {
  if (inherits(target_date, "Date")) {
    target_time <- as.numeric(format(target_date, "%Y")) +
      as.numeric(format(target_date, "%j")) / 365
  } else {
    target_time <- target_date
  }
  unique_dates <- sort(unique(gdf_df$DATE))
  observed_times <- sapply(unique_dates, function(d) {
    as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
  })
  n_obs <- length(observed_times)
  sigma2 <- var(as.vector(Y_matrix), na.rm = TRUE)
  alpha <- diff(range(observed_times)) / 2
  cov_fun <- function(t1, t2) {
    sigma2 * exp(-abs(t1 - t2) / alpha)
  }
  C <- outer(observed_times, observed_times, cov_fun)
  c_vec <- sapply(observed_times, function(t) cov_fun(t, target_time))
  one_vec <- rep(1, n_obs)
  K <- rbind(cbind(C, one_vec),
             c(one_vec, 0))
  rhs <- c(c_vec, 1)
  sol <- solve(K, rhs)
  weights <- sol[1:n_obs]
  order_idx <- order(observed_times)
  Y_matrix_ordered <- Y_matrix[, order_idx]
  pred_curve <- as.vector(Y_matrix_ordered %*% weights)
  
  boot_preds <- replicate(n_boot, {
    boot_idx <- sample(1:n_obs, replace = TRUE)
    C_boot <- outer(observed_times[boot_idx], observed_times[boot_idx], cov_fun)
    c_vec_boot <- sapply(observed_times[boot_idx], function(t) cov_fun(t, target_time))
    one_boot <- rep(1, length(boot_idx))
    K_boot <- rbind(cbind(C_boot, one_boot), c(one_boot, 0))
    rhs_boot <- c(c_vec_boot, 1)
    sol_boot <- tryCatch(solve(K_boot, rhs_boot), error = function(e) rep(NA, length(boot_idx) + 1))
    weights_boot <- sol_boot[1:length(boot_idx)]
    Y_boot <- Y_matrix_ordered[, boot_idx, drop = FALSE]
    as.vector(Y_boot %*% weights_boot)
  })
  
  lower <- apply(boot_preds, 1, quantile, probs = 0.025, na.rm = TRUE)
  upper <- apply(boot_preds, 1, quantile, probs = 0.975, na.rm = TRUE)
  return(list(pred = pred_curve, lower = lower, upper = upper))
}

boot_result <- bootstrap_prediction(target_date, n_boot = 200)
df_boot <- data.frame(
  X = X_common,
  Pred = boot_result$pred,
  Lower = boot_result$lower,
  Upper = boot_result$upper
)

plot_ly(data = df_boot, x = ~X) %>%
  add_ribbons(ymin = ~Lower, ymax = ~Upper,
              fillcolor = "rgba(0,100,80,0.2)",
              line = list(color = "rgba(255,255,255,0)"),
              name = "IC 95% Bootstrap") %>%
  add_lines(y = ~Pred, line = list(color = "blue"), name = "Predicción") %>%
  layout(title = paste("Predicción Funcional con Bootstrap para", target_date),
         xaxis = list(title = "X Relativa"),
         yaxis = list(title = "Y Relativa"))

# ===================================================
# 9. Validación Estadística y Diagnóstico 
# ===================================================
# 9.1 Matriz de Covarianza entre Curvas Funcionales
observed_times <- sapply(sort(unique(gdf_df$DATE)), function(d) {
  as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
})
ordered_Y <- Y_matrix[, order(observed_times)]
cov_matrix <- cov(t(ordered_Y), use = "pairwise.complete.obs")
heatmap(cov_matrix, symm = TRUE, col = viridis(100), 
        main = "Matriz de Covarianza entre Curvas Funcionales")

# 9.2 Autocorrelación Temporal (lag-1)
autocorr_lag1 <- sapply(2:ncol(ordered_Y), function(i) {
  cor(ordered_Y[, i], ordered_Y[, i - 1], use = "complete.obs")
})
mean_autocorr <- mean(autocorr_lag1, na.rm = TRUE)
cat("Autocorrelación media lag-1:", round(mean_autocorr, 3), "\n")

# 9.3 Validación Cruzada Extendida (Errores e IC por Fecha)
unique_dates <- sort(unique(gdf_df$DATE))
cv_resultados <- data.frame(
  Fecha = unique_dates,
  RMSE = NA,
  IC95_Cobertura = NA
)
for (i in seq_along(unique_dates)) {
  fecha_eval <- unique_dates[i]
  target_curve <- Y_matrix[, i]
  Y_reducida <- Y_matrix[, -i]
  fechas_reducidas <- unique_dates[-i]
  pred_result <- tryCatch({
    predict_functional_curve(fecha_eval)
  }, error = function(e) return(NULL))
  if (!is.null(pred_result)) {
    error_curve <- target_curve - pred_result$pred
    rmse_i <- sqrt(mean(error_curve^2, na.rm = TRUE))
    sd_pred <- sqrt(pred_result$variance)
    lower <- pred_result$pred - 1.96 * sd_pred
    upper <- pred_result$pred + 1.96 * sd_pred
    dentro <- (target_curve >= lower) & (target_curve <= upper)
    cobertura <- mean(dentro, na.rm = TRUE)
    cv_resultados$RMSE[i] <- rmse_i
    cv_resultados$IC95_Cobertura[i] <- cobertura
  }
}
print(summary(cv_resultados))

# 9.4 Gráficos de Validación
ggplot(cv_resultados, aes(x = Fecha)) +
  geom_line(aes(y = RMSE), color = "tomato", size = 1.2) +
  geom_point(aes(y = RMSE), color = "tomato") +
  geom_hline(yintercept = mean(cv_resultados$RMSE, na.rm = TRUE), 
             linetype = "dashed", color = "gray40") +
  labs(title = "RMSE por Fecha (Validación Cruzada)",
       y = "RMSE", x = "Fecha") +
  theme_minimal()

ggplot(cv_resultados, aes(x = Fecha)) +
  geom_line(aes(y = IC95_Cobertura), color = "blue", size = 1.2) +
  geom_point(aes(y = IC95_Cobertura), color = "blue") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
  labs(title = "Cobertura Empírica de IC 95% por Fecha",
       y = "Cobertura", x = "Fecha") +
  theme_minimal()
