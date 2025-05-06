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
library(tidyr)
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
# ===================================================
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
# ===================================================
sensitivity_analysis <- function(data, X_common, nbasis_vec = seq(20, 80, 2), 
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
# 3.3 Análisis Funcional con Selección de Modelos por AIC
# ===================================================

# Función para ajustar y seleccionar el mejor modelo por AIC
ajustar_mejor_modelo <- function(x, y) {
  datos <- data.frame(x = x, y = y)
  datos <- datos[complete.cases(datos), ]
  
  if (nrow(datos) < 4) {
    warning("No hay suficientes puntos para ajustar modelos. Usando interpolación lineal.")
    return(approx(datos$x, datos$y, xout = x, rule = 2)$y)
  }
  
  tryCatch({
    # Ajustamos ambos modelos
    m_ps <- mgcv::gam(y ~ s(x, bs = "ps"), data = datos, method = "GCV.Cp")
    m_bs <- mgcv::gam(y ~ s(x, bs = "bs"), data = datos, method = "GCV.Cp")
    
    # Comparamos por AIC
    aic_ps <- AIC(m_ps)
    aic_bs <- AIC(m_bs)
    
    # Seleccionamos el mejor modelo
    if (aic_ps < aic_bs) {
      predict(m_ps, newdata = data.frame(x = x))
    } else {
      predict(m_bs, newdata = data.frame(x = x))
    }
  }, error = function(e) {
    warning("Error en ajuste GAM: ", e$message, ". Usando interpolación lineal.")
    approx(datos$x, datos$y, xout = x, rule = 2)$y
  })
}

# Aplicamos a cada columna de Y_matrix
Y_smooth <- purrr::map_dfc(1:ncol(Y_matrix), function(j) {
  yj <- Y_matrix[, j]
  smooth_curve <- ajustar_mejor_modelo(X_common, yj)
  data.frame(smooth_curve)
}) %>% as.matrix()
# Renombramos columnas para mantener referencia a fechas
colnames(Y_smooth) <- colnames(Y_matrix)

# Análisis de qué modelos fueron seleccionados
model_selection <- sapply(1:ncol(Y_matrix), function(j) {
  yj <- Y_matrix[, j]
  datos <- data.frame(x = X_common, y = yj)
  datos <- datos[complete.cases(datos), ]
  
  m_ps <- mgcv::gam(y ~ s(x, bs = "ps"), data = datos, method = "GCV.Cp")
  m_bs <- mgcv::gam(y ~ s(x, bs = "bs"), data = datos, method = "GCV.Cp")
  
  ifelse(AIC(m_ps) < AIC(m_bs), "P-spline", "B-spline")
})

cat("\nDistribución de modelos seleccionados:\n")
print(table(model_selection))

# Convertir a formato largo para ggplot
Y_smooth_df <- as.data.frame(Y_smooth)
colnames(Y_smooth_df) <- unique(gdf_df$DATE)[1:ncol(Y_smooth)]  # Asegurar misma longitud
Y_smooth_df$X <- X_common

long_data <- Y_smooth_df %>%
  tidyr::pivot_longer(-X, names_to = "DATE", values_to = "Y") %>%
  mutate(DATE = as.Date(DATE))

# 3. Visualización mejorada con facetas por año
library(lubridate)  # Para manejo de fechas

# Añadir columnas de año y mes
long_data <- long_data %>%
  mutate(Year = year(DATE),
         Month = month(DATE, label = TRUE))

# Gráfico por años con colores por mes
ggplot(long_data, aes(X, Y, group = DATE, color = Month)) +
  geom_line(alpha = 0.7) +
  scale_color_viridis_d() +
  facet_wrap(~Year, scales = "free_y") +
  labs(title = "Líneas Costeras Suavizadas por Fecha",
       subtitle = "Separado por años, colores representan meses",
       x = "Coordenada X", y = "Coordenada Y") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
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

estimate_alpha_MLE <- function(observed_times, Y_matrix) {
  Y_centered <- scale(Y_matrix, scale = FALSE)
  
  neg_loglik <- function(params) {
    alpha <- exp(params[1])
    sigma2 <- exp(params[2])
    
    C <- sigma2 * exp(-outer(observed_times, observed_times, function(t1,t2) abs(t1-t2))/alpha)
    diag(C) <- diag(C) + 1e-6 * median(diag(C))
    
    tryCatch({
      -sum(mvtnorm::dmvnorm(Y_centered, mean = rep(0, ncol(Y_centered)), sigma = C, log = TRUE)) # ✅ Suma crítica
    }, error = function(e) 1e10)
  }
  
  init_params <- log(c(diff(range(observed_times))/2, var(as.vector(Y_matrix))))
  opt <- optim(init_params, neg_loglik, method = "Nelder-Mead", control = list(maxit = 10000))
  
  return(list(alpha = exp(opt$par[1]), sigma2 = exp(opt$par[2])))
}


# Función de predicción robusta con regularización y validación
predict_universal_curve_df <- function(target_dates, trend = TRUE, verbose = FALSE) {
  # Conversión de fechas a tiempos normalizados
  unique_dates <- sort(unique(gdf_df$DATE))
  raw_times <- sapply(unique_dates, function(d) as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j"))/365)
  time_range <- range(raw_times)
  observed_times <- (raw_times - time_range[1]) / diff(time_range)
  
  # Estimación MLE
  cov_params <- estimate_alpha_MLE(observed_times, Y_matrix)
  alpha <- cov_params$alpha; sigma2 <- cov_params$sigma2
  if(verbose) message(sprintf("Universal: alpha=%.4f, sigma2=%.4f", alpha, sigma2))
  cov_fun <- function(t1, t2) sigma2 * exp(-abs(t1 - t2)/alpha)
  
  # Construcción de C y nugget extra si hace falta
  C <- outer(observed_times, observed_times, cov_fun)
  diag(C) <- diag(C) + 1e-6
  if(kappa(C) > 1e10) diag(C) <- diag(C) + 1e-4
  
  # Matriz extendida
  n_obs <- length(observed_times)
  if(trend) {
    T_mat <- cbind(1, observed_times)
    K_base <- rbind(cbind(C, T_mat), cbind(t(T_mat), matrix(0,2,2)))
  } else {
    one_vec <- rep(1, n_obs)
    K_base <- rbind(cbind(C, one_vec), c(one_vec, 0))
  }
  
  # Regularización de K_base
  if(any(!is.finite(K_base))) stop("K_base contiene valores no finitos")
  if(kappa(K_base) > 1e12) diag(K_base)[1:n_obs] <- diag(K_base)[1:n_obs] + 1e-3
  
  # Loop de predicción
  out <- lapply(target_dates, function(d) {
    t_star_raw <- as.numeric(format(as.Date(d), "%Y")) + as.numeric(format(as.Date(d), "%j"))/365
    t_star <- (t_star_raw - time_range[1]) / diff(time_range)
    # RHS
    c_vec <- sapply(observed_times, function(t) cov_fun(t, t_star))
    rhs <- if(trend) c(c_vec, 1, t_star) else c(c_vec, 1)
    
    sol <- tryCatch(solve(K_base, rhs), error = function(e) MASS::ginv(K_base) %*% rhs)
    weights <- sol[1:n_obs]
    pred_curve <- as.vector(Y_matrix[, order(observed_times)] %*% weights)
    var_term <- sigma2 - as.numeric(t(c_vec) %*% solve(C, c_vec)) + sol[n_obs + ifelse(trend,2,1)]
    tibble(date = as.Date(d), coord = seq_along(pred_curve), value = pred_curve, variance = var_term)
  })
  bind_rows(out)
}

# ===================================================
# 5.5.1 Función de graficación
# ===================================================

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



# Datos de prueba
test_times <- seq(0, 10, length.out = 50)
C <- 1.0 * exp(-as.matrix(dist(test_times))/1.0)
Y_sim <- rmvnorm(100, mean = rep(0, 50), sigma = C)  # 100x50

# Estimar parámetros
params <- estimate_alpha_MLE(test_times, Y_sim)
cat(sprintf("alpha=%.2f (1.00)\nsigma²=%.2f (1.00)", params$alpha, params$sigma2))

# ===================================================
# 5.5.2  validación cruzada temporal
# ===================================================
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
# 5.6 Kriging Funcional Ordinario y Covarianza  OKFD
# ===================================================
# 5.6.1 Trace-variograma experimental
delta_t <- diff(X_common)[1]
trace_vario <- function(i, j) {
  fi <- eval.fd(X_common, fd_obj$fd)[, i]
  fj <- eval.fd(X_common, fd_obj$fd)[, j]
  0.5 * sum((fi - fj)^2) * delta_t
}

dates     <- sort(unique(gdf_df$DATE))
dates_num <- sapply(dates, function(d) year(d) + yday(d)/365)

pairs <- expand.grid(i = seq_along(dates), j = seq_along(dates)) %>%
  filter(i < j) %>%
  mutate(
    h     = abs(dates_num[i] - dates_num[j]),
    gamma = mapply(trace_vario, i, j)
  )

psill0 <- var(pairs$gamma)
max_h  <- max(pairs$h)

# 5.6.2 Variograma 2D (tiempo vs espacio dummy)
exp_vario <- pairs %>%
  mutate(space = 0,
         time  = h) %>%
  group_by(time, space) %>%
  summarise(gamma = mean(gamma), .groups = "drop")

# Convertir a Spatial
coordinates(exp_vario) <- ~ time + space

# Ajuste de variograma exponencial sin nugget (o con nugget si se desea)
best_model_2d <- vgm(psill = psill0, model = "Exp", range = max_h/3)

# 5.6.3 Preparar matriz de curvas
curve_mat   <- eval.fd(X_common, fd_obj$fd)[, seq_along(dates)]
curve_mat_t <- t(curve_mat)
n_curvas <- nrow(curve_mat_t)

df_scores <- data.frame(
  space = rep(0, n_curvas),
  time  = dates_num,
  curve_mat_t
)
colnames(df_scores) <- c("space", "time", paste0("V", seq_len(ncol(curve_mat_t))))
coordinates(df_scores) <- ~ time + space

# 5.6.4 Función de Kriging Funcional Ordinario
krige_functional <- function(target_date, model2d = best_model_2d) {
  t0     <- year(target_date) + yday(target_date)/365
  newloc <- data.frame(space = 0, time = t0)
  coordinates(newloc) <- ~ time + space
  
  preds <- sapply(seq_len(ncol(curve_mat_t)), function(k) {
    vn <- paste0("V", k)
    kr <- krige(
      formula   = as.formula(paste(vn, "~ 1")),
      locations = df_scores,
      newdata   = newloc,
      model     = model2d,
      nmax      = 20
    )
    kr$var1.pred
  })
  preds
}

# Prueba rápida
pred_ord <- krige_functional(as.Date("2025-01-01"))
plot(X_common, pred_ord, type = "l",
     main = "Kriging Funcional Ordinario (2025-01-01)",
     xlab = "X_common", ylab = "Y_predicha")

# ===================================================
# 5.7. (Kriging Funcional Puro) kriging de traza
# ===================================================

target_date <- as.Date("2025-01-01")

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
# 5.8 Función genérica de cross-validación para los 3 métodos
# ===================================================
cross_validate_all_models <- function(years_ahead = 1:5, trend = TRUE) {
  models <- list(
    universal = function(date) predict_universal_curve_df(date, trend = trend),
    OKFD      = function(date) {
      vec <- krige_functional(date)
      tibble(date = as.Date(date), coord = X_common, value = vec)
    },
    traza     = function(date) {
      pr <- predict_functional_curve(date)
      tibble(date = as.Date(date), coord = X_common, value = pr$pred)
    }
  )
  
  unique_dates <- sort(unique(gdf_df$DATE))
  dates_ordered <- unique_dates  # Ya están ordenadas cronológicamente
  
  results <- list()
  for(h in years_ahead) {
    for(i in seq_along(dates_ordered)) {
      train_dates <- dates_ordered[1:i]
      if(length(train_dates) == 0) next
      
      last_train_date <- max(train_dates)
      target_test_date <- last_train_date + lubridate::years(h)
      test_date_candidates <- which(dates_ordered >= target_test_date)
      
      if(length(test_date_candidates) == 0) next
      test_date <- dates_ordered[test_date_candidates[1]]
      
      # Verificar diferencia mínima de h años
      time_diff <- as.numeric(difftime(test_date, last_train_date, units = "days")) / 365.25
      if(time_diff < h) next
      
      # Subset datos de entrenamiento
      old_gdf <- gdf_df
      old_Y <- Y_matrix
      gdf_df <<- gdf_df[gdf_df$DATE %in% train_dates, ]
      Y_matrix <<- Y_matrix[, 1:i, drop = FALSE]  # Asumiendo columnas ordenadas
      
      # Evaluar modelos
      for(name in names(models)) {
        pred_df <- tryCatch(models[[name]](test_date), error = function(e) NULL)
        if(is.null(pred_df)) next
        
        real_curve <- old_Y[, which(dates_ordered == test_date)]
        rmse <- sqrt(mean((pred_df$value - real_curve)^2, na.rm = TRUE))
        
        results[[length(results)+1]] <- tibble(
          model = name,
          horizon = h,
          date = test_date,
          rmse = rmse
        )
      }
      
      # Restaurar datos globales
      gdf_df <<- old_gdf
      Y_matrix <<- old_Y
    }
  }
  bind_rows(results)
}

# Ejecución de la validación para horizontes 1:8
cv_all <- cross_validate_all_models(years_ahead = 1:8, trend = TRUE)
print(cv_all)
ggplot(cv_all, aes(x = horizon, y = rmse, color = model)) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(title = "Error de predicción por horizonte temporal", x = "Años hacia adelante", y = "RMSE")

cv_summary <- cv_all %>%
  group_by(model, horizon) %>%
  summarise(mean_rmse = mean(rmse, na.rm = TRUE))

ggplot(cv_summary, aes(x = horizon, y = mean_rmse, color = model)) + 
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(title = "RMSE promedio por horizonte temporal", x = "Años hacia adelante", y = "RMSE")


table(cv_all$horizon, cv_all$model)
# ===================================================
# 5.9 Comparación de métodos de Kriging Funcional 
# ===================================================
# Fecha de predicción (ajusta según necesidad)
fecha_obj <- as.Date("2025-06-01")

# Predicciones de los tres métodos
pred_univ_df <- predict_universal_curve_df(fecha_obj, trend = TRUE) %>% arrange(coord)
pred_ord_vec <- krige_functional(fecha_obj)
pred_trace_res <- predict_functional_curve(fecha_obj)

# Combinar resultados en un data.frame
# Usamos el mismo vector X_common como coordenada espacial
df_compare <- data.frame(
  coord     = X_common,
  Universal = pred_univ_df$value,
  Ordinario = pred_ord_vec,
  Traza     = pred_trace_res$pred
)

# Transformar a formato largo para ggplot
library(tidyr)
long_compare <- df_compare %>%
  pivot_longer(-coord, names_to = "Metodo", values_to = "Valor")

# Gráfico comparativo de los tres métodos
library(ggplot2)
ggplot(long_compare, aes(x = coord, y = Valor, color = Metodo, linetype = Metodo)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Comparación de Métodos de Kriging Funcional",
    x = "Coordenada espacial", y = "Valor predicho",
    color = "Método", linetype = "Método"
  ) +
  theme_minimal()

# Cálculo de la diferencia absoluta media (MAD) entre pares de métodos
mad_univ_ord   <- mean(abs(df_compare$Universal - df_compare$Ordinario))
mad_univ_traza <- mean(abs(df_compare$Universal - df_compare$Traza))
mad_ord_traza  <- mean(abs(df_compare$Ordinario - df_compare$Traza))

cat(sprintf(
  "MAD Universal vs Ordinario: %.8f\n" , mad_univ_ord    ), sep = "")
cat(sprintf(
  "MAD Universal vs Traza:   %.8f\n" , mad_univ_traza  ), sep = "")
cat(sprintf(
  "MAD Ordinario vs Traza:   %.8f\n" , mad_ord_traza   ), sep = "")
# Resumen de diferencias absolutas medias (MAD)
mad_summary <- data.frame(
  Pair = c("Universal-Ordinario", "Universal-Traza", "Ordinario-Traza"),
  MAD  = c(mad_univ_ord, mad_univ_traza, mad_ord_traza)
)
print(mad_summary)
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
# --------------------------------------------------
observed_times <- sapply(sort(unique(gdf_df$DATE)), function(d) {
  as.numeric(format(d, "%Y")) + as.numeric(format(d, "%j")) / 365
})

ordered_Y <- Y_matrix[, order(observed_times)]
cov_matrix <- cov(t(ordered_Y), use = "pairwise.complete.obs")

# Visualización con heatmap
heatmap(cov_matrix, 
        symm = TRUE, 
        col = viridis(100), 
        main = "Matriz de Covarianza entre Curvas Funcionales",
        xlab = "Curvas", 
        ylab = "Curvas")

# 9.2 Autocorrelación Temporal (lag-1)
# ------------------------------------
autocorr_lag1 <- sapply(2:ncol(ordered_Y), function(i) {
  cor(ordered_Y[, i], ordered_Y[, i - 1], use = "complete.obs")
})

mean_autocorr <- mean(autocorr_lag1, na.rm = TRUE)
cat("\nAutocorrelación media lag-1:", round(mean_autocorr, 3), "\n")

# 9.3 Validación Cruzada Extendida (Errores e IC por Fecha)
# --------------------------------------------------------
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
  
  pred_result <- tryCatch({
    predict_functional_curve(fecha_eval)
  }, error = function(e) return(NULL))
  
  if (!is.null(pred_result)) {
    # Cálculo de RMSE
    error_curve <- target_curve - pred_result$pred
    rmse_i <- sqrt(mean(error_curve^2, na.rm = TRUE))
    
    # Cálculo de cobertura del IC 95%
    sd_pred <- sqrt(pred_result$variance)
    lower <- pred_result$pred - 1.96 * sd_pred
    upper <- pred_result$pred + 1.96 * sd_pred
    dentro <- (target_curve >= lower) & (target_curve <= upper)
    cobertura <- mean(dentro, na.rm = TRUE)
    
    cv_resultados$RMSE[i] <- rmse_i
    cv_resultados$IC95_Cobertura[i] <- cobertura
  }
}

# Resumen de resultados
cat("\nResumen de Validación Cruzada:\n")
print(summary(cv_resultados))

# 9.4 Gráficos de Validación
# --------------------------
# Gráfico de RMSE
ggplot(cv_resultados, aes(x = Fecha)) +
  geom_line(aes(y = RMSE), color = "tomato", size = 1.2) +
  geom_point(aes(y = RMSE), color = "tomato") +
  geom_hline(yintercept = mean(cv_resultados$RMSE, na.rm = TRUE), 
             linetype = "dashed", color = "gray40") +
  labs(title = "RMSE por Fecha (Validación Cruzada)",
       y = "RMSE", x = "Fecha") +
  theme_minimal()

# Gráfico de cobertura de IC
ggplot(cv_resultados, aes(x = Fecha)) +
  geom_line(aes(y = IC95_Cobertura), color = "blue", size = 1.2) +
  geom_point(aes(y = IC95_Cobertura), color = "blue") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
  labs(title = "Cobertura Empírica de IC 95% por Fecha",
       y = "Cobertura", x = "Fecha") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1))


# ===================================================
# 10. Validación de Supuestos y Métricas Clave
# ===================================================

# 10.1 Calidad y estructura de los datos
# --------------------------------------
# (a) Verificación de columnas requeridas
cols_req <- c("DATE", "X", "Y")
missing_cols <- setdiff(cols_req, names(gdf_df))
if (length(missing_cols) > 0) {
  stop("Faltan columnas requeridas: ", paste(missing_cols, collapse = ", "))
}

# (b) Proporción de valores faltantes
na_prop <- sapply(gdf_df[, cols_req], function(v) mean(is.na(v)))
cat("\nProporción de valores faltantes:\n")
print(na_prop)

# (c) Análisis de fechas
na_dates <- sum(is.na(gdf_df$DATE))
date_range <- range(gdf_df$DATE, na.rm = TRUE)
cat("\nNA en fechas:", na_dates)
cat("\nRango temporal:", format(date_range[1], "%Y-%m-%d"), "a", 
    format(date_range[2], "%Y-%m-%d"), "\n")

# (d) Densidad espacial
library(spatstat)
pp <- ppp(gdf_df$X, gdf_df$Y,
          window = owin(xrange = range(gdf_df$X), yrange = range(gdf_df$Y)))
plot(density(pp), 
     main = "Mapa de densidad de puntos",
     xlab = "Coordenada X",
     ylab = "Coordenada Y")

# 10.2 Construcción de Y_matrix
# -----------------------------
cat("\nDimensiones de Y_matrix:", dim(Y_matrix))
cat("\nProporción de NA en Y_matrix:", round(mean(is.na(Y_matrix)), 4), "\n")

# 10.3 Integridad de la interpolación
# -----------------------------------
na_after_interp <- apply(Y_matrix, 2, function(col) mean(is.na(col)))
cat("\nNA post-interpolación (min, med, max):",
    round(min(na_after_interp), 3), 
    round(median(na_after_interp), 3),
    round(max(na_after_interp), 3), "\n")

# 10.5 Análisis FPCA
# ------------------
# (a) Varianza explicada
cumvar <- cumsum((fpca$values) / sum(fpca$values))
cat("\nComponentes para ≥95% varianza:", which(cumvar >= 0.95)[1], "\n")

# (b) Gráfico de decaimiento de eigenvalores
plot(fpca$values, type = "b", 
     main = "Decaimiento de eigenvalores",
     xlab = "Componente",
     ylab = "Eigenvalor")

# 10.6 Análisis de variogramas
# ----------------------------
var_initial <- var(data_vario[[paste0("score.", componentes_activos[1])]], na.rm = TRUE)
fit0 <- fit_list[[1]]
cat("\nSemivariograma - Componente 1:\n")
cat("  Varianza inicial:", signif(var_initial, 3), "\n")
cat("  Nugget:", signif(fit0$psill[1], 3), "\n")
cat("  Sill:", signif(sum(fit0$psill), 3), "\n")
cat("  Range:", signif(fit0$range[2], 3), "\n")

# 10.7 Validación de Kriging
# --------------------------
# (a) Relación Nugget-to-Sill
nugget <- fit0$psill[1]
sill <- sum(fit0$psill)
cat("\nNugget-to-Sill ratio:", round(nugget/sill, 3), "\n")

# (b) Métricas globales
rmse_overall <- mean(cv_result$rmse, na.rm = TRUE)
coverage_overall <- mean(
  (cv_result$rmse - cv_result$sd <= cv_result$rmse) & 
    (cv_result$rmse + cv_result$sd >= cv_result$rmse), 
  na.rm = TRUE
)

cat("\nMétricas globales de predicción:\n")
cat("  RMSE:", signif(rmse_overall, 3), "\n")
cat("  Cobertura IC95%:", round(coverage_overall*100, 1), "%\n")


# ===================================================
# 11. Ajustes Avanzados
# ===================================================

# 11.1 Matriz de penalización
Pmat_correct <- eval.penalty(final_basis, Lfdobj = 2)
cond_P_correct <- kappa(Pmat_correct)
cat("\nNúmero de condición de P (segundo orden):", 
    format(cond_P_correct, scientific = TRUE, digits = 3), "\n")

# 11.2 Optimización de lambda
# ---------------------------
lambda_grid <- 10^seq(-12, 0, length.out = 50)
gcv_mean <- numeric(length(lambda_grid))

for (i in seq_along(lambda_grid)) {
  lam <- lambda_grid[i]
  fdPar_i <- fdPar(final_basis, lambda = lam)
  sb_i <- smooth.basis(X_common, Y_matrix, fdPar_i)
  gcv_mean[i] <- mean(sb_i$gcv, na.rm = TRUE)
}

# Resultados de optimización
lambda_opt_gcv <- lambda_grid[which.min(gcv_mean)]
cat("\nLambda óptimo (GCV):", format(lambda_opt_gcv, scientific = TRUE), "\n")

# 11.3 RMSE relativo
# ------------------
Y_range <- max(gdf_df$Y, na.rm = TRUE) - min(gdf_df$Y, na.rm = TRUE)
rmse_rel <- rmse_overall / Y_range * 100
cat("\nRMSE relativo:", round(rmse_rel, 2), "% del rango de Y\n")

# 11.4 Ajuste de intervalos de confianza al 90%
# ---------------------------------------------
multiplier <- qnorm(0.95)  # 1.645 para IC 90%

pred_df_adj <- pred_df %>%
  mutate(
    lower_90 = value - multiplier * sqrt(variance),
    upper_90 = value + multiplier * sqrt(variance)
  )

# Visualización
ggplot(pred_df_adj, aes(x = coord, y = value, group = factor(date))) +
  geom_line(aes(color = factor(date))) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = factor(date)), 
              alpha = 0.2) +
  labs(
    title = "Predicciones con Intervalos de Confianza al 90%",
    x = "Coordenada espacial",
    y = "Valor predicho",
    color = "Fecha", 
    fill = "Fecha"
  ) +
  theme_minimal()


# ===================================================



# Carga del paquete fda para construir bases funcionales
library(fda)

# Dominio de definición (aquí [0,1])
rangeval <- c(0, 1)

# Rango de nbasis a explorar (de 20 a 80)
nbasis_vals <- 20:80

# Rango de lambda en escala logarítmica (de 1e-10 a 1e0)
lambda_vals <- 10^seq(-10, 0, length.out = 11)

# Inicializar variables para la mejor combinación encontrada
bestCond <- Inf
bestType <- NULL
bestNbasis <- NA
bestLambda <- NA

# Lista de tipos de bases a considerar
# (omitiendo las gaussianas por no existir función directa en fda)
tipos_base <- c('bspline', 'fourier', 'monomial')

for (tipo in tipos_base) {
  for (nb in nbasis_vals) {
    # Crear objeto de base correspondiente
    if (tipo == 'bspline') {
      basis_obj <- create.bspline.basis(rangeval, nbasis = nb)
    } else if (tipo == 'fourier') {
      # Fourier: asegurar un número impar de funciones base para incluir constante
      nb_use <- ifelse(nb %% 2 == 0, nb + 1, nb)
      basis_obj <- create.fourier.basis(rangeval, nbasis = nb_use)
    } else if (tipo == 'monomial') {
      basis_obj <- create.monomial.basis(rangeval, nbasis = nb)
    }
    # Calcular matriz de Gram y matriz de penalización (segunda derivada)
    Cmat <- inprod(basis_obj, basis_obj, Lfdobj1 = 0, Lfdobj2 = 0, rng = rangeval)
    Rmat <- eval.penalty(basis_obj, Lfdobj = 2, rng = rangeval)
    # Explorar valores de lambda y computar condición
    for (lam in lambda_vals) {
      M <- Cmat + lam * Rmat
      # Calcular número de condición (manejo de posibles errores numéricos)
      cond_val <- try(kappa(M), silent = TRUE)
      if (inherits(cond_val, 'try-error')) next
      # Revisar estabilidad numérica
      if (cond_val < 1e12) {
        # Actualizar mejor combinación si mejora el número de condición
        if (cond_val < bestCond) {
          bestCond <- cond_val
          bestType <- tipo
          bestNbasis <- nb
          bestLambda <- lam
        }
      }
    }
  }
}

# Mostrar resultado final o diagnóstico
if (is.null(bestType)) {
  cat('No se encontró una combinación estable.\n')
} else {
  cat('Mejor base encontrada:', bestType, '\n')
  cat('Nbasis óptimo:', bestNbasis, '\n')
  cat('Lambda óptimo:', bestLambda, '\n')
  cat('Número de condición obtenido:', bestCond, '\n')
}


library(fda)

# Dominio [0, 1] y puntos de evaluación
X <- seq(0, 1, length.out = 200)
Y <- sin(2 * pi * X)  # función a reconstruir
Y_noisy <- Y + rnorm(length(Y), 0, 0.05)  # versión con ruido

# Datos para smooth.basis
rangeval <- c(0, 1)


# Fourier (mejor opción encontrada)
nbasis_fourier <- 20
basis_fourier <- create.fourier.basis(rangeval, nbasis_fourier)
fdPar_f <- fdPar(basis_fourier, lambda = 1e-10)
fit_fourier <- smooth.basis(X, Y_noisy, fdPar_f)

# B-spline (comparativa)
nbasis_bspline <- 70  # puedes ajustar este valor también
basis_bspline <- create.bspline.basis(rangeval, nbasis_bspline)
fdPar_b <- fdPar(basis_bspline, lambda = 1e-6)
fit_bspline <- smooth.basis(X, Y_noisy, fdPar_b)

# Funciones ajustadas
Y_fourier <- eval.fd(X, fit_fourier$fd)
Y_bspline <- eval.fd(X, fit_bspline$fd)

# Error cuadrático medio
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
cat("RMSE Fourier:", rmse(Y, Y_fourier), "\n")
cat("RMSE B-spline:", rmse(Y, Y_bspline), "\n")

# Funciones base de Fourier
matplot(eval.basis(X, basis_fourier), type = 'l', main = "Funciones base Fourier", col = rainbow(nbasis_fourier), lty = 1)

# Funciones base B-spline
matplot(eval.basis(X, basis_bspline), type = 'l', main = "Funciones base B-spline", col = rainbow(nbasis_bspline), lty = 1)

image(eval.penalty(basis_fourier, 2), main = "Penalización Fourier")
image(eval.penalty(basis_bspline, 2), main = "Penalización B-spline")

plot(X, Y, type = 'l', col = 'black', lwd = 2, ylim = range(c(Y, Y_fourier, Y_bspline)), main = "Comparación reconstrucciones")
lines(X, Y_fourier, col = 'blue', lwd = 2, lty = 2)
lines(X, Y_bspline, col = 'red', lwd = 2, lty = 3)
legend("topright", legend = c("Original", "Fourier", "B-spline"),
       col = c("black", "blue", "red"), lty = c(1,2,3), lwd = 2)
# ===================================================



