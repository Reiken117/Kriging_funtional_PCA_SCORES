## ===============================================================
## Kriging Funcional Completo según Giraldo et al. (2011)
## ===============================================================

# 1. Cargar librerías necesarias -----------------------------------------------
library(sf)
library(sp)
library(fda)
library(gstat)
library(ggplot2)
library(dplyr)
library(viridis)
library(purrr)

# 2. Leer y preparar datos ------------------------------------------------------
gdf <- st_read("G:/Univalle/puntos_generados_espo_9377.geojson")
if(nrow(gdf)==0) stop("No hay datos geoposicionados")
# Extraer fecha y coordenadas
gdf_df <- gdf %>%
  mutate(DATE = as.Date(DATE)) %>%
  mutate(
    X_coord = st_coordinates(geometry)[,1],
    Y_coord = st_coordinates(geometry)[,2]
  ) %>%
  st_drop_geometry()

gdf_df <- gdf_df %>%
  group_by(DATE) %>%
  filter(n() >= 2) %>%
  arrange(DATE, X_coord)

X_common <- seq(min(gdf_df$X_coord), max(gdf_df$X_coord), length.out = 100)
Y_matrix <- gdf_df %>%
  group_by(DATE) %>%
  group_map(~ approx(.x$X_coord, .x$Y_coord, xout = X_common, rule = 2)$y) %>%
  discard(~ all(is.na(.x))) %>%
  do.call(cbind, .)
coords <- gdf_df %>%
  distinct(DATE, .keep_all = TRUE) %>%
  summarise(X = mean(X_coord), Y = mean(Y_coord)) %>%
  select(X, Y) %>% as.matrix()
n_sites <- ncol(Y_matrix)
X        <- seq(0, 1, length.out = nrow(Y_matrix))
Y        <- Y_matrix

# --- Función auxiliar para kriging interno en CV ------------------------------
krige_inner <- function(pt, coords_cv, z_cv, basis_cv, vgm_model_cv) {
  n_cv <- ncol(z_cv)
  # dist from pt to each remaining site
  h0 <- sqrt(colSums((t(coords_cv) - pt)^2))
  gamma0 <- variogramLine(vgm_model_cv, dist_vector = h0)$gamma
  # semivariogram between sites
  d0 <- as.matrix(dist(coords_cv))
  Gij <- outer(1:n_cv, 1:n_cv, Vectorize(
    function(i,j) variogramLine(vgm_model_cv, dist_vector = d0[i,j])$gamma
  ))
  A <- rbind(cbind(Gij, 1), c(rep(1, n_cv), 0))
  b <- c(gamma0, 1)
  sol <- qr.solve(A + diag(1e-6, nrow(A)), b)
  lam <- sol[1:n_cv]
  pred_coef <- z_cv %*% lam
  eval.fd(X, fd(pred_coef, basis_cv))
}
# 3. Functional CV para elegir L y eta ------------------------------------------
L_vec   <- seq(15, 20)
eta_vec <- 10^seq(-6, -2, length = 5)
results_fcv <- expand.grid(L = L_vec, eta = eta_vec)
results_fcv$FCV <- NA

for(r in seq_len(nrow(results_fcv))) {
  L   <- results_fcv$L[r]
  eta <- results_fcv$eta[r]
  
  # Base B-spline y penalización
  basis_cv <- create.bspline.basis(range(X), nbasis = L + 4)
  fdPar_cv <- fdPar(basis_cv, lambda = eta)
  W_cv     <- inprod(basis_cv, basis_cv)
  SSE_all  <- 0
  
  for(i in seq_len(n_sites)) {
    # 1) Leave-one-site-out
    coords_cv <- coords[-i, ]
    Y_cv      <- Y[, -i]
    
    # 2) Suavizado funcional
    sm   <- smooth.basis(X, Y_cv, fdPar_cv)
    z_cv <- sm$fd$coefs
    
    # 3) Cálculo del trace-variograma empírico
    d_mat_cv  <- as.matrix(dist(coords_cv))
    combs     <- combn(ncol(z_cv), 2)
    h_vals    <- apply(combs, 2, function(idx) d_mat_cv[idx[1], idx[2]])
    trace_vals <- apply(combs, 2, function(idx) {
      dc <- z_cv[, idx[1]] - z_cv[, idx[2]]
      as.numeric(t(dc) %*% W_cv %*% dc)
    })
    
    # 4) Binning y cálculo de semivarianzas
    bins <- cut(h_vals, breaks = 8, include.lowest = TRUE)
    gamma_emp <- tapply(trace_vals, bins, mean, na.rm = TRUE) / 2
    np_emp    <- tapply(trace_vals, bins, function(x) sum(!is.na(x)))
    h_mid <- sapply(levels(bins), function(lvl) {
      nums <- as.numeric(unlist(regmatches(lvl, gregexpr("[0-9.]+", lvl))))
      if(length(nums) == 2) mean(nums) else NA
    })
    
    # 5) Preparar datos para ajuste del variograma
    df_var <- data.frame(
      np    = as.numeric(np_emp),
      dist  = h_mid,
      gamma = gamma_emp,
      dir.hor = 0, dir.ver = 0, id = "var1"
    ) %>% filter(is.finite(dist), is.finite(gamma))
    class(df_var) <- c("gstatVariogram", "data.frame")
    
    # 6) Ajuste robusto del variograma (según artículo)
    if(nrow(df_var) >= 4 && max(df_var$dist, na.rm = TRUE) > 0) {
      # Modelos sugeridos en el artículo + manejo de errores
      vgm_models <- list(
        vgm("Exp", nugget = 0.5*max(df_var$gamma), 
            psill = 0.5*max(df_var$gamma), 
            range = max(df_var$dist)/2),
        vgm("Sph", nugget = 0.5*max(df_var$gamma), 
            psill = 0.5*max(df_var$gamma), 
            range = max(df_var$dist)/3),
        vgm("Gau", nugget = 0.5*max(df_var$gamma), 
            psill = 0.5*max(df_var$gamma), 
            range = max(df_var$dist)/2)
      )
      
      best_model <- NULL
      min_SS <- Inf
      
      # Probar múltiples modelos y seleccionar el mejor
      for(model in vgm_models) {
        fit <- tryCatch(
          fit.variogram(df_var, model, 
                        fit.method = 6,  # WLS con pesos
                        weights = df_var$np),
          error = function(e) NULL
        )
        if(!is.null(fit)) {  # Corregido: cerrar paréntesis
          current_SS <- attr(fit, "SSErr")
          if(current_SS < min_SS) {
            min_SS <- current_SS
            best_model <- fit
          }
        }
      }
      
      # Fallback a modelo nugget si todos fallan
      if(is.null(best_model)) {
        best_model <- vgm("Nug", psill = mean(df_var$gamma, na.rm = TRUE))
      }
      
      # 7) Predicción y cálculo de SSE
      pred_i <- tryCatch({
        krige_inner(coords[i, ], coords_cv, z_cv, basis_cv, best_model)
      }, error = function(e) {
        rep(NA, nrow(Y))  # Manejo de errores en predicción
      })
      
      SSE_all <- SSE_all + sum((Y[, i] - pred_i)^2, na.rm = TRUE)
      
    } else {
      SSE_all <- Inf
      break
    }
  }
  
  results_fcv$FCV[r] <- SSE_all
  cat("Progreso FCV:", r, "/", nrow(results_fcv), "| L:", L, "eta:", eta, "FCV:", SSE_all, "\n")
}

# Selección de parámetros óptimos
opt <- results_fcv[which.min(results_fcv$FCV), ]
L_opt  <- opt$L
eta_opt <- opt$eta
cat("\nFCV óptimo: L =", L_opt, "eta =", eta_opt, "\n")

# 4. Suavizar con parámetros óptimos -------------------------------------------
basis_opt <- create.bspline.basis(range(X), nbasis = L_opt + 4)
fdPar_opt <- fdPar(basis_opt, lambda = eta_opt)
sm_opt <- smooth.basis(X, Y, fdPar_opt)
z_full <- sm_opt$fd$coefs

# 5. Estimación empírica del trace-variogram ------------------------------------
W_opt <- inprod(basis_opt, basis_opt)
d_mat <- as.matrix(dist(coords))
combs <- combn(n_sites,2)
h_vals <- apply(combs, 2, function(idx) d_mat[idx[1], idx[2]])
trace_vals <- apply(combs, 2, function(idx) {
  dc <- z_full[,idx[1]] - z_full[,idx[2]]
  as.numeric(t(dc) %*% W_opt %*% dc)
})
# binning dinámico
n_pairs <- length(h_vals)
n_bins  <- min(5, max(4, floor(n_pairs / 10)))
bins    <- cut(h_vals, breaks = n_bins, include.lowest = TRUE)
gamma_emp <- tapply(trace_vals, bins, mean, na.rm=TRUE)/2
np_emp    <- tapply(trace_vals, bins, function(x) sum(!is.na(x)))
h_mid     <- sapply(levels(bins), function(lvl) {
  nums <- as.numeric(unlist(regmatches(lvl, gregexpr("[0-9\\.]+", lvl))))
  if(length(nums) == 2) mean(nums) else NA
})

vgm_emp <- data.frame(np = as.numeric(np_emp), dist = h_mid,
                      gamma = gamma_emp, dir.hor=0, dir.ver=0, id="var1")
class(vgm_emp) <- c("gstatVariogram","data.frame")
vgm_emp <- vgm_emp %>% filter(is.finite(dist), is.finite(gamma))

# 6. Ajuste WLS del variograma (fit.method = 2 ya usa np como pesos por defecto) ----

# 6.1) Construir modelo inicial con parámetros de partida
vgm0 <- vgm("Sph",
            nugget = 0,
            psill  = max(vgm_emp$gamma, na.rm = TRUE),
            range  = max(vgm_emp$dist,  na.rm = TRUE) / 2)

# 6.2) Ajuste por Weighted Least Squares (fit.method = 2)
vgm_model <- tryCatch(
  fit.variogram(vgm_emp,
                model      = vgm0,
                fit.method = 2),
  error = function(e) NULL
)

# 6.3) Fallback a Cressie (fit.method = 6) si WLS falla
if (is.null(vgm_model)) {
  vgm_model <- fit.variogram(vgm_emp,
                             model      = vgm0,
                             fit.method = 6)
}

# 6.4) (Opcional) corregir rango negativo si ocurriera
if (!is.null(vgm_model) && vgm_model$range[2] < 0) {
  vgm_model$range[2] <- abs(vgm_model$range[2])
}

# 7. Predicción en un sitio de ejemplo -----------------------------------------
s0 <- coords[1, ]
pred0 <- krige_inner(s0, coords, z_full, basis_opt, vgm_model)
plot(X_common, pred0, type='l', lwd=2,
     xlab='X espacial (normalizado)', ylab='Valor predicho',
     main='Kriging Funcional en sitio s0')

# 8. Validación LOOCV del kriging ----------------------------------------------
pred_loocv <- matrix(NA, nrow=length(X), ncol=n_sites)
for(i in seq_len(n_sites)) {
  pred_loocv[,i] <- krige_inner(coords[i, ], coords, z_full, basis_opt, vgm_model)
}
rmse_loocv <- sqrt(mean((Y - pred_loocv)^2, na.rm=TRUE))
cat("RMSE LOOCV:", rmse_loocv, "\n")

# 9. Visualización final ------------------------------------------------------
matplot(X_common, Y, type='l', col='gray80', lty=1,
        xlab='X normalizado', ylab='Y')
matlines(X_common, pred_loocv, col='red', lty=1)
legend('topright', legend=c('Observado','LOOCV'), col=c('gray','red'), lty=1)


# Visualización
plot(
  vgm_emp$dist, vgm_emp$gamma, type = "b", pch = 16,
  xlab = "Distancia", ylab = "Trace-variograma",
  main = "Variograma Empírico (trace)"
)
plot(vgm_emp, vgm_model, main = "Ajuste del variograma mejor modelo")


fechas <- unique(gdf_df$DATE)
fechas[1]


# ————— 8.bis. RMSE por fecha y horizonte temporal —————

# Partimos de que ya tienes:
# fechas  = vector de fechas únicas (class Date) de longitud n_sites
# Y       = matriz de datos original (dimensiones: longitud X_common × n_sites)
# coords  = matriz de coordenadas (n_sites × 2)
# basis_opt, z_full, vgm_model = objetos que usas en krige_inner

# 1) Calcular RMSE leave‑one‑out por fecha
rmse_per_date <- numeric(n_sites)
for(i in seq_len(n_sites)) {
  # predecir la curva de la fecha i dejando fuera el sitio i
  pred_i <- krige_inner(coords[i, ], coords[-i, ], z_full[, -i], basis_opt, vgm_model)
  # RMSE en la curva
  rmse_per_date[i] <- sqrt(mean((Y[, i] - pred_i)^2, na.rm = TRUE))
}

# 2) Construir data.frame con horizonte temporal (días desde la primera fecha)
fechas_sorted <- sort(unique(gdf_df$DATE))
horizontes   <- as.numeric(fechas_sorted - min(fechas_sorted))
df_error <- data.frame(
  fecha    = fechas_sorted,
  horizonte = horizontes,
  RMSE     = rmse_per_date[order(unique(gdf_df$DATE))]
)

# 3) Gráfico RMSE vs horizonte
library(ggplot2)
ggplot(df_error, aes(x = horizonte, y = RMSE)) +
  geom_line() +
  geom_point() +
  labs(
    x     = "Horizonte temporal (días desde la 1ª muestra)",
    y     = "RMSE de predicción (misma unidad que Y)",
    title = "Error de predicción vs horizonte temporal"
  )

