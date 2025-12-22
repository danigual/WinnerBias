# R/02_gwas_engine.R

#' Ejecutar GWAS (Genome-Wide Association Study)
#'
#' Calcula beta, error estándar y p-valor para cada SNP contra el fenotipo
#' de manera vectorizada (rápida).
#'
#' @param X Matriz de genotipos (n_samples x n_snps).
#' @param y Vector de fenotipos continuo.
#'
#' @return Dataframe con columnas: snp_id, beta_hat, se, pval.
#' @export
run_gwas <- function(X, y) {
  
  # Chequeo de dimensiones
  n_samples <- nrow(X)
  n_snps <- ncol(X)
  
  # 1. Preparación Algebraica
  # Centramos para simplificar cálculos de covarianza
  y_centered <- y - mean(y)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  
  # 2. Calcular Betas (Estimación del efecto)
  # Fórmula OLS rápida: (X'y) / (X'X)
  numerador <- t(X_centered) %*% y_centered
  denominador <- colSums(X_centered^2)
  
  betas <- numerador / denominador
  
  # 3. Calcular Error Estándar (SE)
  # Necesitamos la varianza residual para cada SNP
  se <- numeric(n_snps)
  df <- n_samples - 2 # Grados de libertad (N - 2 parámetros)
  
  # Pre-calculamos varianza total de Y
  sst <- sum(y_centered^2)
  
  for(i in 1:n_snps) {
    b <- betas[i]
    
    # Suma de cuadrados de la regresión y del error
    ssr <- b^2 * denominador[i]
    sse <- sst - ssr
    
    # Evitar errores numéricos si sse < 0
    if(sse < 0) sse <- 0
    
    sigma2 <- sse / df
    se[i] <- sqrt(sigma2 / denominador[i])
  }
  
  # 4. Calcular P-valores
  t_stat <- betas / se
  pvals <- 2 * pt(-abs(t_stat), df)
  
  # 5. Empaquetar resultados
  res <- data.frame(
    snp_id = colnames(X),
    beta_hat = as.vector(betas),
    se = se,
    pval = pvals
  )
  
  return(res)
}