#' Ejecutar GWAS (Genome-Wide Association Study)
#'
#' Calcula beta, error estándar y p-valor de manera vectorizada y segura.
#' Optimizado para memoria (evita scale() masivo) y robusto ante SNPs constantes.
#'
#' @param X Matriz de genotipos (n_samples x n_snps).
#' @param y Vector de fenotipos continuo.
#'
#' @return Dataframe con columnas: snp_id, beta_hat, se, pval.
#' @export
run_gwas <- function(X, y) {
  
  # --- 1. COMPROBACIONES DE ENTRADA (Safety Checks) ---
  if (nrow(X) != length(y)) {
    stop(paste("Error de dimensiones: X tiene", nrow(X), "filas pero y tiene", length(y)))
  }
  
  n_samples <- nrow(X)
  n_snps    <- ncol(X)
  
  # Control de grados de libertad (df)
  # Necesitamos al menos 3 muestras para calcular varianza y tener df > 0
  if (n_samples <= 2) {
    stop("Error: Tamaño de muestra insuficiente (N <= 2). El modelo no es identificable.")
  }
  
  # --- 2. PREPARACIÓN EFICIENTE (Sin scale() para ahorrar RAM) ---
  # En lugar de crear una matriz centrada gigante (X_centered), 
  # usamos las propiedades de la covarianza: Cov(X,Y) = E[XY] - E[X]E[Y]
  
  # Estadísticos de Y
  mean_y <- mean(y)
  sst_y  <- sum((y - mean_y)^2)  # Suma Total de Cuadrados de Y
  
  # Estadísticos de X (Vectorizados)
  # colSums es muy rápido y no copia la matriz
  sum_x  <- colSums(X)
  mean_x <- sum_x / n_samples
  
  # Suma de cuadrados de X (Denominador de la varianza)
  # SX2 = sum(x^2) - n * mean_x^2
  sum_sq_x <- colSums(X^2)
  ssx_denom <- sum_sq_x - n_samples * (mean_x^2)
  
  # --- 3. FILTRO DE MONOMÓRFICOS (Varianza 0) ---
  # Si ssx_denom es casi 0, el SNP es constante. Lo marcamos para evitar NaNs.
  # Usamos una tolerancia pequeña para flotantes.
  is_polymorphic <- ssx_denom > 1e-8
  
  # --- 4. CÁLCULO DE BETAS (Solo en polimórficos) ---
  # Numerador = sum(xy) - n * mean_x * mean_y
  # crossprod es la forma más eficiente de hacer t(X) %*% y
  prod_xy <- as.vector(crossprod(X, y))
  numerador <- prod_xy - n_samples * mean_x * mean_y
  
  # Inicializamos vectores de resultados
  betas <- rep(NA, n_snps)
  se    <- rep(NA, n_snps)
  pvals <- rep(NA, n_snps)
  
  # Calculamos SOLO para los válidos (evitamos división por cero)
  if (any(is_polymorphic)) {
    
    # 4.1 Betas
    betas[is_polymorphic] <- numerador[is_polymorphic] / ssx_denom[is_polymorphic]
    
    # 4.2 Error Estándar (Vectorizado Eficiente)
    # SSR (Explicada) = beta^2 * SSX
    ssr <- (betas[is_polymorphic]^2) * ssx_denom[is_polymorphic]
    
    # SSE (Residual) = SST - SSR
    sse <- sst_y - ssr
    
    # Protección numérica: SSE no puede ser negativo
    sse[sse < 0] <- 0
    
    # Grados de libertad: N - 2 (Intercepto + Pendiente)
    df <- n_samples - 2
    
    # Varianza residual (sigma^2) y SE
    var_residual <- sse / df
    se[is_polymorphic] <- sqrt(var_residual / ssx_denom[is_polymorphic])
    
    # 4.3 P-valores
    t_stats <- betas[is_polymorphic] / se[is_polymorphic]
    pvals[is_polymorphic] <- 2 * pt(-abs(t_stats), df)
  }
  
  # --- 5. EMPAQUETADO ---
  # Añadimos advertencia si hubo monomórficos descartados
  n_drop <- n_snps - sum(is_polymorphic)
  if (n_drop > 0) {
    warning(paste("Se han ignorado", n_drop, "SNPs monomórficos (Varianza 0)."))
  }
  
  res <- data.frame(
    snp_id   = colnames(X),
    beta_hat = betas,
    se       = se,
    pval     = pvals
  )
  
  return(res)
}