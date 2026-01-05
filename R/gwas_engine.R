#' Run GWAS (Genome-wide association study)
#' 
#' Calculates beta, standard error, and p-value in a vectorised and secure manner.
#' Optimised for memory (avoids massive scale) and robust against constant SNPs.
#'
#' @param X Genotype matrix (n_samples x n_snps).
#' @param y Continuous phenotype vector.
#'
#' @return Dataframe with columns: snp_id, beta_hat, se, pval.
#' @export
run_gwas <- function(X, y) {
  
  if (nrow(X) != length(y)) {
    stop(paste("Error de dimensiones: X tiene", nrow(X), "filas pero y tiene", length(y)))
  }
  
  n_samples <- nrow(X)
  n_snps    <- ncol(X)
  
  # --- 1. INPUT CHECKS (Safety Checks) ---
  
  # Control of degrees of freedom (df)
  # We need at least 3 samples to calculate variance and have df > 0
  if (n_samples <= 2) {
    stop("Error: Tamaño de muestra insuficiente (N <= 2). El modelo no es identificable.")
  }
  
  # --- 2. EFFICIENT PREPARATION (Without scale() to save RAM) ---
  
  # Instead of creating a giant centred matrix (X_centred), 
  # We use the properties of covariance: Cov(X,Y) = E[XY] - E[X]E[Y]
  
  # Y statistics
  mean_y <- mean(y)
  sst_y  <- sum((y - mean_y)^2)  # Sum of Squares of Y
  
  # X statistics (vectorised)
  sum_x  <- colSums(X)
  mean_x <- sum_x / n_samples
  
  # Sum of squares of X (Denominator of variance)
  # SX2 = sum(x^2) - n * mean_x^2
  sum_sq_x <- colSums(X^2)
  ssx_denom <- sum_sq_x - n_samples * (mean_x^2)
  
  # --- 3. MONOMORPHIC FILTER (Variance 0) ---
  
  # If ssx_denom is close to 0, the SNP is constant. We flag it to avoid divide by zero
  # We use a small tolerance for floats.
  is_polymorphic <- ssx_denom > 1e-8
  
  # --- 4. CALCULATION OF BETAS (Only in polymorphic cases) ---
  
  # Numerator = sum(xy) - n * mean_x * mean_y
  prod_xy <- as.vector(crossprod(X, y))
  numerador <- prod_xy - n_samples * mean_x * mean_y
  
  # Initialise result vectors
  betas <- rep(NA, n_snps)
  se    <- rep(NA, n_snps)
  pvals <- rep(NA, n_snps)
  
  # We calculate ONLY for valid values (avoiding division by zero).
  if (any(is_polymorphic)) {
    
    # 4.1 Betas
    betas[is_polymorphic] <- numerador[is_polymorphic] / ssx_denom[is_polymorphic]
    
    # 4.2 Standard Error (Efficient Vectorised)
    # SSR (Explained) = beta^2 * SSX
    ssr <- (betas[is_polymorphic]^2) * ssx_denom[is_polymorphic]
    
    # SSE (Residual) = SST - SSR
    sse <- sst_y - ssr
    
    # Numerical protection: SSE cannot be negative
    sse[sse < 0] <- 0
    
    # Degrees of freedom: N - 2 (Intercept + Slope)
    df <- n_samples - 2
    
    # Residual variance (sigma^2) and SE
    var_residual <- sse / df
    se[is_polymorphic] <- sqrt(var_residual / ssx_denom[is_polymorphic])
    
    # 4.3 P-values
    t_stats <- betas[is_polymorphic] / se[is_polymorphic]
    pvals[is_polymorphic] <- 2 * pt(-abs(t_stats), df)
  }
  
  # --- 5. PACKAGING ---
  
  # We add a warning if there were any monomorphic SNPs discarded.
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