# R/bootstrap_correction.R

#' Winner's Curse Correction via Bootstrap Aggregating (Bagging)
#'
#' Implements non-parametric bootstrap to correct upward bias.
#' Includes calculation of Standard Errors (SE) for the corrected estimates.
#'
#' @param winners_df Dataframe with 'snp_id' column.
#' @param X_genotype Matrix (N x M).
#' @param y_phenotype Vector (N).
#' @param n_boot Iterations (Default 200 for stability).
#' @param seed Random seed.
#'
#' @return winners_df with 'beta_corrected' and 'beta_boot_se'
#' @export
run_bootstrap_correction <- function(winners_df, X_genotype, y_phenotype, n_boot = 200, seed = 123) {
  
  # --- 1. Input Validation ---
  if(nrow(winners_df) == 0) {
    warning("No significant variants provided for correction.")
    return(winners_df)
  }
  
  # Check consistency between IDs and Matrix Columns
  snp_ids <- winners_df$snp_id
  missing_snps <- snp_ids[!snp_ids %in% colnames(X_genotype)]
  
  if(length(missing_snps) > 0) {
    stop(paste0("ERROR DE DATOS: Los siguientes SNPs del dataframe no existen en la matriz X:\n", 
                paste(head(missing_snps, 3), collapse=", "), "...\n",
                "Verifica si los nombres tienen formato 'SNP_x' o 'rs...' en ambos lados."))
  }
  
  # --- 2. Setup ---
  message("\n--- Starting Bootstrap Correction (Bagging) ---")
  message(paste(" Target Variants      :", nrow(winners_df)))
  message(paste(" Bootstrap Iterations :", n_boot))
  
  set.seed(seed)
  N <- nrow(X_genotype)
  
  # OPTIMIZATION: Subset Genotype Matrix
  # Only keep the winners to speed up processing
  X_subset <- X_genotype[, snp_ids, drop = FALSE]
  
  # Matrix to store bootstrap estimates [Rows = SNPs, Cols = Bootstraps]
  boot_betas <- matrix(NA, nrow = nrow(winners_df), ncol = n_boot)
  
  # Progress Bar
  pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  
  # --- 3. Bootstrap Loop ---
  for(b in 1:n_boot) {
    
    # A. Resampling (Sampling with replacement)
    idx_boot <- sample(1:N, N, replace = TRUE)
    
    X_b <- X_subset[idx_boot, , drop = FALSE]
    y_b <- y_phenotype[idx_boot]
    
    # B. Fast Vectorized Beta Calculation
    # Beta = Cov(x,y) / Var(x)
    # Centering allows us to ignore the intercept
    y_center <- y_b - mean(y_b)
    X_center <- scale(X_b, center = TRUE, scale = FALSE) 
    
    numerator   <- t(X_center) %*% y_center
    denominator <- colSums(X_center^2)
    
    # SAFETY: Prevent division by zero (monomorphic SNPs in bootstrap sample)
    denominator[denominator < 1e-9] <- NA
    
    # Store result
    boot_betas[, b] <- as.vector(numerator / denominator)
    
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  # --- 4. Aggregation & Metrics ---
  
  # A. Corrected Estimate (Bagging using Median for robustness)
  winners_df$beta_corrected <- apply(boot_betas, 1, median, na.rm = TRUE)
  
  # B. Bootstrap Standard Error (SE)
  # Nos dice cuán estable es nuestra corrección
  winners_df$beta_corrected_se <- apply(boot_betas, 1, sd, na.rm = TRUE)
  
  message("\nCorrection complete. Added 'beta_corrected' and 'beta_corrected_se'.")
  return(winners_df)
}