# R/bootstrap_correction.R

#' Winner's Curse Correction via Bootstrap Aggregating (Bagging)
#'
#' This function implements a non-parametric bootstrap approach to correct
#' the upward bias (Winner's Curse) in GWAS effect size estimates.
#' It resamples the original dataset with replacement and re-estimates the
#' effect sizes for the selected significant variants.
#'
#' @param winners_df A dataframe containing the significant GWAS hits. 
#'        Must include a column 'snp_id' matching the columns in X_genotype.
#' @param X_genotype The genotype matrix (N individuals x M SNPs).
#' @param y_phenotype The phenotype vector (N individuals).
#' @param n_boot Integer. Number of bootstrap iterations. Default is 100.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return The original 'winners_df' with an additional column 'beta_corrected',
#'         containing the average effect size across bootstrap iterations.
#'
#' @export

run_bootstrap_correction <- function(winners_df, X_genotype, y_phenotype, n_boot = 100, seed = 123) {
  
  # 1. Input Validation
  if(nrow(winners_df) == 0) {
    warning("No significant variants provided for correction. Returning original dataframe.")
    return(winners_df)
  }
  
  # 2. Setup
  message("\n--- Starting Bootstrap Correction ---")
  message("Target Variants: ", nrow(winners_df))
  message("Bootstrap Iterations: ", n_boot)
  
  set.seed(seed)
  N <- nrow(X_genotype)
  
  # OPTIMIZATION: Subset Genotype Matrix
  # We only need the columns corresponding to the "winner" SNPs.
  # This drastically reduces computational time (e.g., from 10,000 columns to ~50).
  snp_ids <- winners_df$snp_id
  X_subset <- X_genotype[, snp_ids, drop = FALSE]
  
  # Initialize matrix to store bootstrap estimates
  # Rows: SNPs, Columns: Iterations
  boot_betas <- matrix(NA, nrow = nrow(winners_df), ncol = n_boot)
  
  # Progress bar to keep the user informed
  pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  
  # 3. Bootstrap Loop
  for(b in 1:n_boot) {
    
    # A. Resampling (The core of the Bootstrap)
    # Sample N indices with replacement
    idx_boot <- sample(1:N, N, replace = TRUE)
    
    # Create the bootstrap sample
    X_b <- X_subset[idx_boot, , drop = FALSE]
    y_b <- y_phenotype[idx_boot]
    
    # B. Fast Vectorized Beta Calculation
    # We use linear algebra directly instead of lm() for speed.
    # Formula: Beta = Cov(x,y) / Var(x)
    
    # Center variables to simplify calculation (remove intercept)
    y_center <- y_b - mean(y_b)
    X_center <- scale(X_b, center = TRUE, scale = FALSE)
    
    # Calculate Beta: (X'y) / (X'X)
    numerator <- t(X_center) %*% y_center
    denominator <- colSums(X_center^2)
    
    # Store result for this iteration
    boot_betas[, b] <- as.vector(numerator / denominator)
    
    # Update progress bar
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  # 4. Aggregation (Bagging)
  # The corrected estimate is the mean of the bootstrap estimates.
  # This works because the bootstrap mean tends to be closer to the true population mean
  # than the selected extreme value from the discovery sample.
  winners_df$beta_corrected <- rowMeans(boot_betas, na.rm = TRUE)
  
  message("Correction complete.")
  return(winners_df)
}