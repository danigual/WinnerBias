# R/sim_genetics.R

#' Generate a virtual population for GWAS simulation (Robust & Controlled)
#'
#' @param n_samples Sample size (number of individuals)
#' @param n_snps Total number of SNPs
#' @param n_causal Number of causal SNPs
#' @param h2 Target heritability (0-1)
#' @param mean_beta Standard deviation for causal effects (Effect Magnitude)
#' @param seed seed for reproducibility
#'
#' @return List with X (genotypes), y (phenotype), true_betas, and metadata.
#' @export
generate_gwas_data <- function(n_samples = 50000, 
                               n_snps = 1000, 
                               n_causal = 100, 
                               h2 = 0.5,
                               mean_beta = 0.5,
                               seed = NULL) {     
  
  # --- 1. SETUP AND CHECKS ---
  
  if (!is.null(seed)) set.seed(seed)
  
  # Sanity Checks
  stopifnot(n_causal <= n_snps)
  stopifnot(h2 >= 0 && h2 <= 1)
  stopifnot(n_samples > 10)
  
  # --- 2. GENOTYPES (X) - CONTROLLED---
  
  # We maintain a MAF = prob = 0.3 fixed (without variable MAF) for experimental control.
  X <- matrix(rbinom(n_samples * n_snps, size = 2, prob = 0.3), 
              nrow = n_samples, ncol = n_snps)
  
  # Column names
  colnames(X) <- paste0("SNP_", 1:n_snps)
  
  # --- 3. REAL EFFECTS (Biology) ---
  
  true_betas <- rep(0, n_snps)
  causal_indices <- sample(1:n_snps, n_causal)
  
  # We generate effects centred on 0 with mean_beta dispersion.
  true_betas[causal_indices] <- rnorm(n_causal, mean = 0, sd = mean_beta)
  
  # --- 4. GENETIC COMPONENT (G)  ---
  
  g_component <- as.vector(X %*% true_betas)
  var_g <- var(g_component)
  
  # Protection against zero variance (if n_causal is very low or bad luck)
  if (var_g < 1e-9) var_g <- 1 
  
  # --- 5. ENVIRONMENTAL COMPONENT (E) ---
  
  # We calculate the exact noise required for heritability (h2) = 0.5
  # Ve = Vg * (1 - h2) / h2
  var_e_target <- var_g * (1 - h2) / h2
  sigma_e <- sqrt(var_e_target)
  
  e_component <- rnorm(n_samples, mean = 0, sd = sigma_e)
  
  # --- 6. PHENOTYPE (Y) ---
  
  y <- g_component + e_component
  
  # --- 7. DATA AND OUTPUT ---
  
  h2_empirico <- var(g_component) / var(y)
  
  return(list(
    X = X,
    y = y,
    true_betas = true_betas,
    causal_indices = causal_indices,
    meta = list(
      h2_target = h2,
      h2_real_empirico = h2_empirico,
      var_g = var_g,
      var_e = var(e_component)
    )
  ))
}