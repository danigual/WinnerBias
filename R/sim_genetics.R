
# R/sim_genetics.R

#' Generate a virtual population for GWAS simulation
#'
#' @param n_samples sample size (numbrer of individuals)
#' @param n_snps number of SNPs
#' @param n_causal number of SNPS with true effect (causal SNPs)
#' @param h2 heritability, indicates the genetic background (0-1)
#'
#' @return  List with matrix of genotypes (X), 
#' phenotypes (Y), and true effects (true_betas)
#' @export
generate_gwas_data <- function(n_samples = 1000, 
                               n_snps = 10000, 
                               n_causal = 50, 
                               h2 = 0.5) {
  
  # 1. Genotypes (0, 1, 2)
  # MAF (minor allele frequency) =  prob = 0.3 
  X <- matrix(rbinom(n_samples * n_snps, size = 2, prob = 0.3), 
              nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("SNP_", 1:n_snps)
  
  # 2. Real Effects (Biology) 
  true_betas <- rep(0, n_snps)
  causal_indices <- sample(1:n_snps, n_causal)
  true_betas[causal_indices] <- rnorm(n_causal, mean = 0, sd = 1)
  
  # 3. Genetic Component (G)
  g_component <- as.vector(X %*% true_betas)
  
  # 4. Environmental Component (E) - Adjusted for heritability
  var_g <- var(g_component)
  if(var_g == 0) var_g <- 1
  
  sigma_e <- sqrt(var_g * (1 - h2) / h2)
  e_component <- rnorm(n_samples, mean = 0, sd = sigma_e)
  
  # 5. Phenotype (Y)
  y <- g_component + e_component
  
  return(list(
    X = X,
    y = y,
    true_betas = true_betas,
    causal_indices = causal_indices
  ))
}