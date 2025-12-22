
# R/sim_genetics.R

#' Generar una población virtual para simulación GWAS
#'
#' @param n_samples Número de individuos.
#' @param n_snps Número de variantes genéticas.
#' @param n_causal Número de variantes con efecto real.
#' @param h2 Heredabilidad (0-1).
#'
#' @return Lista con matriz de genotipos (X), fenotipo (y) y efectos reales (true_betas).
#' @export
generate_gwas_data <- function(n_samples = 1000, 
                               n_snps = 10000, 
                               n_causal = 50, 
                               h2 = 0.5) {
  
  # 1. Genotipos (0, 1, 2)
  # MAF = 0.3 (frecuencia del alelo menor)
  X <- matrix(rbinom(n_samples * n_snps, size = 2, prob = 0.3), 
              nrow = n_samples, ncol = n_snps)
  colnames(X) <- paste0("SNP_", 1:n_snps)
  
  # 2. Efectos Reales (Biología)
  true_betas <- rep(0, n_snps)
  causal_indices <- sample(1:n_snps, n_causal)
  true_betas[causal_indices] <- rnorm(n_causal, mean = 0, sd = 1)
  
  # 3. Componente Genético (G)
  g_component <- as.vector(X %*% true_betas)
  
  # 4. Componente Ambiental (E) - Ajustado por heredabilidad
  var_g <- var(g_component)
  if(var_g == 0) var_g <- 1
  
  sigma_e <- sqrt(var_g * (1 - h2) / h2)
  e_component <- rnorm(n_samples, mean = 0, sd = sigma_e)
  
  # 5. Fenotipo (Y)
  y <- g_component + e_component
  
  return(list(
    X = X,
    y = y,
    true_betas = true_betas,
    causal_indices = causal_indices
  ))
}