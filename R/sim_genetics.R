# R/sim_genetics.R

#' Generate a virtual population for GWAS simulation (Robust & Controlled)
#'
#' @param n_samples Sample size (number of individuals)
#' @param n_snps Total number of SNPs
#' @param n_causal Number of causal SNPs
#' @param h2 Target heritability (0-1)
#' @param mean_beta Standard deviation for causal effects (Effect Magnitude)
#' @param seed Optional seed for reproducibility
#'
#' @return List with X (genotypes), y (phenotype), true_betas, and metadata.
#' @export
generate_gwas_data <- function(n_samples = 1000, 
                               n_snps = 10000, 
                               n_causal = 50, 
                               h2 = 0.5,
                               mean_beta = 0.125, # Parametrizado (Tu valor de Winner's Curse)
                               seed = NULL) {     # [MEJORA 1] Semilla Opcional
  
  # --- 1. CONFIGURACIÓN Y CHEQUEOS ---
  if (!is.null(seed)) set.seed(seed)
  
  # [MEJORA 2] Chequeos de seguridad (Sanity Checks)
  stopifnot(n_causal <= n_snps)
  stopifnot(h2 >= 0 && h2 <= 1)
  stopifnot(n_samples > 10)
  
  # --- 2. GENOTIPOS (X) - CONTROLADO ---
  # Mantenemos prob = 0.3 fijo (Sin MAF variable) para control experimental.
  X <- matrix(rbinom(n_samples * n_snps, size = 2, prob = 0.3), 
              nrow = n_samples, ncol = n_snps)
  
  # [MEJORA 4] Nombres de columnas
  colnames(X) <- paste0("SNP_", 1:n_snps)
  
  # --- 3. EFECTOS REALES (Biología) ---
  true_betas <- rep(0, n_snps)
  causal_indices <- sample(1:n_snps, n_causal)
  
  # Generamos efectos centrados en 0 con la dispersión que tú definas
  true_betas[causal_indices] <- rnorm(n_causal, mean = 0, sd = mean_beta)
  
  # --- 4. COMPONENTE GENÉTICO (G) ---
  g_component <- as.vector(X %*% true_betas)
  var_g <- var(g_component)
  
  # Protección contra varianza cero (si n_causal es muy bajo o mala suerte)
  if (var_g < 1e-9) var_g <- 1 
  
  # --- 5. COMPONENTE AMBIENTAL (E) ---
  # Calculamos el ruido exacto necesario para clavar la heredabilidad (h2)
  # Ve = Vg * (1 - h2) / h2
  var_e_target <- var_g * (1 - h2) / h2
  sigma_e <- sqrt(var_e_target)
  
  e_component <- rnorm(n_samples, mean = 0, sd = sigma_e)
  
  # --- 6. FENOTIPO (Y) ---
  y <- g_component + e_component
  
  # --- 7. METADATOS Y SALIDA ---
  # [MEJORA 3] Devolvemos metadatos para verificar que todo ha ido bien
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