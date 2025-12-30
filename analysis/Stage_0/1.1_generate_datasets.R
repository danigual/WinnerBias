# analysis/Stage_0/1.1_generate_datasets.R

# 1. load the generator function
source("R/sim_genetics.R")

# 2. define the parameters of the control scenario
PARAMS <- list(
  n_samples = 50000,  # sample size
  n_snps    = 10000,  # number of genetic variants (SNPs)
  n_causal  = 50,     # number of significant SNPs (causal)
  h2        = 0.5,    # heritability (50% genetic, 50% environmental)
  seed      = 42      
)

# 3. Run the simulation
message("Generando datos simulados (Población Virtual)...")
set.seed(PARAMS$seed)

sim_data <- generate_gwas_data(
  n_samples = PARAMS$n_samples,
  n_snps    = PARAMS$n_snps,
  n_causal  = PARAMS$n_causal,
  h2        = PARAMS$h2
)

# 4. Save data
# we create the folder if it does not exist (for security reasons).
dir_stage1 <- "data/processed/Stage_0"
if(!dir.exists(dir_stage1)) dir.create(dir_stage1, recursive = TRUE)

output_file <- file.path(dir_stage1, "simulacion_wc.rds")

saveRDS(sim_data, file = output_file)

message("Datos generados y guardados en: ", output_file)