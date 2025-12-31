# analysis/Stage_2/1.1_generate_datasets.R

# 1. Load generator function
source("R/sim_genetics.R")

# 2. Define "Winner's Curse" scenario parameters
# ------------------------------------------------------------------------------
# Keeping N=1000 and Beta=0.125 to ensure low power and high inflation.

PARAMS <- list(
  n_samples = 1000,   # Number of patients (Low to provoke bias)
  n_snps    = 10000,  # Total variants
  n_causal  = 100,    # Real variants
  h2        = 0.5,    # Heritability
  mean_beta = 0.125,  # [IMPORTANT] Effect size (Winner's Curse target)
  seed      = 42      # Seed for full reproducibility
)

# 3. Run simulation
# ------------------------------------------------------------------------------
message(paste0("Generating simulated data (N=", PARAMS$n_samples, ", Beta=", PARAMS$mean_beta, ")..."))

sim_data <- generate_gwas_data(
  n_samples = PARAMS$n_samples,
  n_snps    = PARAMS$n_snps,
  n_causal  = PARAMS$n_causal,
  h2        = PARAMS$h2,
  mean_beta = PARAMS$mean_beta,
  seed      = PARAMS$seed
)

# [IMPROVEMENT 1] INJECT PARAMETERS INTO OBJECT
# This makes the .rds file self-descriptive.
sim_data$PARAMS <- PARAMS

# 4. Save data to disk
# ------------------------------------------------------------------------------
dir_stage2 <- "data/processed/Stage_2"
if(!dir.exists(dir_stage2)) dir.create(dir_stage2, recursive = TRUE)

# [DECISION 3] Keep FIXED name to avoid breaking automatic pipeline
output_file <- file.path(dir_stage2, "simulacion_wc.rds") 
saveRDS(sim_data, file = output_file)

# 5. Quality Report (Sanity Check)
# ------------------------------------------------------------------------------
# [IMPROVEMENT 2] Check achieved heritability
h2_achieved <- sim_data$meta$h2_real_empirico

# [FIX] Corrected ifelse syntax (added 'yes' and 'no' return values)
status_check <- ifelse(abs(h2_achieved - PARAMS$h2) < 0.05, "[OK]", "[WARN]")

message("--------------------------------------------------")
message(" DATA GENERATED AND SAVED ")
message("--------------------------------------------------")
message(paste(" File         :", output_file))
message(paste(" Patients     :", nrow(sim_data$X)))
message(paste(" Target h2    :", PARAMS$h2))
message(paste(" Achieved h2  :", round(h2_achieved, 4), status_check))
message(paste(" Config saved in: sim_data$PARAMS"))
message("--------------------------------------------------")