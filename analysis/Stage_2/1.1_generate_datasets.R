# analysis/Stage_2/1.1_generate_datasets.R

# 1. Dynamic Configuration & Tools
# ------------------------------------------------------------------------------
# If run standalone (without Master Script), default to Stage 2
if(!exists("PARAMS")) {
  STAGE_NAME <- "Stage_2"
  PARAMS <- list(
    n_samples = 1000, 
    n_snps    = 10000, 
    n_causal  = 100, 
    h2        = 0.5, 
    mean_beta = 0.125, 
    seed      = 42
  )
  message("⚠ Running in Standalone Mode (Default: Stage_2)")
}

source("R/sim_genetics.R")

# Define Dynamic Paths (This adapts to Stage_0, Stage_1, etc.)
output_dir <- file.path("data/processed", STAGE_NAME)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "simulacion_wc.rds")

# 2. Run simulation
# ------------------------------------------------------------------------------
# We use 'PARAMS' directly (it comes either from Standalone block or Master Script)
message(paste0("Generating simulated data for ", STAGE_NAME, " (N=", PARAMS$n_samples, ", Beta=", PARAMS$mean_beta, ")..."))

sim_data <- generate_gwas_data(
  n_samples = PARAMS$n_samples,
  n_snps    = PARAMS$n_snps,
  n_causal  = PARAMS$n_causal,
  h2        = PARAMS$h2,
  mean_beta = PARAMS$mean_beta,
  seed      = PARAMS$seed
)

# Inject parameters for traceability
sim_data$PARAMS <- PARAMS

# 3. Save data to disk
# ------------------------------------------------------------------------------
# We use the DYNAMIC 'output_file' defined at the top
saveRDS(sim_data, file = output_file)

# 4. Quality Report (Sanity Check)
# ------------------------------------------------------------------------------
h2_achieved <- sim_data$meta$h2_real_empirico
status_check <- ifelse(abs(h2_achieved - PARAMS$h2) < 0.05, "[OK]", "[WARN]")

message("--------------------------------------------------")
message(" DATA GENERATED AND SAVED ")
message("--------------------------------------------------")
message(paste(" Stage        :", STAGE_NAME))
message(paste(" File         :", output_file))
message(paste(" Patients     :", nrow(sim_data$X)))
message(paste(" Target h2    :", PARAMS$h2))
message(paste(" Achieved h2  :", round(h2_achieved, 4), status_check))
message("--------------------------------------------------")