# analysis/Stage_2/1.2_run_gwas.R

# 1. Load tools
source("R/gwas_engine.R")

# 2. Load data
input_file <- "data/processed/Stage_2/simulacion_wc.rds"

# Safety check
if(!file.exists(input_file)) {
  stop("FATAL ERROR! 'simulacion_wc.rds' not found.\n    -> Please run first: analysis/Stage_2/1.1_generate_datasets.R")
}

message(">>> Loading simulation data...")
sim_data <- readRDS(input_file)

# 3. RUN GWAS (Timer included)
# ------------------------------------------------------------------------------
message(paste0(">>> Running GWAS on ", ncol(sim_data$X), " variants..."))

start_time <- Sys.time()

# Call the vectorized engine
gwas_results <- run_gwas(X = sim_data$X, y = sim_data$y)

end_time <- Sys.time()
duration <- round(difftime(end_time, start_time, units = "secs"), 2)

# 4. Add "Answer Key" (Ground Truth)
# ------------------------------------------------------------------------------
# This allows us to calculate bias in step 1.4
gwas_results$true_beta <- sim_data$true_betas
gwas_results$is_causal <- sim_data$true_betas != 0

# Retrieve original parameters if they exist (traceability)
if(!is.null(sim_data$PARAMS)) {
  attr(gwas_results, "simulation_params") <- sim_data$PARAMS
}

# 5. Save results
# ------------------------------------------------------------------------------
output_dir  <- "data/processed/Stage_2"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file <- file.path(output_dir, "gwas_results.rds")

saveRDS(gwas_results, file = output_file)

# 6. IMMEDIATE REPORT (Observability)
# ------------------------------------------------------------------------------
# Count how many passed the GWAS threshold (5e-8)
n_significant <- sum(gwas_results$pval < 5e-8)
n_causal_found <- sum(gwas_results$pval < 5e-8 & gwas_results$is_causal)

message("\n==================================================")
message("    GWAS COMPLETED SUCCESSFULLY ")
message("==================================================")
message(paste(" Execution time   :", duration, "seconds"))
message(paste(" Total Variants   :", nrow(gwas_results)))
message(paste(" Significant Hits (Winners) :", n_significant))
message(paste("    - Real (True Positives)    :", n_causal_found))
message(paste("    - Noise (False Positives)  :", n_significant - n_causal_found))
message("--------------------------------------------------")
message(paste(" Results saved in:", output_file))
message("==================================================\n")

# Early warning if no winners to correct
if(n_significant == 0) {
  warning("WARNING: No significant variants found. Winner's Curse analysis will have no data to correct.")
}