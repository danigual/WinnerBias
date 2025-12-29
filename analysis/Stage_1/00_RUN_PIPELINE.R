# analysis/Stage_1/00_RUN_PIPELINE.R

# ==============================================================================
# MASTER PIPELINE: STAGE 1 (Winner's Curse Simulation & Correction)
# ==============================================================================
# DESCRIPTION:
# This script orchestrates the sequential execution of the Stage 1 analysis.
# It proceeds from synthetic data generation through GWAS to bias correction.
#
# EXECUTION FLOW:
# 1.1 -> Generate Data (simulacion_wc.rds)
# 1.2 -> Run GWAS (gwas_results.rds)
# 1.3 -> Diagnostic Plots (Figures)
# 1.4 -> Bias Quantification (Report CSV)
# 1.5 -> Bootstrap Correction (Final Algorithm)
# ==============================================================================

# 1. ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
# Clear memory to ensure a clean execution environment
rm(list = ls()) 
# Close any open graphics devices
graphics.off()

# Define current stage for logging
CURRENT_STAGE <- "Stage_1"

# Start global timer
start_global <- Sys.time()

message(paste0("\n============================================================"))
message(paste0("   STARTING PIPELINE EXECUTION - ", CURRENT_STAGE))
message(paste0("============================================================"))

# 2. SEQUENTIAL EXECUTION
# ------------------------------------------------------------------------------

# --- STEP 1.1: Data Generation ---
message("\n>>> [1/5] Running 1.1: Dataset Generation...")
source("analysis/Stage_1/1.1_generate_datasets.R")

# --- STEP 1.2: GWAS ---
message("\n>>> [2/5] Running 1.2: GWAS Engine...")
source("analysis/Stage_1/1.2_run_gwas.R")

# --- STEP 1.3: Visualization ---
message("\n>>> [3/5] Running 1.3: Result Plots...")
source("analysis/Stage_1/1.3_plot_results.R")

# --- STEP 1.4: Bias Check ---
message("\n>>> [4/5] Running 1.4: Bias Check and Inflation Report...")
source("analysis/Stage_1/1.4_check_bias.R")

# --- STEP 1.5: Bootstrap Correction ---
message("\n>>> [5/5] Running 1.5: Bootstrap Correction (Winner's Curse)...")
source("analysis/Stage_1/1.5_correct_bias.R")


# 3. FINAL SUMMARY
# ------------------------------------------------------------------------------
end_global <- Sys.time()
total_time <- round(difftime(end_global, start_global, units = "mins"), 2)

message("\n============================================================")
message("   PIPELINE COMPLETED SUCCESSFULLY")
message("============================================================")
message(paste("   Total execution time:", total_time, "minutes."))
message(paste("   Outputs saved in:"))
message(paste("   - Data:    data/processed/", CURRENT_STAGE, "/", sep=""))
message(paste("   - Figures: output/figures/", CURRENT_STAGE, "/", sep=""))
message(paste("   - Reports: output/", CURRENT_STAGE, "/", sep=""))
message("============================================================\n")