# analysis/run_experiment_batch.R
# ==============================================================================
# MASTER SCRIPT: AUTOMATED PIPELINE RUNNER (Stages 0, 1, 2)
# ==============================================================================
# This script orchestrates the full execution of all thesis scenarios.
# It defines parameters for each Stage and injects global variables 
# so the 'generic' scripts (in analysis/Stage_2/) adapt dynamically.

# 1. Define Scenarios (Thesis Parameters)
# ------------------------------------------------------------------------------
SCENARIOS <- list(
  
  # Stage 0
  "Stage_0" = list(
    n_samples = 50000, 
    n_snps    = 1000,
    n_causal  = 100,
    h2        = 0.5,
    mean_beta = 0.5,    
    seed      = 111
  ),
  
  # Stage 1
  "Stage_1" = list(
    n_samples = 7000,
    n_snps    = 1000,
    n_causal  = 100,
    h2        = 0.5,
    mean_beta = 0.5,     
    seed      = 222
  ),
  
  # Stage 2
  "Stage_2" = list(
    n_samples = 2000,
    n_snps    = 1000,
    n_causal  = 100,
    h2        = 0.5,
    mean_beta = 0.5,   
    seed      = 42
  )
)

# 2. Scripts to Execute (The Pipeline)
# ------------------------------------------------------------------------------
# We use the scripts in analysis/Experiment/ as "Master Templates".
# Thanks to parameterization, these scripts will adapt to the active Stage.
scripts_to_run <- c(
  "analysis/Experiment/1.0_power_calculation.R",  # Generate the estimated power
  "analysis/Experiment/1.1_generate_datasets.R",  # Generate simulation
  "analysis/Experiment/1.2_run_gwas.R",           # Run GWAS
  "analysis/Experiment/1.3_plot_results.R",       # Plot Real vs Estimated
  "analysis/Experiment/1.4_check_bias.R",         # Bias Diagnosis
  "analysis/Experiment/1.5_correct_bias.R"        # Bootstrap Correction Attempt
)

# 3. Execution Loop
# ------------------------------------------------------------------------------
message("STARTING BATCH EXECUTION OF THESIS SCENARIOS...")

for(stage_name in names(SCENARIOS)) {
  
  message(paste0("\n\n############################################################"))
  message(paste0("   EXECUTING SCENARIO: ", stage_name))
  message(paste0("############################################################\n"))
  
  # A. Inject Global Variables
  # ----------------------------------------------------------------------------
  # These variables (<<-) will be visible inside the sourced scripts.
  STAGE_NAME <<- stage_name          
  PARAMS     <<- SCENARIOS[[stage_name]] 
  
  # B. Execute Script Sequence
  # ----------------------------------------------------------------------------
  for(script in scripts_to_run) {
    # Safety check
    if(!file.exists(script)) stop(paste("CRITICAL ERROR: Script not found:", script))
    
    message(paste(">>> Running module:", basename(script)))
    
    # 'source' runs the file code as if it were written here.
    # local=FALSE is KEY to sharing the environment and seeing STAGE_NAME.
    source(script, local = FALSE) 
  }
  
  message(paste0(stage_name, " completed successfully."))
}

message("\n============================================================")
message("   ALL DONE! ALL SCENARIOS HAVE BEEN PROCESSED.")
message("============================================================")
message("   Results saved in: data/processed/ and output/figures/")