# analysis/Stage_2/1.5_correct_bias.R

# ==============================================================================
# SCRIPT 1.5: WINNER'S CURSE CORRECTION (BOOTSTRAP)
# ==============================================================================
# Objective: Correct the upward bias in effect size estimates for significant variants
# Method: Parametric Bootstrap (Resampling with replacement)
# Dependencies: Requires 'simulacion_wc.rds' and 'gwas_results.rds'
# ==============================================================================

# 1. Setup & Data Loading
# ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

# Load the modular bootstrap function we created in R/
source("R/bootstrap_correction.R")

message("\n[Step 1] Loading Stage 2 simulation data...")
sim_data_path <- "data/processed/Stage_2/simulacion_wc.rds"
gwas_res_path <- "data/processed/Stage_2/gwas_results.rds"

if(!file.exists(sim_data_path) | !file.exists(gwas_res_path)) {
  stop("CRITICAL ERROR: Input files not found. Please run scripts 1.1 and 1.2 first.")
}

sim_data <- readRDS(sim_data_path)
gwas_res <- readRDS(gwas_res_path)

# 2. Pre-processing: Identify the 'Winners'
# ------------------------------------------------------------------------------
# We apply FDR correction to determine which SNPs are "significant"
gwas_res$qval <- p.adjust(gwas_res$pval, method = "fdr")

# Filter: We only correct variants that passed the significance threshold
winners <- gwas_res %>% 
  filter(qval < 0.05)

n_winners <- nrow(winners)
if(n_winners == 0) stop("No significant variants found. Nothing to correct.")

message(paste("[Step 2] Found", n_winners, "significant variants (Winners)."))

# 3. EXECUTION: Run Bootstrap Correction
# ------------------------------------------------------------------------------
# We use our modular function here. 
# n_boot=200.
message("[Step 3] Running Bootstrap Correction (this may take a moment)...")

winners_corrected <- run_bootstrap_correction(
  winners_df  = winners,
  X_genotype  = sim_data$X,
  y_phenotype = sim_data$y,
  n_boot      = 200,  
  seed        = 123
)

# 4. Evaluation: Calculate Improvement Metrics
# ------------------------------------------------------------------------------
# We compare the error of the Original (Naive) vs. Corrected estimates against the Truth.
winners_corrected <- winners_corrected %>%
  mutate(
    error_naive     = abs(beta_hat - true_beta),
    error_corrected = abs(beta_corrected - true_beta),
    is_improved     = error_corrected < error_naive
  )

# Calculate RMSE (Root Mean Square Error) - The gold standard metric for bias
rmse_naive     <- sqrt(mean(winners_corrected$error_naive^2))
rmse_corrected <- sqrt(mean(winners_corrected$error_corrected^2))
pct_improvement <- mean(winners_corrected$is_improved) * 100

message("\n==================================================")
message("             CORRECTION RESULTS                   ")
message("==================================================")
message(paste(" RMSE (Original) :", round(rmse_naive, 4)))
message(paste(" RMSE (Corrected):", round(rmse_corrected, 4)))
message(paste(" Improvement     :", round(rmse_naive - rmse_corrected, 4), "points"))
message(paste(" Success Rate    :", round(pct_improvement, 2), "% of variants improved"))
message("==================================================\n")

# 5. Visualization: The 'Before & After' Plot
# ------------------------------------------------------------------------------
p_correction <- ggplot(winners_corrected) +
  # Identity line (Perfect estimation)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  
  # Arrows showing the correction path
  geom_segment(aes(x = true_beta, xend = true_beta, 
                   y = beta_hat, yend = beta_corrected),
               arrow = arrow(length = unit(0.2, "cm")), color = "grey60", alpha = 0.6) +
  
  # Points: Red for biased, Green for corrected
  geom_point(aes(x = true_beta, y = beta_hat, color = "Original Estimate"), 
             size = 2.5, alpha = 0.8) +
  geom_point(aes(x = true_beta, y = beta_corrected, color = "Bootstrap Corrected"), 
             size = 2.5, alpha = 0.8) +
  
  scale_color_manual(values = c("Original Estimate" = "#E41A1C", # Red
                                "Bootstrap Corrected" = "#4DAF4A")) + # Green
  
  labs(
    title = "Winner's Curse Correction: Stage 2",
    subtitle = paste0("N = ", nrow(sim_data$X), " | RMSE reduced from ", 
                      round(rmse_naive, 3), " to ", round(rmse_corrected, 3)),
    x = "True Effect Size (Truth)",
    y = "Estimated Effect Size (Observed)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# 6. Save Outputs
# ------------------------------------------------------------------------------
# Save plot
fig_dir <- "output/figures/Stage_2"
if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
ggsave(file.path(fig_dir, "bias_correction_stage2.png"), p_correction, width = 7, height = 6)

# Save data for further analysis
saveRDS(winners_corrected, "data/processed/Stage_2/winners_corrected.rds")

message(paste("[Step 6] Analysis complete. Plot saved to", file.path(fig_dir, "bias_correction_stage2.png")))