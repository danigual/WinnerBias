# analysis/Experiment/1.5_correct_bias.R

library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Dynamic Configuration & Load Tools
# ------------------------------------------------------------------------------
if(!exists("STAGE_NAME")) STAGE_NAME <- "Experiment"

if(!file.exists("R/bootstrap_correction.R")) stop("ATTENTION! R/bootstrap_correction.R is missing")
source("R/bootstrap_correction.R")

# Dynamic Paths
sim_file     <- file.path("data/processed", STAGE_NAME, "simulacion_wc.rds")
gwas_file    <- file.path("data/processed", STAGE_NAME, "gwas_results.rds")
out_fig_dir  <- file.path("output/figures", STAGE_NAME)
out_data_dir <- file.path("data/processed", STAGE_NAME)

# Create output directories if they don't exist
if(!dir.exists(out_fig_dir)) dir.create(out_fig_dir, recursive = TRUE)
if(!dir.exists(out_data_dir)) dir.create(out_data_dir, recursive = TRUE)

# 2. Load Data
# ------------------------------------------------------------------------------
message(paste0(">>> Loading simulation data and GWAS results (", STAGE_NAME, ")..."))

if(!file.exists(sim_file) || !file.exists(gwas_file)) {
  stop(paste0("STOPPED: Data missing in ", STAGE_NAME, ". Run 1.1 and 1.2 first."))
}

sim_data <- readRDS(sim_file)
gwas_res <- readRDS(gwas_file)

# 3. Filter Winners
# ------------------------------------------------------------------------------
GWAS_THRESHOLD <- 5e-8
winners <- gwas_res %>% filter(pval < GWAS_THRESHOLD)

if(nrow(winners) == 0) {
  stop(paste0("No significant winners in ", STAGE_NAME, ". Check N or Betas in 1.1"))
}

message(paste(">>> Winners detected for correction:", nrow(winners)))

# 4. BOOTSTRAP EXECUTION (500 iterations)
# ------------------------------------------------------------------------------
message(">>> Running Bootstrap (500 iterations)...")
winners_corrected <- run_bootstrap_correction(
  winners_df  = winners,
  X_genotype  = sim_data$X,
  y_phenotype = sim_data$y,
  n_boot      = 500, 
  seed        = 123
)

# 5. Improvement Evaluation
# ------------------------------------------------------------------------------
winners_corrected <- winners_corrected %>%
  mutate(
    error_naive     = abs(beta_hat - true_beta),
    error_corrected = abs(beta_corrected - true_beta),
    is_improved     = error_corrected < error_naive
  )

rmse_naive     <- sqrt(mean(winners_corrected$error_naive^2))
rmse_corrected <- sqrt(mean(winners_corrected$error_corrected^2))
improvement    <- rmse_naive - rmse_corrected
pct_improved   <- mean(winners_corrected$is_improved) * 100

# Confidence Intervals (95%)
winners_corrected$ci_lower <- winners_corrected$beta_corrected - 1.96 * winners_corrected$beta_corrected_se
winners_corrected$ci_upper <- winners_corrected$beta_corrected + 1.96 * winners_corrected$beta_corrected_se

message("\n==================================================")
message(paste0("             CORRECTION RESULTS (", STAGE_NAME, ")"))
message("==================================================")
message(paste(" Original RMSE   :", round(rmse_naive, 4)))
message(paste(" Corrected RMSE  :", round(rmse_corrected, 4)))
message(paste(" Net Improvement :", round(improvement, 4)))
message("==================================================")

# 6. Visualization (CLEAN STYLE)
# ------------------------------------------------------------------------------

# Config
x_range <- diff(range(winners_corrected$true_beta))
bar_width_dynamic <- if(x_range == 0) 0.01 else x_range * 0.015
TEXT_SIZE <- 3.5 

# Arrow Logic
MAX_ARROWS <- 50 
if(nrow(winners_corrected) > MAX_ARROWS) {
  arrow_data <- winners_corrected %>%
    mutate(magnitude = abs(beta_hat - beta_corrected)) %>%
    slice_max(order_by = magnitude, n = MAX_ARROWS)
  subtitle_text <- paste0("Arrows shown for Top ", MAX_ARROWS, " largest changes.")
} else {
  arrow_data <- winners_corrected
  subtitle_text <- "Arrows indicate bias reduction towards the dashed line."
}

# Stats Label
stats_label <- paste0(
  "Correction Results:\n",
  "--------------------\n",
  "RMSE Before : ", round(rmse_naive, 4), "\n",
  "RMSE After  : ", round(rmse_corrected, 4), "\n",
  "Improvement : ", round(improvement, 4), "\n",
  "Success Rate: ", round(pct_improved, 1), "%"
)

# Plot
p_corr <- ggplot(winners_corrected) +
  # Identity Line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.6, color = "grey40") +
  
  # Arrows
  geom_segment(data = arrow_data,
               aes(x = true_beta, xend = true_beta, y = beta_hat, yend = beta_corrected),
               arrow = arrow(length = unit(0.2, "cm")), color = "grey60", alpha = 0.5) +
  
  # Points: Original
  geom_point(aes(x = true_beta, y = beta_hat, color = "Original (Inflated)"), 
             size = 3, alpha = 0.7) +
  
  # Points: Corrected + Error Bars
  geom_errorbar(aes(x = true_beta, ymin = ci_lower, ymax = ci_upper), 
                color = "#4DAF4A", width = bar_width_dynamic, alpha = 0.6) +
  
  geom_point(aes(x = true_beta, y = beta_corrected, color = "Corrected (Bootstrap)"), 
             size = 3, alpha = 0.9) +
  
  # Text Box (Adjusted for white background)
  annotate("label", 
           x = -Inf, y = Inf, 
           label = stats_label, 
           hjust = -0.1, vjust = 1.1,  
           size = TEXT_SIZE, 
           family = "sans", alpha = 1, fill = "white", color = "black", label.size = 0.5) +
  
  scale_color_manual(name = "Estimate Type", 
                     values = c("Original (Inflated)" = "#E41A1C", "Corrected (Bootstrap)" = "#4DAF4A")) +
  
  labs(title = paste0("Winner's Curse Correction - ", STAGE_NAME),
       subtitle = subtitle_text,
       x = "TRUE Effect (Simulated)", 
       y = "ESTIMATED Effect (GWAS)") +
  
  # --- CLEAN THEME ---
  theme_classic() + 
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Save
ggsave(file.path(out_fig_dir, "bias_correction.png"), p_corr, width = 8, height = 7, dpi = 300)
saveRDS(winners_corrected, file.path(out_data_dir, "winners_corrected.rds"))

message(paste("Plot saved to", file.path(out_fig_dir, "bias_correction.png")))
print(p_corr)