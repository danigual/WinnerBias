# analysis/Stage_2/1.4_check_bias.R

library(dplyr)
library(ggplot2)

# 1. Load results
gwas_results <- readRDS("data/processed/Stage_2/gwas_results.rds")

# 2. Data Preparation
GWAS_THRESHOLD <- 5e-8
winners <- subset(gwas_results, pval < GWAS_THRESHOLD & is_causal)

if(nrow(winners) == 0) {
  stop("STOPPED: No significant winners. Cannot calculate bias.")
}

# 3. Calculate Metrics
winners$z_score <- abs(winners$beta_hat / winners$se)

winners <- winners %>%
  mutate(
    bias = beta_hat - true_beta,
    bias_sq = (beta_hat - true_beta)^2,
    bias_pct = ((beta_hat - true_beta) / true_beta) * 100,
    type = ifelse(abs(beta_hat) > abs(true_beta), "Inflated (Curse)", "Underestimated")
  )

# 4. REPORTING
mse_global     <- mean(winners$bias_sq)
rmse_global    <- sqrt(mse_global)
mean_abs_inflation <- mean(abs(winners$bias_pct))

message("\n==================================================")
message("    EVALUATION METRICS (Method: A. Forde)         ")
message("==================================================")
message(paste(" Total Winners    :", nrow(winners)))
message(paste(" RMSE (Error)     :", round(rmse_global, 6)))
message(paste(" Mean Inflation   :", round(mean_abs_inflation, 2), "%"))
message("==================================================")

# 5. Save Files
report_dir <- "output/Stage_2"
if(!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

metrics_report <- data.frame(
  Metric = c("MSE", "RMSE", "Mean_Inflation_Pct", "N_Winners"), 
  Value = c(mse_global, rmse_global, mean_abs_inflation, nrow(winners))
)
write.csv(metrics_report, file.path(report_dir, "forde_metrics.csv"), row.names = FALSE)
write.csv(winners, file.path(report_dir, "winners_curse_report.csv"), row.names = FALSE)

# 6. Plot Z-score vs Bias 
# ------------------------------------------------------------------------------
p_forde <- ggplot(winners, aes(x = z_score, y = bias_pct)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Points
  geom_point(aes(color = type), size = 4, alpha = 0.8) +
  
  scale_color_manual(values = c("Inflated (Curse)" = "#E41A1C", "Underestimated" = "#377EB8")) +
  
  labs(
    title = "Winner's Curse Magnitude",
    subtitle = paste0("Observed inflation for the ", nrow(winners), " detected variants."),
    y = "Inflation (%)", 
    x = "Z-score (Signal Strength)",
    color = "Bias Type"
  ) +
  
  # --- CLEAN THEME ---
  theme_classic() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Save
fig_dir <- "output/figures/Stage_2"
if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

ggsave(file.path(fig_dir, "bias_vs_zscore.png"), p_forde, width = 7, height = 5, dpi = 300)
message(paste("Diagnostic plot saved:", file.path(fig_dir, "bias_vs_zscore.png")))