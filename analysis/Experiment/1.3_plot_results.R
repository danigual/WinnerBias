# analysis/Experiment/1.3_plot_results.R

# 1. Load libraries and Dynamic Config
library(ggplot2)
library(dplyr)

# Dynamic Configuration
# If run standalone, defaults to Experiment. If run by Master Script, uses the active Stage.
if(!exists("STAGE_NAME")) STAGE_NAME <- "Experiment"

# 2. Load Results
# Dynamic path: reads from the folder corresponding to the active Stage
results_file <- file.path("data/processed", STAGE_NAME, "gwas_results.rds")

if(!file.exists(results_file)) {
  stop(paste0("ERROR! Results missing (gwas_results.rds) in ", STAGE_NAME, 
              ".\n    -> Please run 1.2 first."))
}

res <- readRDS(results_file)

# 3. Define "Significant" (Consistent with Step 1.2)
# ------------------------------------------------------------------------------
# Using standard GWAS threshold (Bonferroni)
GWAS_THRESHOLD <- 5e-8

res$is_significant <- res$pval < GWAS_THRESHOLD

# 4. Categorize points (Labeling)
# ------------------------------------------------------------------------------
res <- res %>%
  mutate(category = case_when(
    is_significant & is_causal  ~ "Winner (Detected Causal)",
    is_significant & !is_causal ~ "False Positive",
    !is_significant & is_causal ~ "Missed Causal (Low Power)",
    TRUE ~ "Background Noise"
  ))

# Convert to factor to control legend order
res$category <- factor(res$category, levels = c(
  "Background Noise", 
  "Missed Causal (Low Power)", 
  "False Positive", 
  "Winner (Detected Causal)"
))

# Statistics for subtitle
n_winners <- sum(res$category == "Winner (Detected Causal)")
n_fp      <- sum(res$category == "False Positive")

# 5. PLOT: The Anatomy of Winner's Curse (CLEAN STYLE)
# ------------------------------------------------------------------------------
p1 <- ggplot(res, aes(x = true_beta, y = beta_hat)) +
  
  # A. Noise Layer (Grey, very transparent)
  geom_point(data = subset(res, category == "Background Noise"), 
             color = "grey90", alpha = 0.3, size = 1) + 
  
  # B. Missed Causal Layer (Blue)
  geom_point(data = subset(res, category == "Missed Causal (Low Power)"), 
             color = "#377EB8", alpha = 0.6, size = 2, shape = 1) +
  
  # C. False Positive Layer (Orange)
  geom_point(data = subset(res, category == "False Positive"), 
             color = "orange", alpha = 0.8, size = 3, shape = 17) + 
  
  # D. WINNERS Layer (Red)
  geom_point(data = subset(res, category == "Winner (Detected Causal)"), 
             color = "#E41A1C", size = 3, alpha = 0.9) +
  
  # E. Visual References
  # Identity Line (The perfect truth)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = 0.45, y = 0.4, label = "Ideal Estimation", angle = 45, hjust=0, size=3, color="grey30") +
  
  # Trend Line for Winners (Shows Bias)
  geom_smooth(data = subset(res, category == "Winner (Detected Causal)"), 
              method = "lm", se = FALSE, color = "#E41A1C", linetype = "dotted", fullrange=TRUE, linewidth=0.8) +
  
  # Labels
  labs(
    title = paste0("Winner's Curse Effect: Truth vs. Estimation (", STAGE_NAME, ")"),
    subtitle = paste0("Threshold P < 5e-8 | Winners: ", n_winners, " | False Positives: ", n_fp, 
                      "\nRed points (Winners) tend to be inflated away from the diagonal."),
    x = "TRUE Effect (True Beta)",
    y = "ESTIMATED Effect (GWAS Beta)",
    color = "Category"
  ) +
  
  # --- CLEAN THEME APPLIED ---
  theme_classic() + 
  #coord_fixed(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6)) + 
  coord_fixed() + 
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold"),
    legend.background = element_rect(fill = "white", color = NA)
  )

# 6. Save and Report
# ------------------------------------------------------------------------------
# Dynamic output directory
output_dir <- file.path("output/figures", STAGE_NAME)

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "winners_curse_plot.png"), plot = p1, width = 8, height = 7, dpi = 300)

# Summary table to console
message("\n--- CLASSIFICATION SUMMARY ---")
print(table(res$category))
message("--------------------------------")
message(paste("Plot saved to:", file.path(output_dir, "winners_curse_plot.png")))