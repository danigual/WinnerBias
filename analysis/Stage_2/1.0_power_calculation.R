# analysis/Stage_2/1.0_power_calculation.R

#' Calculate the statistical power of Stage_2 before running the experiment.
#' Justifies why we expect Winner's Curse (Low Power).

library(ggplot2)
library(dplyr)

source("R/statistical_power.R")

# 1. Scenario Configuration (Stage 2: Low Power Scenario)
n_samples_stage2 <- 1000  # Low N to induce Winner's Curse
n_causal_total   <- 100
s_desired        <- 20    # Number of causal SNPs we ideally want to detect

PARAMS <- list(
  n_samples = n_samples_stage2,
  n_snps    = 10000,
  n_causal  = n_causal_total,
  h2        = 0.5,
  alpha     = 5e-8,
  target_power = s_desired / n_causal_total 
)

# 2. Inverse Calculation (Beta required for 80% Power)
# We specifically use MAF = 0.3 because that represents our simulated causal variants
maf_simulation <- 0.3

get_beta_for_power <- function(target_power, n, maf, alpha) {
  f <- function(beta) {
    get_gwas_power(n = n, beta = beta, maf = maf, alpha = alpha) - target_power
  }
  tryCatch(uniroot(f, interval = c(0, 5))$root, error = function(e) NA)
}

beta_80pct <- get_beta_for_power(0.80, PARAMS$n_samples, maf = maf_simulation, alpha = PARAMS$alpha)

# [CORRECTION] Added explicit mention of MAF in the console output
message(paste0("For 80% power (N=", PARAMS$n_samples, ", MAF=", maf_simulation, "), required Beta: ", round(beta_80pct, 3)))

# 3. Build Curves Table
betas_to_test <- seq(0, 0.6, by = 0.01)
mafs_to_test  <- c(0.05, 0.1, 0.3, 0.5)

power_data <- expand.grid(
  beta = betas_to_test,
  maf  = mafs_to_test
)

power_data$power <- mapply(
  get_gwas_power, 
  n     = PARAMS$n_samples, 
  beta  = power_data$beta, 
  maf   = power_data$maf, 
  alpha = PARAMS$alpha
)

# Factor for ordered legend
power_data$maf_label <- factor(
  paste("MAF =", power_data$maf),
  levels = paste("MAF =", sort(unique(power_data$maf)))
)

# 4. Plotting

beta_real_simulation <- 0.125

# Calculate exact power for your specific simulation point (to plot the dot)
power_simulation <- get_gwas_power(PARAMS$n_samples, beta_real_simulation, maf_simulation, PARAMS$alpha)

p_power <- ggplot(power_data, aes(x = beta, y = power, color = maf_label)) +
  # 1. Risk Zone (Background)
  annotate("rect", xmin = 0, xmax = 0.25, ymin = 0, ymax = 1, alpha = 0.1, fill = "red") +
  annotate("text", x = 0.02, y = 0.95, label = "WINNER'S CURSE\nZONE", 
           color = "#D73027", hjust=0, fontface = "bold", size = 3.5) +
  
  # 2. Reference Lines
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgray", linewidth = 0.6) + 
  annotate("text", x = 0.45, y = 0.82, label = "Ideal Power (80%)", color = "darkgray", size = 3.5) +
  
  # 3. Power Curves
  geom_line(linewidth = 1.2) +
  
  # [MEJORA CLAVE] Highlight specific Simulation Point
  # Vertical line
  geom_vline(xintercept = beta_real_simulation, linetype = "dotted", color = "blue", linewidth = 1) +
  
  # The specific DOT where your simulation lives (Beta 0.125, MAF 0.3)
  annotate("point", x = beta_real_simulation, y = power_simulation, color = "blue", size = 4) +
  
  # Label explaining the dot
  annotate("text", x = beta_real_simulation + 0.02, y = 0.45, 
           label = paste0("\n(Beta=", beta_real_simulation, ", MAF=", maf_simulation, ")\nPower ~ ", round(power_simulation*100,1), "%"), 
           color = "blue", hjust = 0, fontface = "italic", size = 3.5) +
  
  # 4. Aesthetics
  labs(
    title = paste0("Stage 2 Power Analysis (N = ", PARAMS$n_samples, ")"),
    subtitle = paste0("At Beta=", beta_real_simulation, " and MAF=", maf_simulation, ", statistical power is very low."),
    x = "Effect Size (True Beta)",
    y = "Statistical Power",
    color = "Frequency (MAF)"
  ) +
  
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 0.6), expand = c(0, 0)) +
  
  theme_classic() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.background = element_rect(fill = "white", color = NA)
  )

# 5. Save
output_dir <- "output/figures/Stage_2"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "power_curves.png"), p_power, width = 8, height = 6)

print(p_power)