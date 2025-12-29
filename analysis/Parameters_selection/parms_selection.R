# analysis/parms_selection.R

#' A graph is created showing how statistical power changes with 
#' sample size, number of SNPs tested and effect size.
#' Taking into account an MAF = 0.3 and an alpha = 0.05.
#' 
#' These graph will help us determine the parameters to use in each scenario.

# load the required libraries and functions

library(ggplot2)
library(dplyr)

source("R/statistical_power.R")

# parameters
Ns_to_test     <- c(500, 1000, 5000, 10000)
betas_to_test  <- seq(0, 1, by = 0.1)
Nsnps_to_test  <- c(1e2, 1e3, 1e4)
maf_fixed      <- 0.3

# create table of combinations
grid <- expand.grid(
  N     = Ns_to_test,
  beta  = betas_to_test,
  Nsnps = Nsnps_to_test
)

# calculate Bonferroni alpha
grid$alpha <- 0.05 / grid$Nsnps

# calculate statistical power
grid$power <- mapply(
  get_gwas_power,
  n     = grid$N,
  beta  = grid$beta,
  maf   = maf_fixed,
  alpha = grid$alpha
)

# plot
p_power <- ggplot(grid, aes(x = beta, y = power, color = as.factor(Nsnps))) +
  geom_line(size = 1.2) +
  facet_wrap(~ N, nrow = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "Potencia estadística según N, β y número total de SNPs",
    subtitle = "MAF fijo = 0.3, umbral Bonferroni = 0.05 / N_SNPs",
    x = "Beta (tamaño del efecto)",
    y = "Potencia",
    color = "N SNPs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  )

print(p_power)

# Save

if(!dir.exists("output")) dir.create("output")
ggsave("output/figures/power_curves.png", p_power, width = 8, height = 6)

print(p_power)
message("Gráfico de potencia en función de N, beta y N_SNPs")

