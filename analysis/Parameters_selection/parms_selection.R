# analysis/parms_selection.R

#' A graph is created showing how statistical power changes with 
#' sample size, minor allele frequency, and effect size.
#' 
#' These graph help us determine the parameters to use in each scenario.

# load the required libraries and functions

library(ggplot2)
library(dplyr)

source("R/statistical_power.R")

# create table and calculate power

Ns_to_test <- c(1000, 5000, 10000, 20000)
betas_to_test <- seq(0, 1, by = 0.1)
mafs_to_test <- c(0.1, 0.3, 0.5)
alpha <- 1e-4

grid <- expand.grid(
  N = Ns_to_test,
  beta = betas_to_test,
  maf = mafs_to_test
)

grid$power <- mapply(
  get_gwas_power,
  n = grid$N,
  beta = grid$beta,
  maf = grid$maf,
  alpha = alpha
)

# plot

p_power<- ggplot(grid, aes(x = beta, y = power, color = as.factor(N))) +
  geom_line(size = 1.1) +
  facet_wrap(~ maf, nrow = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "Potencia en función de N, β y MAF",
    x = "Beta (magnitud del efecto)",
    y = "Potencia Estadística",
    color = "Tamaño muestral"
  ) +
  theme_minimal()

# Save

if(!dir.exists("output")) dir.create("output")
ggsave("output/figures/power_curves.png", p_power, width = 8, height = 6)

print(p_power)
message("Gráfico de potencia en función de N, beta y MAF")

