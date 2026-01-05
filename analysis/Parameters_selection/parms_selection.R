# analysis/parms_selection.R

#' A graph is created showing how statistical power, of each SNP in a GWAS, changes with 
#' sample size and effect size.
#' Taking into account a MAF = 0.3 and an alpha = 5e-8.
#' 
#' These graph will help us determine the parameters to use in each scenario.

# load the required libraries and functions

library(ggplot2)
library(dplyr)

source("R/statistical_power.R")

# parameters
Ns_to_test     <- c(2000, 7000, 50000)
betas_to_test  <- seq(0, 1, by = 0.1)
maf_fixed      <- 0.3

# create table of combinations
grid <- expand.grid(
  N     = Ns_to_test,
  beta  = betas_to_test
)

# set alpha
grid$alpha <- 5e-8 

# calculate statistical power
grid$power <- mapply(
  get_gwas_power,
  n     = grid$N,
  beta  = grid$beta,
  maf   = maf_fixed,
  alpha = grid$alpha
)

# plot
p_power <- ggplot(grid, aes(x = beta, y = power)) +
  geom_line(size = 1.2, color="blue") +
  facet_wrap(~ N, nrow = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "Potencia estadística de cada SNP en función de N y β ",
    subtitle = "MAF = 0.3, alpha = 5e-8 ",
    x = "Beta (tamaño del efecto)",
    y = "Potencia")

print(p_power)

# save DATA

if(!dir.exists("output")) dir.create("output")
ggsave("output/figures/power_curves.png", p_power, width = 8, height = 6)

print(p_power)
message("Gráfico de potencia en función de N, beta y N_SNPs")

