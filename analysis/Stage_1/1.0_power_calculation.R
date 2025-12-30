
# analysis/Stage_1/1.0_power_calculation.R

# load the required libraries and functions

source("R/statistical_power.R")

# 1. create parameters
PARAMS <- list(
  n_samples = 5000,
  n_snps    = 10000,
  n_causal  = 50,
  h2        = 0.5,
  seed      = 42
  maf <- 0.3
)

# calculate beta from h2 and n_causal
beta <- sqrt(PARAMS$h2 / (PARAMS$n_causal * 2 * maf * (1 - maf)))

# Umbral GWAS
alpha <- 0.0001

# Poder estadístico
power <- get_gwas_power(
  n     = PARAMS$n_samples,
  beta  = beta,
  maf   = maf,
  alpha = alpha
)

power