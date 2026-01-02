# analysis/Stage_2/1.1_power_calculation.R

#' Before conducting the experiment, the statistical power is calculated in 
#' order to identify its quality. if it is 0.8, 
#' the experiment has good statistical power.

# load required libraries and functions

source("R/statistical_power.R")

# 1. Dynamic Configuration & Tools
# ------------------------------------------------------------------------------
# If run standalone (without Master Script), define defaults (Stage 2)

if (!exists("PARAMS")) {
  STAGE_NAME <- "Stage_2"
  
  PARAMS <- list(
    n_samples = 2000, 
    n_snps    = 1000, 
    n_causal  = 100, 
    h2        = 0.5, 
    mean_beta = 0.5, 
    seed      = 42,
    maf       = 0.3
  )
  message("⚠ Running in Standalone Mode (Default: Stage_2)")
}

# 2. Set parameters
#-------------------------------------------------------------------------------
# set seed, alpha and MAF
set.seed(PARAMS$seed)
alpha_gwas <- 5e-8 / PARAMS$n_snps

# 3. simulate data
#-------------------------------------------------------------------------------
# simulate heterogeneous betas
beta_raw <- rnorm(PARAMS$n_causal, mean = 0, sd = PARAMS$mean_beta)

# genetic variance induced by these betas
genetic_var_raw <- sum(
  2 * 0.3 * (1 - 0.3) * beta_raw^2
)

# scale factor so that the total h2 is as desired
scale_factor <- sqrt(PARAMS$h2 / max(genetic_var_raw, 1e-12))
beta_causal <- beta_raw * scale_factor

# 4. Calculate statistical power
#-------------------------------------------------------------------------------
# power by causal SNP
power_per_snp <- mapply(
  get_gwas_power,
  n     = PARAMS$n_samples,
  beta  = abs(beta_causal),
  maf   = 0.3,
  alpha = alpha_gwas
)

# 5. Output
#-------------------------------------------------------------------------------

mean_power    <- mean(power_per_snp)
expected_hits <- sum(power_per_snp)

output <- list(
  beta_causal    = beta_causal,
  power_per_snp  = power_per_snp,
  mean_power     = mean_power,
  expected_hits  = expected_hits
)

output

# global output

output2 <- c(
  mean_power   = output$mean_power,
  expected_hits = output$expected_hits
)

global_df <- as.data.frame(t(output2))

# 6. Save data
#-------------------------------------------------------------------------------

# Save CSV in output/
output_dir <- file.path("output", STAGE_NAME)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.csv(global_df, file.path(output_dir, "statistical_power.csv"), row.names = FALSE)

# Save full RDS in data/processed/
data_dir <- file.path("data", "processed", STAGE_NAME)
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
saveRDS(output, file.path(data_dir, "power_details.rds"))


