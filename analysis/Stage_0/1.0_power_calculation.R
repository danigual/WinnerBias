# analysis/Stage_0/1.1_generate_datasets.R

#' Before conducting the experiment, the statistical power is calculated in 
#' order to identify its quality. if it is 0.8, 
#' the experiment has good statistical power.

# load required libraries and functions

source("R/statistical_power.R")

# set the parameters for Stage_1

PARAMS <- list(
  n_samples = 5000,
  n_snps    = 10000,
  n_causal  = 50,
  h2        = 0.5,
  seed      = 42
)

# set seed, alpha
set.seed(PARAMS$seed)
alpha_gwas <- 0.05 / PARAMS$n_snps
maf_causal <- 0.3

# 1. simulate heterogeneous betas
beta_raw <- rnorm(PARAMS$n_causal, mean = 0, sd = 1)

# genetic variance induced by these betas
genetic_var_raw <- sum(
  2 * maf_causal * (1 - maf_causal) * beta_raw^2
)

# scale factor so that the total h2 is as desired
scale_factor <- sqrt(PARAMS$h2 / genetic_var_raw)

beta_causal <- beta_raw * scale_factor

# 2. power by causal SNP
power_per_snp <- mapply(
  get_gwas_power,
  n     = PARAMS$n_samples,
  beta  = abs(beta_causal),
  maf   = maf_causal,
  MoreArgs = list(alpha = alpha_gwas)
)

# output 

mean_power    <- mean(power_per_snp)
expected_hits <- sum(power_per_snp)

output <- list(
  maf_causal     = maf_causal,
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

# 3. Save data

# RDS format
# we create the folder if it does not exist (for security reasons).
dir_stage1 <- "data/processed/Stage_0"
if(!dir.exists(dir_stage1)) dir.create(dir_stage1, recursive = TRUE)

output_file <- file.path(dir_stage1, "statistical_power.rds")
saveRDS(output, file = output_file)


# CSV format (only the global output)
# we create the folder if it does not exist (for security reasons).
output_dir <- "output/Stage_0"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(global_df, file.path(output_dir, "statistical_power.csv"), row.names = FALSE)

message(paste0(
  "Datos generados y guardados en: ", output_file,
  " en formato RDS y en: ", output_dir,
  " en formato CSV (solo datos globales)"
))
