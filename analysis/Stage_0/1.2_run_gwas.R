# analysis/Stage_0/1.2_run_gwas.R

# This script allows you to perform the GWAS with the simulated data. 

# 1. load the GWAS engine
source("R/gwas_engine.R")

# 2. load data
input_file <- "data/processed/Stage_0/simulacion_wc.rds"

# security check
if(!file.exists(input_file)) stop("¡ALERTA! No encuentro simulacion_wc.rds. ¿Has ejecutado el script 01?")

message("Cargando datos...")
sim_data <- readRDS(input_file)

# 3. run GWAS
message("Ejecutando GWAS...")
gwas_results <- run_gwas(X = sim_data$X, y = sim_data$y)

# 4. results
gwas_results$true_beta <- sim_data$true_betas
gwas_results$is_causal <- sim_data$true_betas != 0

# 5. save results
output_file <- "data/processed/Stage_0/gwas_results.rds"
saveRDS(gwas_results, file = output_file)
message(paste("Análisis completado. Resultados guardados en", output_file))