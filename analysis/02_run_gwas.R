# analysis/02_run_gwas.R

# 1. Cargar herramientas
source("R/gwas_engine.R")

# 2. Cargar los datos
input_file <- "data/processed/simulacion_wc.rds"

# Comprobación de seguridad
if(!file.exists(input_file)) stop("¡ALERTA! No encuentro simulacion_wc.rds. ¿Has ejecutado el script 01?")

message("Cargando datos...")
sim_data <- readRDS(input_file)

# 3. EJECUTAR EL GWAS
message("Ejecutando GWAS en 10,000 variantes...")
gwas_results <- run_gwas(X = sim_data$X, y = sim_data$y)

# 4. Añadir la "Hoja de Respuestas" (Para ver el Winner's Curse luego)
gwas_results$true_beta <- sim_data$true_betas
gwas_results$is_causal <- sim_data$true_betas != 0

# 5. Guardar resultados
saveRDS(gwas_results, file = "data/processed/gwas_results.rds")

message("Análisis completado. Resultados guardados en data/processed/gwas_results.rds")



