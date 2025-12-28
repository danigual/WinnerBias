# analysis/01_generate_datasets.R

# 1. Cargar la función generadora

source("R/sim_genetics.R")

# 2. Definir los parámetros del escenario "Winner's Curse"
# Usamos una muestra moderada (1000) con muchos SNPs (10000) para dificultar
# que pasen el filtro de significación, forzando el sesgo.
PARAMS <- list(
  n_samples = 5000,   # Número de pacientes
  n_snps    = 10000,  # Número de variantes genéticas
  n_causal  = 50,     # Número de variantes que realmente funcionan
  h2        = 0.5,    # Heredabilidad (50% genética, 50% ambiente)
  seed      = 42      
)

# 3. Ejecutar la simulación
message("Generando datos simulados (Población Virtual)...")
set.seed(PARAMS$seed)

sim_data <- generate_gwas_data(
  n_samples = PARAMS$n_samples,
  n_snps    = PARAMS$n_snps,
  n_causal  = PARAMS$n_causal,
  h2        = PARAMS$h2
)

# 4. Guardar los datos en el disco
# Creamos la carpeta si no existe (por seguridad)
if(!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)

output_file <- "data/processed/simulacion_wc.rds"
saveRDS(sim_data, file = output_file)

message("Datos generados y guardados en: ", output_file)