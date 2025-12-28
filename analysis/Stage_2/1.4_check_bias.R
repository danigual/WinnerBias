# analysis/Stage_2/1.4_check_bias.R
library(dplyr)

# 1. Cargar resultados
gwas_results <- readRDS("data/processed/Stage_2/gwas_results.rds")

# 2. Recalcular FDR y Filtrar Winners
# Ajustamos p-valores
gwas_results$qval <- p.adjust(gwas_results$pval, method = "fdr")

# Nos quedamos con los que son significativos (Winners) Y además son reales (Causales)
winners <- subset(gwas_results, qval < 0.05 & is_causal)

if(nrow(winners) == 0) stop("No hay winners para analizar.")

# 3. Calcular el Sesgo (Bias) y el Porcentaje de Inflación
winners <- winners %>%
  mutate(
    # Diferencia absoluta
    bias = beta_hat - true_beta,
    
    # Porcentaje de error respecto al valor real
    inflacion_pct = round(((beta_hat - true_beta) / true_beta) * 100, 2),
    
    # Etiqueta para entenderlo rápido
    tipo_error = ifelse(abs(beta_hat) > abs(true_beta), "Inflado (Curse)", "Subestimado")
  )

# 4. Generar Reporte Limpio
# Seleccionamos las columnas más importantes para ver
reporte <- winners %>%
  select(snp_id, true_beta, beta_hat, inflacion_pct, tipo_error) %>%
  arrange(desc(abs(inflacion_pct))) # Ordenar por los más exagerados

message("--- REPORTE DETALLADO DE WINNER'S CURSE ---")
message("Total de Causales detectados: ", nrow(winners))

# Conteo de casos
print(table(winners$tipo_error))

# Inflación media (solo de los inflados para ser más justos con el concepto de Curse, 
# o de todos para ver el sesgo general. Aquí ponemos la media global del valor absoluto).
mean_inflation <- mean(abs(winners$inflacion_pct))
message("\nError porcentual promedio (en valor absoluto): ", round(mean_inflation, 2), "%")

message("\n--- TABLA DE DATOS ---")
print(reporte)

# Guardar esta tabla en un CSV si quiere usarse luego
output_dir <- "output/Stage_2"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(reporte, file.path(output_dir, "winners_curse_report.csv"), row.names = FALSE)
message("\nReporte guardado en: ", file.path(output_dir, "winners_curse_report.csv"))
