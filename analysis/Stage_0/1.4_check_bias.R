# analysis/Stage_0/1.4_check_bias.R


# 1. load libraries and results
library(dplyr)
gwas_results <- readRDS("data/processed/Stage_0/gwas_results.rds")

# 2. recalculate FDR and Filter Winners
# We adjust p-values
gwas_results$qval <- p.adjust(gwas_results$pval, method = "fdr")

# We keep the ones that are significant (Winners) and also real (Causal)
winners <- subset(gwas_results, qval < 0.05 & is_causal)

if(nrow(winners) == 0) stop("No hay winners para analizar.")

# 3. calculate Bias and Inflation Percentage
winners <- winners %>%
  mutate(
    # absolute difference
    bias = beta_hat - true_beta,
    
    # percentage error relative to actual value
    inflacion_pct = round(((beta_hat - true_beta) / true_beta) * 100, 2),
    
    # label for quick understanding
    tipo_error = ifelse(abs(beta_hat) > abs(true_beta), "Inflado (Curse)", "Subestimado")
  )

# 4. generate clean report
# select the most important columns to view
reporte <- winners %>%
  select(snp_id, true_beta, beta_hat, inflacion_pct, tipo_error) %>%
  arrange(desc(abs(inflacion_pct))) # sort by most exaggerated

message("--- REPORTE DETALLADO DE WINNER'S CURSE ---")
message("Total de Causales detectados: ", nrow(winners))

# case count
print(table(winners$tipo_error))

# average inflation (only for those inflated to be fairer to Curse's concept, 
# or for all to see the general bias. Here we give the overall average of the absolute value)
mean_inflation <- mean(abs(winners$inflacion_pct))
message("\nError porcentual promedio (en valor absoluto): ", round(mean_inflation, 2), "%")

message("\n--- TABLA DE DATOS ---")
print(reporte)

# save data
output_dir <- "output/Stage_0"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write.csv(reporte, file.path(output_dir, "winners_curse_report.csv"), row.names = FALSE)
message("\nReporte guardado en: ", file.path(output_dir, "winners_curse_report.csv"))
