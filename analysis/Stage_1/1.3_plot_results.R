# analysis/Stage_1/1.3_plot_results.R

# 1. Cargar librerías y resultados
library(ggplot2)
library(dplyr)

results_file <- "data/processed/Stage_1/gwas_results.rds"
if(!file.exists(results_file)) stop("Faltan los resultados del paso 02.")

res <- readRDS(results_file)

# Usamos FDR para ajustar por multiple testing
# Calculamos el q-valor (p-valor ajustado) para cada variante
res$qval <- p.adjust(res$pval, method = "fdr")

# Definimos "Significativo" si el FDR (q-valor) es < 0.05
# (Esperamos menos de un 5% de falsos positivos en este grupo)
res$is_significant <- res$qval < 0.05

# 3. Categorizar los puntos para el gráfico
res <- res %>%
  mutate(category = case_when(
    is_significant & is_causal ~ "Causal Detectado (Winner)",
    is_significant & !is_causal ~ "Falso Positivo",
    !is_significant & is_causal ~ "Causal No Detectado",
    TRUE ~ "Ruido de Fondo"
  ))

# Calculamos cuántos Winners tenemos para ponerlo en el subtítulo
n_winners <- sum(res$category == "Causal Detectado (Winner)")

# 4. GRÁFICO: La Verdad (X) vs La Estimación (Y)
p1 <- ggplot(res, aes(x = true_beta, y = beta_hat)) +
  # Pintamos primero el ruido en gris (capa del fondo)
  geom_point(data = subset(res, category == "Ruido de Fondo"), 
             color = "grey90", alpha = 0.3, size = 1) +
  # Pintamos los causales perdidos en azul
  geom_point(data = subset(res, category == "Causal No Detectado"), 
             color = "blue", alpha = 0.5, size = 2) +
  # Pintamos los falsos positivos en naranja
  geom_point(data = subset(res, category == "Falso Positivo"), 
             color = "orange", alpha = 0.6, size = 2) +
  # Pintamos los WINNERS en ROJO (Capa superior)
  geom_point(data = subset(res, category == "Causal Detectado (Winner)"), 
             color = "red", size = 2.5) +
  
  # Línea de identidad (La verdad teórica)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  
  labs(
    title = "Visualización del Winner's Curse (Corrección FDR < 0.05)",
    subtitle = paste("Los puntos rojos tienden a estar 'por fuera' de la línea discontinua.\nWinners detectados:", n_winners),
    x = "Efecto REAL (True Beta)",
    y = "Efecto ESTIMADO (GWAS Beta)"
  ) +
  theme_minimal() +
  coord_fixed() # Para que los ejes tengan la misma escala visual

# 5. Guardar gráfico
output_dir <- "output/figures/Stage_1"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "winners_curse_plot.png"), plot = p1, width = 8, height = 7)

print(p1)
message("Gráfico generado. Guardado en output/figures/Stage_1/winners_curse_plot.png")