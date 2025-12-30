# analysis/Stage_0/1.3_plot_results.R

# This script allows you to visualize the data using a graph.This 
# makes it easier to visualize the winners 
# (SNPs that have a greater effect than they should)

# 1. load libraries and results
library(ggplot2)
library(dplyr)

results_file <- "data/processed/Stage_0/gwas_results.rds"
if(!file.exists(results_file)) stop("Faltan los resultados del paso 02.")

res <- readRDS(results_file)

# we use FDR to adjust for multiple testing
# we calculate the adjusted p-value for each variant
res$qval <- p.adjust(res$pval, method = "fdr")

# we define ‘Significant’ if the FDR (p-value) is < 0.05
res$is_significant <- res$qval < 0.05

# 3. categorise the points for the graph
res <- res %>%
  mutate(category = case_when(
    is_significant & is_causal ~ "Causal Detectado (Winner)",
    is_significant & !is_causal ~ "Falso Positivo",
    !is_significant & is_causal ~ "Causal No Detectado",
    TRUE ~ "Ruido de Fondo"
  ))

# we calculate how many Winners we have to put in the subtitle
n_winners <- sum(res$category == "Causal Detectado (Winner)")

# 4. plot: true betas (x) vs estimated betas (y)
p1 <- ggplot(res, aes(x = true_beta, y = beta_hat)) +
  # first, we paint the noise in grey (background layer)
  geom_point(data = subset(res, category == "Ruido de Fondo"), 
             color = "grey90", alpha = 0.3, size = 1) +
  # We colour the lost causal SNPs in blue
  geom_point(data = subset(res, category == "Causal No Detectado"), 
             color = "blue", alpha = 0.5, size = 2) +
  # pintamos los falsos positivos en naranja
  geom_point(data = subset(res, category == "Falso Positivo"), 
             color = "orange", alpha = 0.6, size = 2) +
  # We paint the winners in red (top coat)
  geom_point(data = subset(res, category == "Causal Detectado (Winner)"), 
             color = "red", size = 2.5) +
  
  # identity line (the theoretical truth)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  
  labs(
    title = "Visualización del Winner's Curse (Corrección FDR < 0.05)",
    subtitle = paste("Los puntos rojos tienden a agruparse sobre la línea discontinua.\nWinners detectados:", n_winners),
    x = "Efecto REAL (True Beta)",
    y = "Efecto ESTIMADO (GWAS Beta)"
  ) +
  theme_minimal() +
  coord_fixed() # So that the axes have the same visual scale

# 5. save graph
output_dir <- "output/figures/Stage_0"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "winners_curse_plot.png"), plot = p1, width = 8, height = 7)

print(p1)
message("Gráfico generado. Guardado en output/figures/Stage_0/winners_curse_plot.png")