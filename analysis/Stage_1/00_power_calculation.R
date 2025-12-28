# analysis/00_power_calculation.R

library(ggplot2)
library(dplyr)

# 1. Función de Potencia 
get_gwas_power <- function(n, beta, maf, alpha) {
  var_g <- 2 * maf * (1 - maf)
  
  # ¡OJO! Aquí elevamos al cuadrado (beta^2), por lo que el signo negativo desaparece.
  # Un beta de -2 tiene el mismo R2 que un beta de 2.
  r2 <- (beta^2 * var_g)
  
  if(r2 >= 1) return(1)
  
  ncp <- n * r2 / (1 - r2)
  chi_crit <- qchisq(1 - alpha, df = 1)
  power <- 1 - pchisq(chi_crit, df = 1, ncp = ncp)
  
  return(power)
}

# 2. Configuración del escenario
N_simulado <- 5000
alpha_fdr_approx <- 1e-4 

# Como la potencia es simétrica, graficamos la MAGNITUD desde 0 hasta 4.
# Esto cubre tanto los efectos pequeños (0.1) como los gigantes (3.0).
betas_to_test <- seq(0, 4, by = 0.1) 

mafs_to_test <- c(0.1, 0.3, 0.5)

# 3. Construir tabla
power_data <- expand.grid(
  beta = betas_to_test,
  maf = mafs_to_test
)

power_data$power <- mapply(get_gwas_power, 
                           n = N_simulado, 
                           beta = power_data$beta, 
                           maf = power_data$maf, 
                           alpha = alpha_fdr_approx)

power_data$maf_label <- paste("MAF =", power_data$maf)

# 4. Graficar
p_power <- ggplot(power_data, aes(x = beta, y = power, color = as.factor(maf_label))) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") + 
  
  # Añadimos un rectángulo sombreado para mostrar la "Zona de Peligro"
  # (Donde la potencia es baja)
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 1, 
           alpha = 0.1, fill = "red") +
  annotate("text", x = 0.25, y = 0.5, label = "Zona de\nRiesgo", color = "red", angle = 90) +
  
  labs(
    title = paste("Potencia Estadística (N =", N_simulado, " | FDR aprox.)"),
    subtitle = "El eje X representa la MAGNITUD del efecto (positivo o negativo).\nLos betas > 2 tienen 100% de potencia (se detectan siempre).",
    x = "Magnitud del Efecto (|Beta Real|)",
    y = "Potencia",
    color = "Frecuencia Alelo"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))

# 5. Guardar
if(!dir.exists("output")) dir.create("output")
ggsave("output/power_curves.png", p_power, width = 8, height = 6)

print(p_power)
message("Gráfico de potencia cubriendo el rango de 0 a 4.")
