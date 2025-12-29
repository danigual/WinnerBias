
# analysis/Stage_1/1.0_power_calculation.R

#Este script calcula y visualiza la potencia estadística de un experimento de GWAS (Genome-Wide Association Study) 
#imulado, en función de varios parámetros clave, para ayudarte a planificar el estudio y entender 
#qué tan probable es detectar efectos genéticos bajo distintas condiciones

#' We calculate the statistical power of stage_1 before running the experiment

# load the required libraries and functions

library(ggplot2)
library(dplyr)

source("R/statistical_power.R")

# 1. Stage setup
S_deseado <- 40
n_causal <- 50

PARAMS <- list(
  n_samples = 5000,   # Número de pacientes
  n_snps    = 10000,  # Número de variantes genéticas
  n_causal  = 50,     # Número de variantes que realmente funcionan
  h2        = 0.5,    # Heredabilidad (50% genética, 50% ambiente)
  seed      = 42,
  S_deseado   = S_deseado,            # p.ej. quieres 40 causales significativos
  P_target    = S_deseado / n_causal   # potencia objetivo por SNP causal
)

# calcular el beta para saber poder estadistico

source("R/statistical_power.R")

get_beta_for_power <- function(target_power, n, maf, alpha) {
  f <- function(beta) {
    get_gwas_power(n = n, beta = beta, maf = maf, alpha = alpha) - target_power
  }
  uniroot(f, interval = c(0, 5))$root
}

# Variables para simulación
N <- PARAMS$n_samples
maf <- 0.1  # ejemplo
alpha <- 5e-8

beta_needed <- get_beta_for_power(
  target_power = PARAMS$P_target,
  n = N,
  maf = maf,
  alpha = alpha
)

print (beta_needed)


# 3. Construir tabla

# Valores para probar en el gráfico
betas_to_test <- seq(0, 4, by = 0.1)
mafs_to_test <- c(0.01, 0.05, 0.1, 0.3)
N_simulado <- N
alpha_fdr_approx <- alpha

#Table
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
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") + 
  # Añadimos un rectángulo sombreado para mostrar la "Zona de Peligro"
  # (Donde la potencia es baja)
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 1, alpha = 0.1, fill = "red") +
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
ggsave("output/Stage_1/power_curves.png", p_power, width = 8, height = 6)

print(p_power)
message("Gráfico de potencia cubriendo el rango de 0 a 4.")
