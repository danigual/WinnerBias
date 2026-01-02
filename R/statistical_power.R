# R/statistical_power.R

#' Calculate the statistical power of a study before conducting it.
#' We must enter the parameters we expect to obtain.
#'
#' @param n sample size
#' @param beta expected effect
#' @param maf expected minor allele frequency
#' @param alpha significance threshold
#'
#' @return statistical power, a number between 0 and 1 
#' (the higher the number, the greater the statistical power)
#' @export

get_gwas_power <- function(n, beta, maf= 0.3, alpha= 5e-8) {
  var_g <- 2 * maf * (1 - maf)
  
  # Here we square (beta^2), so the negative sign disappears.
  # A beta of -2 has the same effect as a beta of 2.
  r2 <- (beta^2 * var_g)
  # Protección frente a valores de r2 mayores de 1 ya que eso no puede existir
  if(r2 >= 1) return(1)
  # Calculo de valor de no centralidad
  ncp <- n * r2 / (1 - r2)
  # Chi cuadrado
  chi_crit <- qchisq(1 - alpha, df = 1)
  # Valor de la potencia 
  power <- 1 - pchisq(chi_crit, df = 1, ncp = ncp)
  
  return(power)
}