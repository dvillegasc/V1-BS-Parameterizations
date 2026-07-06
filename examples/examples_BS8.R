#Verificacion de la dBS8

integrate(dBS8, lower=0, upper=999, mu=10, sigma=1.5) 

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual = function(y, mu, sigma) {
  a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
  b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
  db_dm <- b0 / (2 * mu)
  
  term1 <- (1 / (y + b0)) * db_dm
  term2 <- -1 / (2 * b0) * db_dm
  term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3
  return(result)
}

dldd_manual = function(y, mu, sigma) { 
  a0 <- (2 * sqrt(sigma - 1)) / sqrt(5)
  b0 <- sqrt(5 * mu) / (2 * sqrt(sigma * (sigma - 1)))
  da_ds <- 1 / sqrt(5 * (sigma - 1))
  db_ds <- -b0 * (2 * sigma - 1) / (2 * sigma * (sigma - 1))
  
  term1 <- (-1 / a0) * da_ds
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  term3 <- (1 / (y + b0)) * db_ds
  term4 <- (-1 / (2 * b0)) * db_ds
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  
  result <- term1 + term2 + term3 + term4 + term5
  return(result)
}


#Derivadas computacionales

dldm_compu <- function(y, mu, sigma) {
  
  dm <- gamlss::numeric.deriv(
    expr = dBS8(y, mu, sigma, log = TRUE), 
    theta = "mu",                           
    delta = 1e-04)
  
  # Extrae el gradiente
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}


dldd_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS8(y, mu, sigma, log = TRUE), 
    theta = "sigma",                        
    delta = 1e-04)
  
  dldd <- as.vector(attr(ds, "gradient"))
  return(dldd)
}


# PRUEBA

y_test     <- c(1, 2, 5, 15)
mu_test    <- 10 
sigma_test <- 1.5 

cat("--- Verificación de dldm (derivada de mu) ---\n")
manual_mu <- dldm_manual(y = y_test, mu = mu_test, sigma = sigma_test)
compu_mu  <- dldm_compu(y = y_test, mu = mu_test, sigma = sigma_test)

print(data.frame(y = y_test, manual = manual_mu, computacional = compu_mu))


cat("\n--- Verificación de dldd (derivada de sigma) ---\n")
manual_sigma <- dldd_manual(y = y_test, mu = mu_test, sigma = sigma_test)
compu_sigma  <- dldd_compu(y = y_test, mu = mu_test, sigma = sigma_test)

print(data.frame(y = y_test, manual = manual_sigma, computacional = compu_sigma))



#-------------------------- Validation Familia ----------------------------------



n <- 1000

# True parameters are:
true_mu <- 1
true_si <- 5

y <- rBS8(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS8,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------------------

curve(dBS8(x, mu = 10, sigma= 1.05), from = 0.0000001, to = 25,
      ylim = c(0, 0.25),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS8(x, mu = 10, sigma= 1.1), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS8(x, mu = 10, sigma= 1.3), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS8(x, mu = 10, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS8(x, mu = 10, sigma= 1.7), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS8(x, mu = 10, sigma= 1.9), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("γ = 1.05","γ = 1.1", "γ = 1.3", "γ = 1.5", "γ = 1.7", "γ = 1.9"))


#------------------------ Grafica 2 ------------------------------------


curve(dBS8(x, mu = 5/10, sigma= 1.5), from = 0.0000001, to = 4, #Requirio una transformacion sobre mu
      ylim = c(0, 0.85),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS8(x, mu = 10/10, sigma= 1.5), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS8(x, mu = 15/10, sigma= 1.5), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS8(x, mu = 20/10, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS8(x, mu = 25/10, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS8(x, mu = 30/10, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("σ² = 5", "σ² = 10","σ² = 15", "σ² = 20", "σ² = 25", "σ² = 30"))


