#Verificacion de la dBS5

integrate(dBS5, lower=0, upper=99, mu=1.7, sigma=2.3)

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual <- function(y, mu, sigma) {
  a0 <- sqrt(2 / sigma)
  b0 <- (sigma * mu) / (sigma + 1)
  db_dm <- b0 / mu
  
  term1 <- (1 / (y + b0)) * db_dm
  term2 <- -1 / (2 * b0) * db_dm
  term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3
  return(result)
}

dldd_manual <- function(y, mu, sigma) { 
  a0 <- sqrt(2 / sigma)
  b0 <- (sigma * mu) / (sigma + 1)
  
  da_ds <- -a0 / (2 * sigma)
  db_ds <- mu / ((sigma + 1)^2)
  
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
    expr = dBS5(y, mu, sigma, log = TRUE), 
    theta = "mu",                           
    delta = 1e-04)
  
  # Extrae el gradiente
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}


dldd_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS5(y, mu, sigma, log = TRUE), 
    theta = "sigma",                        
    delta = 1e-04)
  
  dldd <- as.vector(attr(ds, "gradient"))
  return(dldd)
}


# PRUEBA

y_test     <- c(1, 2, 5, 15)
mu_test    <- 0.7
sigma_test <- 0.75

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

y <- rBS5(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS5,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------------------ Grafica 1 -----------------------------

curve(dBS5(x, mu = 1, sigma= 2), from = 0.0000001, to = 2,
      ylim = c(0, 3),
      col = "black",        
      lwd = 2,              
      las = 1,
      lty = 1,              
      type= "l",
      ylab = "f(t)",      
      xlab = "t")          

curve(dBS5(x, mu = 1, sigma= 5),   add = TRUE, col = "black", lty = 2, lwd = 2) 
curve(dBS5(x, mu = 1, sigma= 10),  add = TRUE, col = "black", lty = 3, lwd = 2) 
curve(dBS5(x, mu = 1, sigma= 25),  add = TRUE, col = "gray",  lty = 1, lwd = 2) 
curve(dBS5(x, mu = 1, sigma= 50),  add = TRUE, col = "gray",  lty = 2, lwd = 2) 
curve(dBS5(x, mu = 1, sigma= 100), add = TRUE, col = "gray",  lty = 3, lwd = 2) 

legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       lwd = 2,
       bty="n",
       cex = 0.9,        
       legend = c("δ = 2","δ = 5", "δ = 10", "δ = 25", "δ = 50", "δ = 100"))


#---------------------------------- Grafica 2 ------------------------------

curve(dBS5(x, mu = 1, sigma= 5), from = 0.0000001, to = 4,
      ylim = c(0, 1),
      col = "black",        
      lwd = 2,              
      las = 1,
      lty = 1,              
      type= "l",
      ylab = "f(t)",      
      xlab = "t")          

curve(dBS5(x, mu = 1.5, sigma= 5), add = TRUE, col = "black", lty = 2, lwd = 2) 
curve(dBS5(x, mu = 2, sigma= 5),   add = TRUE, col = "black", lty = 3, lwd = 2) 
curve(dBS5(x, mu = 2.5, sigma= 5), add = TRUE, col = "gray",  lty = 1, lwd = 2) 
curve(dBS5(x, mu = 3, sigma= 5),   add = TRUE, col = "gray",  lty = 2, lwd = 2) 
curve(dBS5(x, mu = 3.5, sigma= 5), add = TRUE, col = "gray",  lty = 3, lwd = 2) 

legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       lwd = 2,
       bty="n",
       cex = 0.9,        
       legend = c("μ = 1","μ = 1.5", "μ = 2", "μ = 2.5", "μ = 3", "μ = 3.5"))

#-------------------------------- Grafica 3 --------------------------------

varBS5 <- function(mu, sigma) {
  return ( (mu^2) * (2 * sigma + 5) / ((sigma + 1)^2))
}


mu <- 2
sigma <- seq(0, 20, length.out = 20) 

var_values <- varBS5(mu = mu, sigma = sigma)


plot(sigma, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(0, 20),      
     xlim = c(0, 20),      
     xlab = expression(delta),
     ylab = "Var[T]",      
     main = "",            
     las = 1)             


legend(x= 8, y= 11,
       bty="n",
       cex = 0.9,       
       legend = c("μ = 2"))

