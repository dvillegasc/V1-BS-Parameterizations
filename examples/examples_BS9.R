#Verificacion de la dBS9

integrate(dBS9, lower=0, upper=999, mu=10, sigma=1.5)

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual <- function(y, mu, sigma) {
  a0 <- sqrt(2 * (sigma - 1))
  b0 <- mu / sigma
  db_dm <- 1 / sigma
  
  term1 <- (1 / (y + b0)) * db_dm
  term2 <- -1 / (2 * b0) * db_dm
  term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  return(term1 + term2 + term3)
}

dldd_manual <- function(y, mu, sigma) { 
  a0 <- sqrt(2 * (sigma - 1))
  b0 <- mu / sigma
  
  da_ds <- 1 / a0
  db_ds <- -b0 / sigma
  
  term1 <- (-1 / a0) * da_ds
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  term3 <- (1 / (y + b0)) * db_ds
  term4 <- (-1 / (2 * b0)) * db_ds
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_ds
  
  return(term1 + term2 + term3 + term4 + term5)
}


# Computacionales

dldm_compu <- function(y, mu, sigma) {
  dm <- gamlss::numeric.deriv(expr = dBS9(y, mu, sigma, log = TRUE), theta = "mu", delta = 1e-04)
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}

dldd_compu <- function(y, mu, sigma) {
  ds <- gamlss::numeric.deriv(expr = dBS9(y, mu, sigma, log = TRUE), theta = "sigma", delta = 1e-04)
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
true_mu <- 0.1
true_si <- 50

y <- rBS9(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS9,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------------------

curve(dBS9(x, mu = 2, sigma= 1.1), from = 0.0000001, to = 3,
      ylim = c(0, 0.92),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS9(x, mu = 2, sigma= 1.3), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 2, sigma= 1.5), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS9(x, mu = 2, sigma= 1.7), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS9(x, mu = 2, sigma= 1.9), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 2, sigma= 2.1), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("Φ = 1.1", "Φ = 1.3", "Φ = 1.5", "Φ = 1.7", "Φ = 1.9", "Φ = 2.1"))


#------------------------ Grafica 2 ------------------------------------


curve(dBS9(x, mu = 0.75, sigma= 1.5), from = 0.0000001, to = 2,
      ylim = c(0, 1.5),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS9(x, mu = 1, sigma= 1.5), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 1.25, sigma= 1.5), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS9(x, mu = 1.5, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS9(x, mu = 1.75, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS9(x, mu = 2, sigma= 1.5), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("μ = 0.75", "μ = 1","μ = 1.25", "μ = 1.5", "μ = 1.75", "μ = 2"))

#-------------------------------- Grafica 3 --------------------------------

varBS9 <- function(mu, sigma) {
  return (((mu^2) * (sigma - 1) * (5*sigma - 3)) / sigma^2)
}


mu <- 2
sigma <- seq(1.1, 250, length.out = 1000) 

var_values <- varBS9(mu = mu, sigma = sigma)


plot(sigma, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(0.1, 20),      
     xlim = c(0.1, 250),      
     xlab = expression(phi),
     ylab = "Var[T]",      
     main = "",            
     las = 1)             


legend(x= 100, y= 10,
       bty="n",
       cex = 0.9,       
       legend = c("μ = 2"))
