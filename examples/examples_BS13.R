#Verificacion de la dBS13

integrate(dBS13, lower=0, upper=999, mu=1.5, sigma=1.5)

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual = function(y, mu, sigma) {
  b0 <- mu / sigma
  a0 <- 1 / sqrt(sigma)
  db_dm <- 1 / sigma
  da_dm <- 0 # <--
  
  term3 <- (1 / (y + b0)) * db_dm
  term4 <- (-1 / (2 * b0)) * db_dm
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  result <- term3 + term4 + term5
  return(result)
}

dldd_manual = function(y, mu, sigma) { 
  b0 <- mu / sigma
  a0 <- 1 / sqrt(sigma)
  db_ds <- -b0 / sigma
  da_ds <- -1 / (2 * sigma * sqrt(sigma))
  
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
    expr = dBS13(y, mu, sigma, log = TRUE), 
    theta = "mu",                           
    delta = 1e-04)
  
  # Extrae el gradiente
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}


dldd_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS13(y, mu, sigma, log = TRUE), 
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

y <- rBS13(n=n, mu=true_mu, sigma=true_si)


library(gamlss)
mod <- gamlss(y ~ 1, family=BS13,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------------------

curve(dBS13(x, mu = 5, sigma= 10), from = 0.0000001, to = 0.7,
      ylim = c(0, 16.5),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS13(x, mu = 5, sigma= 15), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS13(x, mu = 5, sigma= 20), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS13(x, mu = 5, sigma= 25), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS13(x, mu = 5, sigma= 30), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS13(x, mu = 5, sigma= 35), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("Ψ = 10", "Ψ = 15", "Ψ = 20", "Ψ = 25", "Ψ = 30", "Ψ = 35"))


#------------------------ Grafica 2 ------------------------------------


curve(dBS13(x, mu = 3, sigma= 10), from = 0, to = 1.5,
      ylim = c(0, 4.5),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS13(x, mu = 3.5, sigma= 10), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS13(x, mu = 4, sigma= 10), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS13(x, mu = 5, sigma= 10), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS13(x, mu = 7, sigma= 10), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS13(x, mu = 10, sigma= 10), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("ω = 3", "ω = 3.5", "ω = 3", "ω = 4", "ω = 7", "ω = 10"))

#-------------------------------- Grafica 3 --------------------------------

varBS13 <- function(mu, sigma) {
  return ((mu^2)/(sigma^3) + (5*mu^2)/(4*sigma^4))
}


mu1 <- 2
sigma1 <- seq(0.20, 0.50, length.out = 100) 

var_values <- varBS13(mu = mu1, sigma = sigma1)


plot(sigma1, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(0, 3600),
     xlim = c(0.20, 0.50),      
     xlab = expression(psi),
     ylab = "Var[T]",      
     main = "",            
     las = 1)      


legend(x= 0.30, y= 1900,
       bty="n",
       cex = 1,       
       legend = c("ω = 2"))




