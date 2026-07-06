#Verificacion de la dBS12

integrate(dBS12, lower=0, upper=999, mu=1.5, sigma=1.5)

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual = function(y, mu, sigma) {
  b0 <- mu
  a0 <- 1 / sqrt(sigma)
  db_dm <- 1
  da_dm <- 0 # <--
  
  term3 <- (1 / (y + b0)) * db_dm
  term4 <- (-1 / (2 * b0)) * db_dm
  term5 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  result <- term3 + term4 + term5
  return(result)
}

dldd_manual = function(y, mu, sigma) { 
  b0 <- mu
  a0 <- 1 / sqrt(sigma)
  db_ds <- 0 # <--
  da_ds <- -1 / (2 * sigma * sqrt(sigma))
  
  term1 <- (-1 / a0) * da_ds
  term2 <- (1 / a0^3) * ((y / b0) + (b0 / y) - 2) * da_ds
  
  result <- term1 + term2
  return(result)
}


#Derivadas computacionales

dldm_compu <- function(y, mu, sigma) {
  
  dm <- gamlss::numeric.deriv(
    expr = dBS12(y, mu, sigma, log = TRUE), 
    theta = "mu",                           
    delta = 1e-04)
  
  # Extrae el gradiente
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}


dldd_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS12(y, mu, sigma, log = TRUE), 
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

y <- rBS12(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS12,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------------------

curve(dBS12(x, mu = 1, sigma= 0.5), from = 0.0000001, to = 2,
      ylim = c(0, 2.2),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS12(x, mu = 1, sigma= 2), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS12(x, mu = 1, sigma= 5), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS12(x, mu = 1, sigma= 10), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS12(x, mu = 1, sigma= 20), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS12(x, mu = 1, sigma= 30), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("Ψ = 0.5", "Ψ = 2", "Ψ = 5", "Ψ = 10", "Ψ = 20", "Ψ = 30"))


#------------------------ Grafica 2 ------------------------------------


curve(dBS12(x, mu = 0.75, sigma= 25), from = 0.25, to = 2.5,
      ylim = c(0, 2.7),
      col = "black",        
      lwd = 2,              
      las = 1,
      type= "l",
      ylab = "f(t)",      
      xlab = "t",
      main = "")          

curve(dBS12(x, mu = 0.95, sigma= 25), add = TRUE, col = "black", type= "l", lty=2, lwd = 2)

curve(dBS12(x, mu = 1, sigma= 25), add = TRUE, col = "black", type= "l", lty=3, lwd = 2)

curve(dBS12(x, mu = 1.15, sigma= 25), add = TRUE, col = "gray", type= "l", lty=1, lwd = 2)

curve(dBS12(x, mu = 1.25, sigma= 25), add = TRUE, col = "gray", type= "l", lty=2, lwd = 2)

curve(dBS12(x, mu = 1.5, sigma= 25), add = TRUE, col = "gray", type= "l", lty=3, lwd = 2)


legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       bty="n",
       cex = 0.9,        
       legend = c("β = 0.75", "β = 0.95", "β = 1", "β = 1.15", "β = 1.25", "β = 1.5"))

#-------------------------------- Grafica 3 --------------------------------

varBS12 <- function(mu, sigma) {
  return ((mu^2)/(sigma) + (5*mu^2)/(4*sigma^2))
}


mu1 <- 1
sigma1 <- seq(0.20, 0.50, length.out = 100) 

var_values <- varBS12(mu = mu1, sigma = sigma1)


plot(sigma1, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(7, 36),
     xlim = c(0.20, 0.50),      
     xlab = expression(omega),
     ylab = "Var[T]",      
     main = "",            
     las = 1)      


legend(x= 0.31, y= 22,
       bty="n",
       cex = 1,       
       legend = c("β = 1"))
    



