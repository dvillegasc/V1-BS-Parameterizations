#Verificacion de la dBS7

integrate(dBS7, lower=0, upper=99, mu=1.7, sigma=2.3)

#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)

dldm_manual = function(y, mu, sigma) {
  a0 <- sigma
  b0 <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
  db_dm <- b0 / (2 * mu)
  
  term1 <- (1 / (y + b0)) * db_dm
  term2 <- -1 / (2 * b0) * db_dm
  term3 <- (1 / (2 * a0^2)) * ((y / b0^2) - (1 / y)) * db_dm
  
  result <- term1 + term2 + term3
  return(result)
}

dldd_manual = function(y, mu, sigma) { 
  a0 <- sigma
  b0 <- (2 * sqrt(mu)) / (sigma * sqrt(4 + 5 * sigma^2))
  da_ds <- 1
  db_ds <- -b0 * ((4 + 10 * sigma^2) / (sigma * (4 + 5 * sigma^2)))
  
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
    expr = dBS7(y, mu, sigma, log = TRUE), 
    theta = "mu",                           
    delta = 1e-04)
  
  # Extrae el gradiente
  dldm <- as.vector(attr(dm, "gradient"))
  return(dldm)
}


dldd_compu <- function(y, mu, sigma) {
  
  ds <- gamlss::numeric.deriv(
    expr = dBS7(y, mu, sigma, log = TRUE), 
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

y <- rBS7(n=n, mu=true_mu, sigma=true_si)

library(gamlss)
mod <- gamlss(y ~ 1, family=BS7,
              n.cyc = 100)

exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)

#------------------------ Grafica 1 ------------------------

curve(dBS7(x, mu = 10, sigma= 0.1), from = 0.0000001, to = 40,
      ylim = c(0, 0.35),
      col = "black",        
      lwd = 2,              
      las = 1,
      lty = 1,            
      type= "l",
      ylab = "f(t)",      
      xlab = "t")          

curve(dBS7(x, mu = 10, sigma= 0.2),  add = TRUE, col = "black", lty = 2, lwd = 2) 
curve(dBS7(x, mu = 10, sigma= 0.4),  add = TRUE, col = "black", lty = 3, lwd = 2) 
curve(dBS7(x, mu = 10, sigma= 0.5),  add = TRUE, col = "gray",  lty = 1, lwd = 2) 
curve(dBS7(x, mu = 10, sigma= 0.75), add = TRUE, col = "gray",  lty = 2, lwd = 2) 
curve(dBS7(x, mu = 10, sigma= 1),    add = TRUE, col = "gray",  lty = 3, lwd = 2)

legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       lwd = 2,
       bty="n",
       cex = 0.9,        
       legend = c("α = 0.1","α = 0.2", "α = 0.4", "α = 0.5", "α = 0.75", "α = 1"))


#------------------------ Grafica 2 ------------------------

curve(dBS7(x, mu = 0.1, sigma = 0.5), from = 0.0000001, to = 5, #El paper se equivoco, la grafica funciona con sigma 0.5 en lugar de 0.1
      ylim = c(0, 2),
      col = "black",        
      lwd = 2,              
      las = 1,
      lty = 1,              
      type= "l",
      ylab = "f(t)",      
      xlab = "t")          

curve(dBS7(x, mu = 0.2, sigma = 0.5),  add = TRUE, col = "black", lty = 2, lwd = 2) 
curve(dBS7(x, mu = 0.5, sigma = 0.5),  add = TRUE, col = "black", lty = 3, lwd = 2) 
curve(dBS7(x, mu = 0.75, sigma = 0.5), add = TRUE, col = "gray",  lty = 1, lwd = 2) 
curve(dBS7(x, mu = 1, sigma = 0.5),    add = TRUE, col = "gray",  lty = 2, lwd = 2) 
curve(dBS7(x, mu = 1.5, sigma = 0.5),  add = TRUE, col = "gray",  lty = 3, lwd = 2) 

legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       lwd = 2,
       bty="n",
       cex = 0.9,        
       legend = c("μ = 0.1", "μ = 1","μ = 1.5", "μ = 2", "μ = 2.5", "μ = 3"))

