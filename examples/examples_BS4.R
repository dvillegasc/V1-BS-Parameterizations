#Verificacion de la dBS4

integrate(dBS4, lower=0, upper=99, mu=1.7, sigma=2.3)


#Verificacion de las derivadas


#Derivadas manuales

library(gamlss)


dldm_manual <- function(y, mu, sigma) {
  result <- (y / (sigma + mu * y)) + sigma - (mu * y)
  return(result)
}

dldd_manual <- function(y, sigma, mu) {
  result <- (1 / (sigma + mu * y)) + mu - (sigma / y)
  return(result)
}


#Derivadas computacionales

    dldm_compu <- function(y, mu, sigma) {
      
      dm <- gamlss::numeric.deriv(
        expr = dBS4(y, mu, sigma, log = TRUE), 
        theta = "mu",                          
        delta = 1e-04)
      
      # Extrae el gradiente
      dldm <- as.vector(attr(dm, "gradient"))
      return(dldm)
    }
    
    
    dldd_compu <- function(y, mu, sigma) {
      
      ds <- gamlss::numeric.deriv(
        expr = dBS4(y, mu, sigma, log = TRUE), 
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
true_mu <- 4
true_si <- 5
    
y <- rBS4(n=n, mu=true_mu, sigma=true_si)
    
library(gamlss)
mod <- gamlss(y ~ 1, family=BS4,
              n.cyc = 100)
    
exp(coef(mod, what="mu"))
exp(coef(mod, what="sigma"))


summary(mod)


#---------------------------------- Grafica 1 ---------------------------------

curve(dBS4(x, mu = 1.5, sigma= 5), from = 0.0000001, to = 10,
      ylim = c(0, 1.7),
      col = "black",        
      lwd = 2,              
      las = 1,
      lty = 1,              
      type= "l",
      ylab = "f(t)",      
      xlab = "t")          

curve(dBS4(x, mu = 2, sigma= 5),   add = TRUE, col = "black", lty = 2, lwd = 2) 
curve(dBS4(x, mu = 3, sigma= 5),   add = TRUE, col = "black", lty = 3, lwd = 2) 
curve(dBS4(x, mu = 3.5, sigma= 5), add = TRUE, col = "gray",  lty = 1, lwd = 2) 
curve(dBS4(x, mu = 4, sigma= 5),   add = TRUE, col = "gray",  lty = 2, lwd = 2) 
curve(dBS4(x, mu = 4.5, sigma= 5), add = TRUE, col = "gray",  lty = 3, lwd = 2) 

legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       lwd = 2,
       bty="n",
       cex = 0.9,        
       legend = c(expression(mu[A] == 1.5), 
                  expression(mu[A] == 2), 
                  expression(mu[A] == 3), 
                  expression(mu[A] == 3.5), 
                  expression(mu[A] == 4), 
                  expression(mu[A] == 4.5)))


#----------------------------------- Grafica 2 ----------------------------

curve(dBS4(x, mu = 5, sigma= 0.4), from = 0.0000001, to = 1,
      ylim = c(0, 9),
      col = "black",        
      lwd = 2,              
      las = 1,
      lty = 1,              
      type= "l",
      ylab = "f(t)",      
      xlab = "t")          

curve(dBS4(x, mu = 5, sigma= 0.5), add = TRUE, col = "black", lty = 2, lwd = 2) 
curve(dBS4(x, mu = 5, sigma= 1),   add = TRUE, col = "black", lty = 3, lwd = 2) 
curve(dBS4(x, mu = 5, sigma= 1.5), add = TRUE, col = "gray",  lty = 1, lwd = 2) 
curve(dBS4(x, mu = 5, sigma= 2),   add = TRUE, col = "gray",  lty = 2, lwd = 2) 
curve(dBS4(x, mu = 5, sigma= 3),   add = TRUE, col = "gray",  lty = 3, lwd = 2) 

legend("topright",
       col = c("black", "black", "black", "gray", "gray", "gray"),
       lty = c(1, 2, 3, 1, 2, 3),
       lwd = 2,
       bty="n",
       cex = 0.9,        
       legend = c(expression(lambda[A] == 0.4), 
                  expression(lambda[A] == 0.5),
                  expression(lambda[A] == 1), 
                  expression(lambda[A] == 1.5), 
                  expression(lambda[A] == 2.5), 
                  expression(lambda[A] == 3)))

#-------------------------------- Grafica 3 --------------------------------

varBS4 <- function(mu, sigma) {
  return ( ((sigma*mu) + 5/4) / (mu^4) )
}


sigma <- 2
mu <- seq(1, 4, length.out = 5) 

var_values <- varBS4(mu = mu, sigma = sigma)


plot(mu, var_values, 
     type = "l",           
     lwd = 2,              
     ylim = c(0, 3.5),      
     xlim = c(1, 4),      
     xlab = expression(mu[A]),
     ylab = "Var[T]",      
     main = "",            
     las = 1)             


legend(x= 2.1, y= 1.7,
       bty="n",
       cex = 0.9,       
       legend = expression(lambda[A] == 2))


    