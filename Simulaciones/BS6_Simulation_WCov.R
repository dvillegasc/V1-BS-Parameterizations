# To perform the simulation -----------------------------------------------
if (!dir.exists("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls")) {
  dir.create("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls")
}


library("parSim")

parSim(
  ### SIMULATION CONDITIONS
  n = c(200, 600, 1000, 1400),
  mu = c(1, 2),
  sigma = c(0.1, 1),
  
  reps = 1000,                                # repetitions
  write = TRUE,                               # Writing to a file
  name = "C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls/BS6_Sim_WCov",  # Name of the file
  nCores = 6,                                 # Number of cores to use
  progressbar = TRUE,               # Progress bar
  export = c("dBS", "pBS", "qBS", "rBS", "hBS", "dBS6", "pBS6", "qBS6", "rBS6", "hBS6", "BS6"),
  
  expression = {
    library(gamlss)
    library(gamlss2)
    
    # True parameter values
    y <- rBS6(n=n, mu, sigma)
    
    f   <- y ~ 1
    mod <- try(suppressMessages(gamlss2(f, family=BS6)), silent = TRUE)
    
    if (class(mod)[1] == "try-error") {
      mu_hat    <- NA
      sigma_hat <- NA
    }
    else {
      mu_hat <- exp(coef(mod, what="mu")["(Intercept)"])
      sigma_hat <- exp(coef(mod, what="sigma")["(Intercept)"])
    }
    
    # Results list:
    Results <- list(
      mu_hat = mu_hat,
      sigma_hat = sigma_hat
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------

archivos <- list.files(pattern = "^BS6_Sim_WCov.*\\.txt$", 
                       path="C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls",
                       full.names = TRUE)

archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)

datos$case <- with(datos, 
                   ifelse(mu==1 & sigma==0.1, 1, 
                          ifelse(mu==1 & sigma==1, 2,
                                 ifelse(mu==2 & sigma==0.1, 3, 4))))
datos$case <- as.factor(datos$case)

# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.03 # percentage of values to be trimmed

dat <- datos %>% group_by(n, mu, case) %>% 
  summarise(nobs = n(),
            mean_mu = mean(mu_hat, trim=trim, na.rm=TRUE),
            mean_si = mean(sigma_hat, trim=trim, na.rm=TRUE),
            mse_mu = mean((mu_hat - mu)^2, trim=trim, na.rm=TRUE),
            mse_si = mean((sigma_hat - sigma)^2, trim=trim, na.rm=TRUE),
            bias_mu = mean(mu_hat-mu, trim=trim, na.rm=TRUE),
            bias_si = mean(sigma_hat-sigma, trim=trim, na.rm=TRUE),
  )

dat


# Plots

if (!dir.exists("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs")) {
  dir.create("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs")
}


p1 <- ggplot(dat, aes(x=n, y=bias_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("Bias for ", mu)))
p1

p2 <- ggplot(dat, aes(x=n, y=bias_si, colour=case)) +
  geom_line() + 
  ylab(expression(paste("Bias for ", sigma)))
p2

p1_final <- p1 + theme_bw(base_size = 13)
p2_final <- p2 + theme_bw(base_size = 13)


ggsave(filename="C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs/bias_BS6_Sim_WCov.pdf", width=12, height=6,
       plot=p1_final+p2_final)


p3 <- ggplot(dat, aes(x=n, y=mse_mu, colour=case)) +
  geom_line() + 
  ylab(expression(paste("MSE for ", mu)))
p3

p4 <- ggplot(dat, aes(x=n, y=mse_si, colour=case)) +
  geom_line() + 
  ylab(expression(paste("MSE for ", sigma)))
p4

p3_final <- p3 + theme_bw(base_size = 13)
p4_final <- p4 + theme_bw(base_size = 13)

ggsave(filename="C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs/mse_BS6_Sim_WCov.pdf", width=12, height=6,
       plot=p3_final+p4_final)

