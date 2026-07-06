# To perform the simulation -----------------------------------------------
if (!dir.exists("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls")) {
  dir.create("C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls")
}


library("parSim")

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(0.2 + 0.3 * x1) # 1.4 approximately
  sigma <- exp(-0.2 + 0.4 * x2) #  1  approximately
  y <- rBS4(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

parSim(
  ### SIMULATION CONDITIONS
  
  n = c(200, 600, 1000, 1400),
  
  reps = 1000,                     # repetitions
  write = TRUE,                     # Writing to a file
  name = "C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls/BS4_Sim_Cov",# Name of the file
  nCores = 6,                       # Number of cores to use
  progressbar = TRUE,               # Progress bar
  export = c("gendat", "dBS", "pBS", "qBS", "rBS", "hBS", "dBS4", "pBS4", "qBS4", "rBS4", "hBS4", "BS4"),
  
  expression = {
    library(gamlss2)
    
    # True parameter values
    dat <- gendat(n=n)
    
    f   <- y ~ x1 | x2 
    mod <- try(suppressMessages(
      gamlss2(f, data = dat, family = BS4, 
              control = gamlss2_control(trace = FALSE, eps = 1e-05, maxit = 300),
              optimizer = RS_CG)
    ), silent = TRUE)
    
    if (class(mod)[1] == "try-error") {
      beta_0_hat  <- NA
      beta_1_hat  <- NA
      gamma_0_hat  <- NA
      gamma_1_hat  <- NA
    }
    else {
      coefs_mu <- coef(mod, what="mu")
      coefs_sigma <- coef(mod, what="sigma")
      
      beta_0_hat  <- coefs_mu["(Intercept)"]
      beta_1_hat  <- coefs_mu["x1"]
      gamma_0_hat  <- coefs_sigma["(Intercept)"]
      gamma_1_hat  <- coefs_sigma["x2"]
    }
    
    
    # Results list:
    Results <- list(
      beta_0_hat = beta_0_hat,
      beta_1_hat = beta_1_hat,
      gamma_0_hat = gamma_0_hat,
      gamma_1_hat = gamma_1_hat
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------

archivos <- list.files(pattern = "^BS4_Sim_Cov.*\\.txt$", 
                       path="C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Simuls",
                       full.names = TRUE)
archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)


prop.table(table(datos$error == TRUE))


# To analize the results --------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

trim <- 0.03

dat <- datos %>% group_by(n) %>% 
  summarise(nobs = n(),
            
            bias_b0 = mean(beta_0_hat - (0.2), trim=trim, na.rm=TRUE),
            bias_b1 = mean(beta_1_hat - (0.3), trim=trim, na.rm=TRUE),
            bias_g0 = mean(gamma_0_hat - (-0.2), trim=trim, na.rm=TRUE),
            bias_g1 = mean(gamma_1_hat - (0.4), trim=trim, na.rm=TRUE),
            
            mse_b0 = mean((beta_0_hat - (0.2))^2, trim=trim, na.rm=TRUE),
            mse_b1 = mean((beta_1_hat - (0.3))^2, trim=trim, na.rm=TRUE),
            mse_g0 = mean((gamma_0_hat - (-0.2))^2, trim=trim, na.rm=TRUE),
            mse_g1 = mean((gamma_1_hat - (0.4))^2, trim=trim, na.rm=TRUE)
            
  )

dat


# Legend and colores
leyenda <- c(expression(hat(beta)[0]), 
             expression(hat(beta)[1]), 
             expression(hat(gamma)[0]), 
             expression(hat(gamma)[1]))

colores <- c("#F8766D", "#D39200", "#00BA38", "#619CFF")


d <- pivot_longer(data=dat, 
                  cols=c("bias_b0", "bias_b1", 
                         "bias_g0", "bias_g1"),
                  names_to="Estimator",
                  values_to="value")

# Plots
p1 <- ggplot(d, aes(x=n, y=value, colour=Estimator)) +
  geom_line() + 
  ylab("Bias") + 
  scale_color_manual(labels=leyenda,
                     values=colores)

p1

d <- pivot_longer(data=dat, 
                  cols=c("mse_b0", "mse_b1", 
                         "mse_g0", "mse_g1"),
                  names_to="Estimator",
                  values_to="value")

p2 <- ggplot(d, aes(x=n, y=value, colour=Estimator)) +
  geom_line() + 
  ylab("MSE") + 
  scale_color_manual(labels=leyenda,
                     values=colores)

p2


p1_final <- p1 + theme_bw(base_size = 13)
p2_final <- p2 + theme_bw(base_size = 13)

# Guardar el archivo PDF con las dimensiones correctas
ggsave(filename = "C:/Users/davil/Desktop/BS-Parametrizations/Simulaciones/Figs/bias_mse_BS4_Sim_Cov.pdf", 
       plot = p1_final + p2_final, 
       width = 7.5, 
       height = 3.5, 
       units = "in")


