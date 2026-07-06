library(gamlss)
library(parSim)
library(dplyr)
library(ggplot2)
library(gridExtra) 

# config
sim_results <- parSim(
  n = c(5, 10, 20, 30, 50),  
  mu = c(100),
  sigma = c(0.5, 1.0),       
  
  reps = 100,                # Repeticiones externas
  write = FALSE,
  nCores = 1,
  bar = TRUE,
  
  B_boot = 500,              # Bootstrap interno
  
  expression = {
    y <- rBS6(n = n, mu = mu, sigma = sigma)
    
    mod <- try(gamlss(y ~ 1, family = BS6, 
                      control = gamlss.control(trace = FALSE)), 
               silent = TRUE)
    
    if (class(mod)[1] == "try-error") {
      
      LCI <- NA; 
      LCS <- NA; 
      Width <- NA; 
      Cubre_real <- NA
      mu_est <- NA
      
    } else {
      mu_est <- fitted(mod, "mu")[1]
      sigma_est <- fitted(mod, "sigma")[1]
      
      # Bootstrap
      boots <- matrix(rBS6(n * B_boot, mu = mu_est, sigma = sigma_est),
                      nrow = B_boot, ncol = n)
      
      # n * B_boot números aleatorios de una vez.
      # mu = mu_est, mu estimado antes
      # sigma = sigma_est, sigma estimado antes
      # nrow = 500 filas
      # ncol = c(5, 10, 20, 30, 50) filas
      
      medias_boot <- rowMeans(boots) # Promedio de cada fila
      
      LCI <- quantile(medias_boot, 0.0027) 
      LCS <- quantile(medias_boot, 0.9973)
      Width <- LCS - LCI
      
      Cubre_real <- (mu >= LCI & mu <= LCS)
    }
    
    list(mu_hat = mu_est, Width = Width, Cubre = Cubre_real)
  }
)




# Resultados
sim_results$case <- as.factor(paste0("Sigma=", sim_results$sigma))

resumen <- sim_results %>% 
  group_by(n, case) %>% 
  summarise(
    # Precisión (Anchura de los limites)
    Ancho_Promedio = mean(Width, na.rm = TRUE),
    # Exactitud (Que tanto se aleja del valor verdadero)
    Sesgo_Porc = mean((mu_hat - 100)/100 * 100, na.rm = TRUE),
    # Confiabilidad (Cobertura del valor verdadero)
    Cobertura_Pct = mean(Cubre, na.rm = TRUE) * 100
  )

# ---------------------------------- Graficas -------------------------

# Gráfica A: Ancho
p1 <- ggplot(resumen, aes(x = n, y = Ancho_Promedio, color = case)) +
  geom_line(size = 1.2) + geom_point(size = 3) +
  labs(title = "Precisión (Ancho de Límites)", 
       subtitle = "Menor ancho = Mayor sensibilidad",
       y = "Ancho", x = "") +
  theme_bw() + theme(legend.position = "none")


# Gráfica B: Sesgo
p2 <- ggplot(resumen, aes(x = n, y = Sesgo_Porc, color = case)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(size = 1.2) + geom_point(size = 3) +
  labs(title = "Exactitud de estimación", 
       subtitle = "Cercanía al valor real (0%)",
       y = "Sesgo relativo (%)", x = "") +
  theme_bw() + theme(legend.position = "none")


# Gráfica C: Cobertura
p3 <- ggplot(resumen, aes(x = n, y = Cobertura_Pct, color = case)) +
  geom_hline(yintercept = 99.73, linetype = "dashed", color = "red") + # Meta Teórica
  geom_line(size = 1.2) + geom_point(size = 3) +
  scale_y_continuous(limits = c(90, 100)) +
  labs(title = "Confiabilidad (Cobertura)", 
       subtitle = "Meta teórica: 99.73% (Línea roja)",
       y = "Cobertura real (%)", x = "Tamaño de muestra (n)",
       color = "Variabilidad") +
  theme_bw() + theme(legend.position = "bottom")

# Unión
grid.arrange(p1, p2, p3, nrow = 3)

