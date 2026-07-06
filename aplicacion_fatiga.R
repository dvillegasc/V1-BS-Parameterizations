library(ggplot2)
library(gamlss)
# library(RelDists) 

# ==============================================================================
# 1. DATOS
# ==============================================================================
fatiga_data <- c(70, 90, 96, 97, 99, 100, 103, 104, 104, 105, 107, 108, 108, 108, 109, 
                 109, 112, 112, 113, 114, 114, 114, 116, 119, 120, 120, 120, 121, 121, 
                 123, 124, 124, 124, 124, 124, 128, 128, 129, 129, 130, 130, 130, 131, 
                 131, 131, 131, 131, 132, 132, 132, 133, 134, 134, 134, 134, 134, 136, 
                 136, 137, 138, 138, 138, 139, 139, 141, 141, 142, 142, 142, 142, 142, 
                 142, 144, 144, 145, 146, 148, 148, 149, 151, 151, 152, 155, 156, 157, 
                 157, 157, 157, 158, 159, 162, 163, 163, 164, 166, 166, 168, 170, 174, 
                 196, 212)

# ==============================================================================
# 2. MODELO Y EXTRACCIÓN LIMPIA DE PARÁMETROS
# ==============================================================================
mod_fatiga <- gamlss(fatiga_data ~ 1, family = BS, trace = FALSE)

beta_hat  <- fitted(mod_fatiga, "mu")[1]     
alpha_hat <- fitted(mod_fatiga, "sigma")[1]  

val_a <- sprintf("%.3f", alpha_hat)
val_b <- sprintf("%.1f", beta_hat)

# Formulación matemática robusta para R
lbl_a <- sprintf("bold('Forma')~'('*hat(alpha)*')' == '%s'", val_a)
lbl_b <- sprintf("bold('Escala')~'('*hat(beta)*')' == '%s'", val_b)

# ==============================================================================
# 3. DISEÑO GRÁFICO (ROJO CLÁSICO + ESTRUCTURA LIMPIA)
# ==============================================================================
p_estimacion <- ggplot(data.frame(x = fatiga_data), aes(x = x)) +
  
  # Histograma rojo pastel original
  geom_histogram(aes(y = after_stat(density)), fill = "#FF9999", color = "white", 
                 bins = 18, alpha = 0.9) +
  
  # Curva azul oscuro contrastante
  stat_function(fun = dBS, args = list(mu = beta_hat, sigma = alpha_hat), 
                color = "#2c3e50", linewidth = 1.5) +
  
  # Aumentamos el techo de la gráfica sutilmente para acomodar el texto sin que se corte
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  # Textos inyectados directamente como anotaciones matemáticas (sin cajas de borde)
  annotate("text", x = 165, y = 0.026, label = lbl_a, parse = TRUE, 
           color = "#2c3e50", size = 5.5, hjust = 0) +
  
  annotate("text", x = 165, y = 0.024, label = lbl_b, parse = TRUE, 
           color = "#2c3e50", size = 5.5, hjust = 0) +
  
  labs(title = "Modelación del tiempo de fatiga de láminas de aluminio",
       x = "Ciclos de vibración (miles)", y = "Densidad de probabilidad") +
  
  # Tema limpio sin distractores
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank() # Mantiene la limpieza profesional sin rayas verticales
  )

# Exportación directa y segura
ggsave("aplicacion_fatiga.pdf", plot = p_estimacion, width = 7, height = 5, units = "in", device = cairo_pdf)

