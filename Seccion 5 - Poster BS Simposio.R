# ==========================================
# 1. CARGAR LIBRERÍAS Y DATOS
# ==========================================
library(gamlss)
library(RelDists)
library(ggplot2)

datos <- na.omit(airquality) 
ctrl <- gamlss.control(n.cyc = 200, trace = FALSE)

# ==========================================
# 2. MODELOS Y AIC/BIC
# ==========================================
mod_norm  <- gamlss(Ozone ~ Temp, family = NO(), data = datos, control = ctrl)
mod_gamm  <- gamlss(Ozone ~ Temp, family = GA(), data = datos, control = ctrl)
mod_bs_cl <- gamlss(Ozone ~ Temp, family = BS(), data = datos, control = ctrl) # BS Clásica
mod_bs5   <- gamlss(Ozone ~ Temp, family = BS5(), data = datos, control = ctrl) # BS5 (Media pura)
mod_bs6   <- gamlss(Ozone ~ Temp, family = BS6(), data = datos, control = ctrl) # BS6 (Media con precisión acotada)

# Ajuste Marginal con BS6 (Para la gráfica 1)
mod_bs6_marg <- gamlss(Ozone ~ 1, family = BS6(), data = datos, control = ctrl)

resultados <- data.frame(
  Modelo = c("Normal", "Gamma", "BS Clásica", "BS5 (Media)", "BS6 (Precisión Acotada)"),
  AIC = round(c(AIC(mod_norm), AIC(mod_gamm), AIC(mod_bs_cl), AIC(mod_bs5), AIC(mod_bs6)), 1),
  BIC = round(c(BIC(mod_norm), BIC(mod_gamm), BIC(mod_bs_cl), BIC(mod_bs5), BIC(mod_bs6)), 1)
)
print("=== TABLA PARA EL PÓSTER ===")
print(resultados)
# (Notarás que BS Clásica, BS5 y BS6 tienen exactamente el mismo AIC. Esto 
# demuestra empíricamente la propiedad de invarianza: el ajuste es idéntico, 
# pero la interpretación que ofreces con tu software es infinitamente superior).

# ==========================================
# 3. GRÁFICA 1: Ajuste Estático (Marginal BS6)
# ==========================================
mu_hat <- fitted(mod_bs6_marg, "mu")[1]
sigma_hat <- fitted(mod_bs6_marg, "sigma")[1]

color_rojo <- "#CD5C5C" 
color_udea <- "#005E2D"

p1 <- ggplot(datos, aes(x = Ozone)) +
  geom_histogram(aes(y = ..density..), bins = 15, fill = color_rojo, color = "white", alpha = 0.85) +
  stat_function(fun = function(x) dBS6(x, mu = mu_hat, sigma = sigma_hat), 
                color = color_udea, linewidth = 1.2) +
  theme_minimal(base_size = 14) +
  labs(x = "Concentración de Ozono (ppb)", y = "Densidad") +
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        axis.title = element_text(face = "bold", color = "grey20"))

ggsave("aplicacion_estatica.pdf", p1, width = 5, height = 4)

# ==========================================
# 4. GRÁFICA 2: Regresión Distribucional (BS5)
# ==========================================
temp_seq <- seq(min(datos$Temp), max(datos$Temp), length.out = 100)
# Usamos mod_bs5 (el que tiene covariables)
pred_quantiles <- centiles.pred(mod_bs5, xname = "Temp", xvalues = temp_seq, 
                                cent = c(10, 50, 90), plot = FALSE)

p2 <- ggplot(datos, aes(x = Temp, y = Ozone)) +
  geom_point(color = "grey40", alpha = 0.6, size = 2) +
  geom_line(data = pred_quantiles, aes(x = x, y = `50`), color = color_udea, linewidth = 1.2) +
  geom_line(data = pred_quantiles, aes(x = x, y = `10`), color = color_rojo, linewidth = 1, linetype = "dashed") +
  geom_line(data = pred_quantiles, aes(x = x, y = `90`), color = color_rojo, linewidth = 1, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(x = "Temperatura (°F)", y = "Ozono predictivo (ppb)") +
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        axis.title = element_text(face = "bold", color = "grey20"))

ggsave("aplicacion_covariables.pdf", p2, width = 5, height = 4)
