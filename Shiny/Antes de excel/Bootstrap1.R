library(shiny)
library(ggplot2)
library(dplyr)
library(gamlss)

# Se deben tener cargadas dBS6 y BS6


# --- Interfaz ---
ui <- fluidPage(
  titlePanel("Carta de control BS6 (Media) - Bootstrap"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Parámetros reales del proceso"),
      numericInput("target_mu", "Media real (µ):", value = 100, step = 10),
      numericInput("target_sigma", "Variabilidad real (Sigma):", value = 0.5, step = 0.1),
      
      hr(),
      h4("2. Diseño de la carta"),
      numericInput("n_sample", "Tamaño de subgrupo (n):", value = 5, min = 1),
      numericInput("k_groups", "Cantidad de subgrupos (m):", value = 30, min = 20),
      
      hr(),
      h4("3. Configuración bootstrap"),
      numericInput("B_boot", "Réplicas bootstrap:", value = 10000, min = 1000),
      sliderInput("conf_level", "Nivel de confianza (%):", 
                  min = 90, max = 99.9, value = 99.73),
      
      hr(),
      actionButton("run_sim", "Generar carta", class = "btn-primary btn-lg")
    ),
    
    mainPanel(
      plotOutput("control_chart", height = "500px"),
      verbatimTextOutput("resumen_estadistico")
    )
  )
)

# --- SERVER (LÓGICA) ---
server <- function(input, output) {
  
  # Calculamos solo al pulsar el botón
  resultados <- eventReactive(input$run_sim, {
    
    # 1. Inputs
    n <- input$n_sample
    k <- input$k_groups
    mu_real <- input$target_mu
    sigma_real <- input$target_sigma
    B <- input$B_boot
    
    # Generar datos (Simulacion)
    datos_matriz <- matrix(NA, nrow = k, ncol = n)
    
    for(i in 1:k) {
      datos_matriz[i, ] <- rBS6(n, mu = mu_real, sigma = sigma_real)
    }
    
    # Puntos a graficar
    promedios_subgrupos <- rowMeans(datos_matriz)
    
    #  Estimacion
    todos_los_datos <- as.vector(datos_matriz)
    
    # Ajuste con GAMLSS
    modelo <- try(gamlss(todos_los_datos ~ 1, family = BS6, 
                         control = gamlss.control(trace = FALSE)), 
                  silent = TRUE)
    
    if (class(modelo)[1] == "try-error") {
      mu_est <- mean(todos_los_datos)
      sigma_est <- sd(todos_los_datos) 
    } else {
      mu_est <- fitted(modelo, "mu")[1]
      sigma_est <- fitted(modelo, "sigma")[1]
    }
    
    #  Bootstrap
    boots <- matrix(rBS6(n * B, mu = mu_est, sigma = sigma_est),
                    nrow = B, ncol = n)
    
    medias_boot <- rowMeans(boots)
    
    #  Calcular limites
    alpha <- 1 - (input$conf_level / 100)
    LCI <- quantile(medias_boot, probs = alpha / 2)
    LCS <- quantile(medias_boot, probs = 1 - (alpha / 2))
    LC  <- mu_est 
    
    list(
      datos_grafica = data.frame(Grupo = 1:k, Media = promedios_subgrupos),
      limites = c(LCI = LCI, LC = LC, LCS = LCS),
      params_est = c(mu = mu_est, sigma = sigma_est)
    )
  })
  
  # --- Grafica ---
  output$control_chart <- renderPlot({
    res <- resultados()
    df <- res$datos_grafica
    lims <- res$limites
    
    # Puntos fuera de control
    df$Color <- ifelse(df$Media > lims[3] | df$Media < lims[1], "Fuera", "Dentro")
    
    ggplot(df, aes(x = Grupo, y = Media)) +
      # Límites 
      geom_hline(yintercept = lims[1], color = "red", linetype = "dashed", size = 0.8) +
      geom_hline(yintercept = lims[3], color = "red", linetype = "dashed", size = 0.8) +
      # Línea Central
      geom_hline(yintercept = lims[2], color = "darkgreen", size = 0.8) +
      
      # Puntos y Conectores
      geom_line(color = "gray50") +
      geom_point(aes(color = Color), size = 2.5) +
      scale_color_manual(values = c("Dentro" = "black", "Fuera" = "red")) +
      
      # Etiquetas
      annotate("text", x = max(df$Grupo)+0.5, y = lims[3], label = "LCS", color = "red", hjust = 0) +
      annotate("text", x = max(df$Grupo)+0.5, y = lims[1], label = "LCI", color = "red", hjust = 0) +
      annotate("text", x = max(df$Grupo)+0.5, y = lims[2], label = "LC", color = "darkgreen", hjust = 0) +
      
      labs(title = "Carta de control X-barra (Bootstrap)", 
           subtitle = paste("Límites basados en", input$B_boot, "réplicas simuladas"),
           y = "Media muestral", x = "Subgrupo") +
      
      theme_bw() + 
      theme(legend.position = "none") +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) 
  })
  
  # --- Texto ---
  output$resumen_estadistico <- renderPrint({
    res <- resultados()
    
    cat("--- RESULTADOS ---\n")
    cat(sprintf("Mu estimado (Fase I):    %.4f\n", res$params_est['mu']))
    cat(sprintf("Sigma estimado (Fase I): %.4f\n", res$params_est['sigma']))
    cat("------------------\n")
    cat("Límites calculados:\n")
    print(res$limites)
  })
}

shinyApp(ui, server)
