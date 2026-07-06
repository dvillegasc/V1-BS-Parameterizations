library(shiny)
library(ggplot2)
library(dplyr)
library(gamlss)
library(readxl)

# Asegúrate de tener cargadas tus funciones rBS6 y dBS6 antes de correr la app
# source("tus_funciones_BS6.R") 

# --- UI ---
ui <- fluidPage(
  titlePanel("Carta de control BS (P4) - Fase I"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Carga de datos"),
      
      # Botón de Ayuda
      actionLink("ayuda_formato", " ¿Cómo debe ser mi archivo?", icon = icon("question-circle")),
      br(), br(),
      
      fileInput("file1", "Sube tu archivo (.xlsx o .csv)",
                multiple = FALSE,
                accept = c(".csv", ".xlsx")),
      
      checkboxInput("header", "El archivo tiene títulos en la primera fila", TRUE),
      
      # Panel dinámico para CSV (Siempre visible si es CSV)
      uiOutput("opciones_csv"), 
      
      # Panel de Transformación (Solo aparece si detectamos 1 sola columna)
      uiOutput("panel_agrupacion"),
      
      hr(),
      h4("2. Configuración"),
      numericInput("B_boot", "Réplicas bootstrap:", value = 10000, min = 1000),
      sliderInput("conf_level", "Nivel de confianza (%):", 
                  min = 90, max = 99.9, value = 99.73),
      
      hr(),
      actionButton("analyze", "Analizar datos", class = "btn-success btn-lg")
    ),
    
    mainPanel(
      # Mensaje de estado
      uiOutput("mensaje_estado"),
      
      tabsetPanel(
        tabPanel("Carta de control", 
                 br(),
                 plotOutput("control_chart", height = "500px"),
                 verbatimTextOutput("resumen_params")),
        
        tabPanel("Vista previa de datos", 
                 br(),
                 h4(textOutput("info_dimensiones")),
                 helpText("Nota: Se muestran máximo las primeras 15 filas y 10 columnas."),
                 br(),
                 tableOutput("tabla_datos"))
      )
    )
  )
)

# --- SERVER ---
server <- function(input, output, session) {
  
  # 1. Ventana de Ayuda
  observeEvent(input$ayuda_formato, {
    showModal(modalDialog(
      title = "Formatos aceptados",
      h4("Opción A: Subgrupos (Estándar)"),
      p("El archivo tiene varias columnas (Muestra1, Muestra2...). Cada fila es un subgrupo."),
      h4("Opción B: Columna única (Individuales o para Agrupar)"),
      p("El archivo tiene una sola columna de datos. El programa le dará la opción de cortar esta columna en subgrupos automáticamente (ej. de 5 en 5)."),
      easyClose = TRUE,
      footer = modalButton("Entendido")
    ))
  })
  
  # 2. Opciones CSV (Estándar Latino)
  output$opciones_csv <- renderUI({
    req(input$file1)
    ext <- tools::file_ext(input$file1$name)
    
    if (ext == "csv") {
      tagList(
        wellPanel(
          h5("Configuración CSV"),
          radioButtons("sep", "Separador de columnas:",
                       choices = c("Punto y coma (;)" = ";", 
                                   "Coma (,)" = ",", 
                                   "Tabulación" = "\t"),
                       selected = ";"),
          radioButtons("dec", "Separador decimal:",
                       choices = c("Coma (,)" = ",", 
                                   "Punto (.)" = "."),
                       selected = ",") 
        )
      )
    } else {
      return(NULL) 
    }
  })
  
  # 3. Lectura de Datos iniciales
  datos_crudos <- reactive({
    req(input$file1)
    archivo <- input$file1
    ext <- tools::file_ext(archivo$name)
    df <- NULL
    
    tryCatch({
      if (ext == "xlsx") {
        df <- read_excel(archivo$datapath, col_names = input$header)
      } else if (ext == "csv") {
        mi_sep <- if(is.null(input$sep)) ";" else input$sep
        mi_dec <- if(is.null(input$dec)) "," else input$dec 
        
        df <- read.csv(archivo$datapath, 
                       header = input$header, 
                       sep = mi_sep,
                       dec = mi_dec)
      }
    }, error = function(e) return(NULL))
    
    if(is.null(df) || ncol(df) == 0) return(NULL)
    df_num <- df %>% select(where(is.numeric))
    if(ncol(df_num) == 0) return(NULL)
    return(df_num)
  })
  
  # 4. Panel UI para Agrupar
  output$panel_agrupacion <- renderUI({
    req(datos_crudos())
    df <- datos_crudos()
    
    if (ncol(df) == 1) {
      wellPanel(
        h5("Transformación de datos"),
        p("Se detectó una sola columna."),
        checkboxInput("activar_agrupacion", "Agrupar datos en subgrupos (n > 1)", value = FALSE),
        
        conditionalPanel(
          condition = "input.activar_agrupacion == true",
          numericInput("n_deseado", "Tamaño de subgrupo (n):", value = 5, min = 2)
        )
      )
    } else {
      return(NULL)
    }
  })
  
  # 5. Lógica de Transformación
  datos_finales <- reactive({
    req(datos_crudos())
    df <- datos_crudos()
    
    if (ncol(df) == 1 && isTRUE(input$activar_agrupacion)) {
      datos_vec <- df[[1]]
      n_size <- input$n_deseado
      total_obs <- length(datos_vec)
      
      residuo <- total_obs %% n_size
      
      if (residuo > 0) {
        datos_vec <- datos_vec[1:(total_obs - residuo)]
      }
      
      matriz_nueva <- matrix(datos_vec, ncol = n_size, byrow = TRUE)
      df_transformado <- as.data.frame(matriz_nueva)
      colnames(df_transformado) <- paste0("Muestra_", 1:n_size)
      
      return(list(data = df_transformado, warning = residuo))
      
    } else {
      return(list(data = df, warning = 0))
    }
  })
  
  # 6. Mensajes de Estado
  output$mensaje_estado <- renderUI({
    if (is.null(input$file1)) return(NULL)
    
    if (is.null(datos_crudos())) {
      div(style = "color: red; font-weight: bold; background-color: #ffe6e6; padding: 10px; border-radius: 5px; margin-top: 10px;",
          "⚠️ Error de lectura: No se encontraron datos numéricos. Verifique el separador decimal.")
    } else {
      res <- datos_finales()
      df <- res$data
      sobrantes <- res$warning
      
      if (sobrantes > 0) {
        div(style = "color: #856404; background-color: #fff3cd; padding: 10px; border: 1px solid #ffeeba; border-radius: 5px; margin-top: 10px;",
            HTML(paste0("<b>⚠️ Advertencia:</b> Se descartaron las últimas <b>", sobrantes, 
                        "</b> filas del archivo original porque no completaban un subgrupo de tamaño ", input$n_deseado, ".")))
      } else {
        div(style = "color: #155724; background-color: #d4edda; padding: 10px; border-radius: 5px; border: 1px solid #c3e6cb; margin-top: 10px;",
            paste("✅ Archivo listo. Se analizarán", nrow(df), "subgrupos (filas) y", ncol(df), "muestras (columnas)."))
      }
    }
  })
  
  # 7. Análisis (CORREGIDO PARA EVITAR NA)
  resultados <- eventReactive(input$analyze, {
    req(datos_finales())
    df_matrix <- datos_finales()$data
    
    n <- ncol(df_matrix)
    m <- nrow(df_matrix)
    B <- input$B_boot
    
    promedios_reales <- rowMeans(df_matrix, na.rm = TRUE)
    
    todos_los_datos <- as.vector(as.matrix(df_matrix))
    todos_los_datos <- todos_los_datos[!is.na(todos_los_datos)]
    
    modelo <- try(gamlss(todos_los_datos ~ 1, family = BS6, 
                         control = gamlss.control(trace = FALSE)), 
                  silent = TRUE)
    
    if (class(modelo)[1] == "try-error") {
      mu_est <- mean(todos_los_datos)
      sigma_est <- sd(todos_los_datos)
      msg <- "Aviso: Estimación simple utilizada (GAMLSS falló)."
    } else {
      # --- CORRECCIÓN CLAVE: as.numeric para limpiar nombres ---
      mu_est <- as.numeric(fitted(modelo, "mu")[1])
      sigma_est <- as.numeric(fitted(modelo, "sigma")[1])
      
      # Doble chequeo por si acaso devuelve NA
      if(is.na(mu_est) || is.na(sigma_est)){
        mu_est <- mean(todos_los_datos)
        sigma_est <- sd(todos_los_datos)
        msg <- "Aviso: Estimación simple utilizada (GAMLSS dio NA)."
      } else {
        msg <- "Estimación con GAMLSS realizada con éxito."
      }
    }
    
    boots <- matrix(rBS6(n * B, mu = mu_est, sigma = sigma_est), nrow = B, ncol = n)
    medias_boot <- rowMeans(boots)
    
    alpha_risk <- 1 - (input$conf_level / 100)
    LCI <- quantile(medias_boot, probs = alpha_risk / 2)
    LCS <- quantile(medias_boot, probs = 1 - (alpha_risk / 2))
    LC  <- mu_est 
    
    list(
      df_grafica = data.frame(ID = 1:m, Media = promedios_reales),
      limites = c(LCI = LCI, LC = LC, LCS = LCS),
      params = c(mu_est, sigma_est), # Guardamos vector simple sin nombres
      info = c(n = n, m = m),
      msg = msg
    )
  })
  
  # 8. Gráfica
  output$control_chart <- renderPlot({
    res <- resultados()
    df <- res$df_grafica
    lims <- res$limites
    m_val <- res$info['m']
    
    df$Color <- ifelse(df$Media > lims[3] | df$Media < lims[1], "Fuera", "Dentro")
    
    if(m_val <= 20) {
      mis_breaks <- 1:m_val 
    } else if (m_val <= 50) {
      mis_breaks <- seq(1, m_val, by = 2) 
    } else {
      mis_breaks <- pretty(1:m_val, n = 20)
    }
    
    ggplot(df, aes(x = ID, y = Media)) +
      geom_hline(yintercept = lims[1], color = "red", linetype = "dashed") +
      geom_hline(yintercept = lims[3], color = "red", linetype = "dashed") +
      geom_hline(yintercept = lims[2], color = "darkgreen") +
      geom_line(color = "gray50") +
      geom_point(aes(color = Color), size = 3) +
      scale_color_manual(values = c("Dentro" = "black", "Fuera" = "red")) +
      scale_x_continuous(breaks = mis_breaks, expand = expansion(mult = c(0.05, 0.1))) +
      scale_y_continuous(n.breaks = 10) +
      labs(title = paste("Carta de control (n =", res$info['n'], ")"),
           subtitle = res$msg,
           y = "Media muestral", x = "Subgrupo (m)") +
      theme_bw() + theme(legend.position = "none") +
      annotate("text", x = max(df$ID)+0.5, y = lims[3], label = "LCS", color="red", hjust=0) +
      annotate("text", x = max(df$ID)+0.5, y = lims[1], label = "LCI", color="red", hjust=0) + 
      annotate("text", x = max(df$ID)+0.5, y = lims[2], label = "LC", color="darkgreen", hjust=0)
  })
  
  # 9. Resumen de Parámetros (CORREGIDO PARA EVITAR NA)
  output$resumen_params <- renderPrint({
    res <- resultados()
    cat("--- DETALLES DEL PROCESO ---\n")
    cat("Cantidad de subgrupos (m):", res$info['m'], "\n")
    cat("Tamaño de subgrupo (n):   ", res$info['n'], "\n")
    cat("\n--- PARÁMETROS ESTIMADOS ---\n")
    # --- CORRECCIÓN CLAVE: Acceso por posición [1] y [2] ---
    cat(sprintf("Media estimada (Mu):   %.4f\n", res$params[1]))
    cat(sprintf("Alpha (Forma/Sigma):   %.4f\n", res$params[2]))
    cat("\n--- LÍMITES DE CONTROL ---\n")
    print(res$limites)
  })
  
  # 10. Info y Tabla
  output$info_dimensiones <- renderText({
    req(datos_finales())
    df <- datos_finales()$data
    paste("Total de datos cargados:", nrow(df), "filas (subgrupos) x", ncol(df), "columnas (muestras).")
  })
  
  output$tabla_datos <- renderTable({
    req(datos_finales())
    df <- datos_finales()$data
    filas_a_mostrar <- min(nrow(df), 15)
    cols_a_mostrar <- min(ncol(df), 10)
    df[1:filas_a_mostrar, 1:cols_a_mostrar]
  })
}

shinyApp(ui, server)