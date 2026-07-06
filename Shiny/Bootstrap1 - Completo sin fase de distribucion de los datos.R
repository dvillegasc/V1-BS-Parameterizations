library(shiny)
library(ggplot2)
library(dplyr)
library(gamlss)
library(readxl)

# --- Interfaz ---

ui <- fluidPage(
  titlePanel("Carta de control BS (P4)"),
  
  sidebarLayout(
    sidebarPanel(
      # --- codicionales ---
      
      # Diseño y ajuste
      conditionalPanel(
        condition = "input.fases_principales == 'Diseño y ajuste'",
        
        h4("1. Carga de datos"),
        # Ayuda formato
        actionLink("ayuda_f1", " Ayuda en formato de carga de datos", icon = icon("question-circle")),
        br(), br(),
        
        fileInput("file1", "Sube tus datos (.xlsx / .csv)", multiple = FALSE, accept = c(".csv", ".xlsx")),
        checkboxInput("header1", "El archivo tiene títulos", TRUE),
        
        # Solo para CSV o columna única
        uiOutput("opciones_csv_1"),      
        uiOutput("panel_agrupacion_1"),  
        
        hr(),
        h4("2. Configuración"),
        numericInput("B_boot", "Réplicas bootstrap:", value = 10000, min = 1000, max = 20000),
        sliderInput("conf_level", "Nivel de confianza (%):", min = 90, max = 99.9, value = 99.73),
        
        hr(),
        actionButton("analyze", "Diseñar carta", class = "btn-success btn-lg"),
        br(), br(),
        actionButton("reset_puntos", "Restaurar puntos", class = "btn-warning btn-sm", icon = icon("refresh")),
        
        hr(style="border-top: 2px solid #333;"),
        h4("3. Finalizar diseño"),
        p("Al validar los límites, proceda a la fase de monitoreo."), 
        actionButton("pasar_fase2", "Fijar límites y monitorear", class = "btn-primary", icon = icon("arrow-right"))
      ),
      
      # Monitoreo
      conditionalPanel(
        condition = "input.fases_principales == 'Monitoreo'",
        h4("1. Carga de datos"),
        helpText("Estos datos se compararán con los límites fijados en la fase de diseño."), 
        
        fileInput("file2", "Sube nuevos datos (.xlsx / .csv)", multiple = FALSE, accept = c(".csv", ".xlsx")),
        checkboxInput("header2", "El archivo tiene títulos", TRUE),
        uiOutput("opciones_csv_2"),
        uiOutput("panel_agrupacion_2")
      )
    ),
    
    mainPanel(
      # Pestañas principales
      tabsetPanel(id = "fases_principales",
                  
                  # --- Fase 1 ----
                  tabPanel("Diseño y ajuste", 
                           br(),
                           tabsetPanel(id = "tabs_fase1",
                                       tabPanel("Carta de control", 
                                                br(),
                                                # Instrucciones visuales
                                                div(class = "alert alert-info", icon("mouse-pointer"), 
                                                    " Modo interactivo: Haga clic en un punto para eliminarlo/restaurarlo y recalcular los limites."),
                                                # Gráfica interactiva
                                                plotOutput("chart_fase1", height = "500px", click = "click_fase1"),
                                                verbatimTextOutput("resumen_fase1"),
                                                verbatimTextOutput("info_eliminados")
                                       ),
                                       tabPanel("Vista previa de datos", 
                                                br(),
                                                uiOutput("titulo_tabla_1"),
                                                helpText("Nota: Se muestran máximo las primeras 15 filas y 10 columnas."),
                                                tableOutput("tabla_fase1")
                                       )
                           )
                  ),
                  
                  # --- Fase 2 ----
                  tabPanel("Monitoreo", 
                           br(),
                           uiOutput("estado_fase2"), # Mensajes de alerta o exito
                           
                           tabsetPanel(id = "tabs_fase2",
                                       tabPanel("Carta de monitoreo", 
                                                br(),
                                                plotOutput("chart_fase2", height = "500px")
                                       ),
                                       tabPanel("Vista previa de datos", 
                                                br(),
                                                uiOutput("titulo_tabla_2"),
                                                helpText("Nota: Se muestran máximo las primeras 15 filas y 10 columnas."),
                                                tableOutput("tabla_fase2")
                                       )
                           )
                  )
      )
    )
  )
)

# -------------------- Servidor -----------------------
server <- function(input, output, session) {
  
  # Memoria temporal
  valores <- reactiveValues(
    omitidos = numeric(0),       # Lista de puntos eliminados en fase 1
    limites_fijos = NULL,        # Limites fijos para fase 2
    n_fase1 = NULL               # Tamaño de subgrupo original
  )
  
  #  ------------------- Fase 1 - Diseño y ajuste -------------------
  
  # 1. Ventana de ayuda
  observeEvent(input$ayuda_f1, {
    showModal(modalDialog(
      title = "Ayuda en formato de carga de datos",
      h4("Opción A: Subgrupos (Estándar)"),
      p("El archivo tiene varias columnas (Muestra1, Muestra2...). Cada fila es un subgrupo."),
      h4("Opción B: Columna única"),
      p("El archivo tiene una sola columna. El programa le dará la opción de cortar esta columna en subgrupos automáticamente."),
      easyClose = TRUE,
      footer = modalButton("Entendido")
    ))
  })
  
  # 2. Reseteo de memoria al cargar nuevos archivos
  observeEvent(input$file1, { 
    valores$omitidos <- numeric(0)
    valores$limites_fijos <- NULL
    updateTabsetPanel(session, "tabs_fase1", selected = "Vista previa de datos")
  })
  
  observeEvent(input$reset_puntos, { valores$omitidos <- numeric(0) })
  
  observeEvent(input$analyze, {
    updateTabsetPanel(session, "tabs_fase1", selected = "Carta de control")
  })
  
  # 3. Opciones de CSV (Solo aparecen si es .csv)
  output$opciones_csv_1 <- renderUI({
    req(input$file1); if(tools::file_ext(input$file1$name) != "csv") return(NULL)
    wellPanel(radioButtons("sep1", "Separador de columnas:", c("Punto y coma (;)"=";", "Coma (,)"=",", "Tabulación"="\t"), ";"),
              radioButtons("dec1", "Separador decimal:", c("Coma (,)"=",", "Punto (.)"="."), ","))
  })
  
  # 4. Lectura de datos
  datos_crudos_1 <- reactive({
    req(input$file1); archivo <- input$file1; ext <- tools::file_ext(archivo$name); df <- NULL
    tryCatch({
      if (ext == "xlsx") df <- read_excel(archivo$datapath, col_names = input$header1)
      else if (ext == "csv") {
        # separadores por defecto
        mi_sep <- if(is.null(input$sep1)) ";" else input$sep1; mi_dec <- if(is.null(input$dec1)) "," else input$dec1 
        df <- read.csv(archivo$datapath, header = input$header1, sep = mi_sep, dec = mi_dec)
      }
    }, error = function(e) return(NULL))
    
    # Validacion
    if(is.null(df) || ncol(df)==0) return(NULL)
    df_num <- df %>% select(where(is.numeric))
    if(ncol(df_num) == 0) return(NULL) 
    return(df_num)
  })
  
  # 5. Agrupacion cuando es 1 sola columna
  output$panel_agrupacion_1 <- renderUI({
    req(datos_crudos_1()); df <- datos_crudos_1()
    if (ncol(df) == 1) {
      wellPanel(checkboxInput("agrupar1", "Agrupar datos en subgrupos", FALSE),
                conditionalPanel("input.agrupar1 == true", numericInput("n1", "Tamaño (n):", 5, min=2)))
    } else return(NULL)
  })
  
  # 6. Transformacion de datos - matrices
  # Convierte una columna larga en una matriz ancha si el usuario lo pide
  datos_finales_1 <- reactive({
    df <- datos_crudos_1()
    if(is.null(df)) return(NULL)
    
    if (ncol(df) == 1 && isTRUE(input$agrupar1)) {
      vec <- df[[1]]; n <- input$n1; res <- length(vec) %% n
      if (res > 0) vec <- vec[1:(length(vec)-res)]
      mat <- matrix(vec, ncol = n, byrow = TRUE); dft <- as.data.frame(mat); colnames(dft) <- paste0("Muestra_", 1:n)
      return(dft)
    } else return(df)
  })
  
  # 7. Interaccion en la grafica
  observeEvent(input$click_fase1, {
    res <- resultados_fase1() 
    if(!is.null(res)) {
      # nearPoints encuentra el dato más cercano al clic del mouse
      punto <- nearPoints(res$df_grafica, input$click_fase1, xvar="ID", yvar="Media", threshold=10, maxpoints=1)
      if (nrow(punto) > 0) {
        id <- punto$ID
        # Si ya estaba borrado, lo restaura. Si no, lo borra
        if (id %in% valores$omitidos) valores$omitidos <- setdiff(valores$omitidos, id)
        else valores$omitidos <- c(valores$omitidos, id)
      }
    }
  })
  
  # 8. Calculo estadistico
  # Se recalcula cada vez que se presiona "Diseñar" o se eliminan puntos
  resultados_fase1 <- eventReactive(c(input$analyze, valores$omitidos), {
    # validate() muestra mensajes de error en la pantalla
    validate(need(datos_finales_1(), "⚠️ Error de lectura: No se encontraron datos numéricos. Sugerencia: Verifique si el separador de columnas (Punto y coma / Coma) es el correcto en el panel izquierdo."))
    
    if (input$analyze == 0) return(NULL)
    
    df <- datos_finales_1(); n <- ncol(df); m <- nrow(df); B <- input$B_boot
    promedios <- rowMeans(df, na.rm=TRUE)
    validos <- setdiff(1:m, valores$omitidos)
    if(length(validos) < 2) return(NULL)
    
    # Ajuste GAMLSS
    
    datos_est <- as.vector(as.matrix(df[validos, ]))
    datos_est <- as.numeric(datos_est[!is.na(datos_est)]) 
    
    # Verificamos si existe la familia BS6 para evitar crash si no se cargó el archivo externo
    if(!exists("BS6")) {
      mod <- try(stop("Faltan funciones BS6"), silent=TRUE)
    } else {
      mod <- try(gamlss(datos_est ~ 1, family = BS6, control = gamlss.control(trace=FALSE)), silent=TRUE)
    }
    
    if (class(mod)[1] == "try-error") { 
      mu <- mean(datos_est); sigma <- sd(datos_est) 
    } else { 
      mu <- as.numeric(fitted(mod,"mu")[1]); sigma <- as.numeric(fitted(mod,"sigma")[1]) 
      if(is.na(mu)) {mu<-mean(datos_est); sigma<-sd(datos_est)} 
    }
    
    # Simulación Bootstrap
    if(exists("rBS6") && class(mod)[1] != "try-error") {
      boots <- matrix(rBS6(n*B, mu, sigma), B, n)
    } else {
      boots <- matrix(rnorm(n*B, mu, sigma), B, n)
    }
    
    medias_boot <- rowMeans(boots)
    alpha <- 1 - (input$conf_level/100)
    LCI <- quantile(medias_boot, alpha/2); LCS <- quantile(medias_boot, 1-alpha/2)
    
    dfg <- data.frame(ID=1:m, Media=promedios, Estado="Dentro")
    dfg$Estado[dfg$Media > LCS | dfg$Media < LCI] <- "Fuera"
    dfg$Estado[dfg$ID %in% valores$omitidos] <- "Omitido"
    
    list(df_grafica=dfg, limites=c(LCI=LCI, LC=mu, LCS=LCS), params=c(mu, sigma), n=n)
  })
  
  # 9. Salidas graficas y texto fase 1
  output$chart_fase1 <- renderPlot({
    validate(need(datos_finales_1(), "⚠️ Error de lectura: No se encontraron datos numéricos. Sugerencia: Verifique si el separador de columnas (Punto y coma / Coma) es el correcto en el panel izquierdo."),
             need(input$analyze > 0, "Presione 'Diseñar carta' para generar el gráfico."))
    
    req(resultados_fase1()); res <- resultados_fase1(); df <- res$df_grafica; lims <- res$limites
    m_val <- res$n
    if(nrow(df) <= 20) brk <- 1:nrow(df) else brk <- pretty(1:nrow(df), 20)
    
    ggplot(df, aes(ID, Media)) +
      geom_hline(yintercept=lims[c(1,3)], col="red", linetype="dashed") + geom_hline(yintercept=lims[2], col="darkgreen") +
      geom_line(col="gray80") + geom_point(aes(col=Estado, shape=Estado, size=Estado)) +
      scale_color_manual(values=c("Dentro"="black","Fuera"="red","Omitido"="gray90")) +
      scale_shape_manual(values=c("Dentro"=19,"Fuera"=19,"Omitido"=4)) +
      scale_size_manual(values=c("Dentro"=3,"Fuera"=4,"Omitido"=3)) + theme_bw() +
      scale_x_continuous(breaks = brk) +
      scale_y_continuous(n.breaks = 15) + # Más números en el eje Y
      labs(title="Carta de la parametrización Birnbaum Saunders basada en la media", subtitle="Haga clic para omitir puntos.", y="Media muestral", x="Subgrupos") +
      theme(legend.position = "bottom")
  })
  
  output$resumen_fase1 <- renderPrint({ req(resultados_fase1()); cat("--- PARÁMETROS FINALES ---\n"); print(resultados_fase1()$limites) })
  output$info_eliminados <- renderText({ if(length(valores$omitidos)==0) return(""); paste("Subgrupos omitidos:", paste(sort(valores$omitidos), collapse=", ")) })
  
  output$titulo_tabla_1 <- renderUI({
    if(is.null(datos_finales_1())) return(h4("⚠️ Error de datos", style="color:gray;"))
    h4(paste("Dimensiones:", nrow(datos_finales_1()), "x", ncol(datos_finales_1())))
  })
  
  output$tabla_fase1 <- renderTable({ 
    validate(need(datos_finales_1(), "No hay datos válidos para mostrar. Verifique el separador")); 
    head(datos_finales_1(), 15) 
  })
  
  # 10. Transicion fase 2
  observeEvent(input$pasar_fase2, {
    req(resultados_fase1())
    res <- resultados_fase1()
    # Guardamos los limites
    valores$limites_fijos <- res$limites
    valores$n_fase1       <- res$n
    updateTabsetPanel(session, "fases_principales", selected = "Monitoreo") 
    showNotification("✅ Límites fijados. Fase de monitoreo activa.", type="message", duration=4)
  })
  
  
  # ---------------- Fase 2: Monitoreo ------------------------
  
  observeEvent(input$file2, {
    updateTabsetPanel(session, "tabs_fase2", selected = "Vista previa de datos")
  })
  
  output$opciones_csv_2 <- renderUI({
    req(input$file2); if(tools::file_ext(input$file2$name) != "csv") return(NULL)
    wellPanel(radioButtons("sep2", "Separador de columnas:", c("Punto y coma (;)"=";", "Coma (,)"=",", "Tabulación"="\t"), ";"),
              radioButtons("dec2", "Separador decimal:", c("Coma (,)"=",", "Punto (.)"="."), ","))
  })
  
  # Lectura de datos fase 2 (lo mismo que en la fase 1)
  datos_crudos_2 <- reactive({
    req(input$file2); archivo <- input$file2; ext <- tools::file_ext(archivo$name); df <- NULL
    tryCatch({
      if (ext == "xlsx") df <- read_excel(archivo$datapath, col_names = input$header2)
      else if (ext == "csv") {
        mi_sep <- if(is.null(input$sep2)) ";" else input$sep2; mi_dec <- if(is.null(input$dec2)) "," else input$dec2 
        df <- read.csv(archivo$datapath, header = input$header2, sep = mi_sep, dec = mi_dec)
      }
    }, error = function(e) return(NULL))
    
    if(is.null(df) || ncol(df)==0) return(NULL)
    df_num <- df %>% select(where(is.numeric))
    if(ncol(df_num) == 0) return(NULL) # Error si falla separador
    return(df_num)
  })
  
  output$panel_agrupacion_2 <- renderUI({
    req(datos_crudos_2()); df <- datos_crudos_2()
    if (ncol(df) == 1) {
      wellPanel(checkboxInput("agrupar2", "Agrupar datos en subgrupos", FALSE),
                conditionalPanel("input.agrupar2 == true", numericInput("n2", "Tamaño (n):", 5, min=2)))
    } else return(NULL)
  })
  
  datos_finales_2 <- reactive({
    # CORRECCIÓN IMPORTANTE 3: También aquí if(is.null) return(NULL)
    df <- datos_crudos_2()
    if(is.null(df)) return(NULL)
    
    if (ncol(df) == 1 && isTRUE(input$agrupar2)) {
      vec <- df[[1]]; n <- input$n2; res <- length(vec) %% n
      if (res > 0) vec <- vec[1:(length(vec)-res)]
      mat <- matrix(vec, ncol = n, byrow = TRUE); dft <- as.data.frame(mat); colnames(dft) <- paste0("Muestra_", 1:n)
      return(dft)
    } else return(df)
  })
  
  # Mensajes de validacion
  output$estado_fase2 <- renderUI({
    # Prioridad 1: completar fase 1
    if (is.null(valores$limites_fijos)) {
      return(div(class="alert alert-danger", "Primero debe completar 'Diseño y ajuste' y fijar los límites."))
    }
    # Prioridad 2: cargar archivo
    if (is.null(input$file2)) {
      return(div(class="alert alert-warning", "⚠️ Esperando carga de nuevos datos..."))
    }
    # Prioridad 3: error de lectura (separador)
    if (is.null(datos_finales_2())) {
      return(div(class="alert alert-danger", "⚠️ Error de lectura: Verifique el separador de columnas."))
    }
    # Prioridad 4: tamaño de subgrupo correcto
    n_nuevo <- ncol(datos_finales_2()); n_viejo <- valores$n_fase1
    if (n_nuevo != n_viejo) {
      return(div(class="alert alert-warning", paste("⚠️ ADVERTENCIA: Diferente tamaño de subgrupo (n). Diseño:", n_viejo, "vs Monitoreo:", n_nuevo)))
    }
    
    div(class="alert alert-success", "✅ Monitoreando proceso con límites fijos.")
  })
  
  output$chart_fase2 <- renderPlot({
    validate(
      need(valores$limites_fijos, "Complete la Fase 1 primero."),
      need(datos_finales_2(), "⚠️ Error de lectura: No se encontraron datos numéricos. Verifique el separador.")
    )
    
    df <- datos_finales_2(); m <- nrow(df); promedios <- rowMeans(df, na.rm=TRUE)
    lims <- valores$limites_fijos
    
    dfg <- data.frame(ID=1:m, Media=promedios)
    dfg$Color <- ifelse(dfg$Media > lims[3] | dfg$Media < lims[1], "Fuera", "Dentro")
    
    if(m <= 20) brk <- 1:m else brk <- pretty(1:m, 20)
    
    ggplot(dfg, aes(ID, Media)) +
      geom_hline(yintercept=lims[c(1,3)], col="red", linetype="dashed") + geom_hline(yintercept=lims[2], col="darkgreen") +
      geom_line(col="gray50") + geom_point(aes(col=Color), size=3) +
      scale_color_manual(values=c("Dentro"="blue", "Fuera"="red")) +
      scale_x_continuous(breaks = brk) +
      scale_y_continuous(n.breaks = 15) + # Más números en el eje Y
      theme_bw() + theme(legend.position="none") +
      labs(title="Carta de la parametrización Birnbaum Saunders basada en la media", subtitle=paste("Límites fijos (n =", valores$n_fase1, ")"), y="Media muestral", x="Subgrupos")
  })
  
  output$titulo_tabla_2 <- renderUI({
    if(is.null(datos_finales_2())) return(h4("⚠️ Error de datos", style="color:gray;"))
    h4(paste("Dimensiones:", nrow(datos_finales_2()), "x", ncol(datos_finales_2())))
  })
  
  output$tabla_fase2 <- renderTable({ 
    validate(need(datos_finales_2(), "No hay datos válidos.")); 
    head(datos_finales_2(), 15) 
  })
}

shinyApp(ui, server)