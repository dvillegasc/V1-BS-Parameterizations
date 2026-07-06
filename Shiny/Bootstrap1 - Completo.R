library(shiny)
library(ggplot2)
library(dplyr)
library(gamlss)
library(readxl)

# ------------- IMPORTANTE -----------------
# Tener cargado BS6 y dBS6. Esto es lo que yo implemente. Aqui:

# https://github.com/dvillegasc/BS-Parametrizations/blob/main/BS6.R
# https://github.com/dvillegasc/BS-Parametrizations/blob/main/dBS6.R

#---------------------
# --- Interfaz ---

ui <- fluidPage(
  titlePanel("Carta de control Birnbaum Saunders - Basada en la media"),
  
  sidebarLayout(
    sidebarPanel(
      # --- Fase 0: Analisis ---
      conditionalPanel(
        condition = "input.fases_principales == 'Análisis de distribución'",
        h4("1. Carga de datos iniciales"),
        helpText("Sube datos para verificar si la distribución BS es la adecuada."),
        
        fileInput("file0", "Sube tus datos (.xlsx / .csv)", multiple = FALSE, accept = c(".csv", ".xlsx")),
        checkboxInput("header0", "El archivo tiene títulos", TRUE),
        
        uiOutput("opciones_csv_0"),
        uiOutput("panel_agrupacion_0"),
        
        hr(),
        actionButton("check_dist", "Verificar ajuste y comparar", class = "btn-primary btn-lg"),
        
        br(),
        tags$small("Nota: El ajuste puede tardar unos momentos. Por favor espere.", 
                   style = "color: #666; font-style: italic; display:block; margin-top:5px;"),
        
        hr(style = "margin-top: 5px; margin-bottom: 15px;"),
        actionButton("ir_fase1", "Pasar al diseño y ajuste", class = "btn-success", icon = icon("arrow-right"))
      ),
      
      # --- Fase 1: Diseño ---
      conditionalPanel(
        condition = "input.fases_principales == 'Diseño y ajuste'",
        h4("1. Carga de datos"),
        actionLink("ayuda_f1", " Ayuda en formato de carga de datos", icon = icon("question-circle")),
        br(), br(),
        
        fileInput("file1", "Sube tus datos (.xlsx / .csv)", multiple = FALSE, accept = c(".csv", ".xlsx")),
        checkboxInput("header1", "El archivo tiene títulos", TRUE),
        uiOutput("opciones_csv_1"),      
        uiOutput("panel_agrupacion_1"),  
        
        hr(),
        h4("2. Configuración"),
        numericInput("B_boot", "Réplicas bootstrap:", value = 10000, min = 1000, max = 20000),
        sliderInput("conf_level", "Nivel de confianza (%):", min = 90, max = 99.9, value = 99.73),
        
        hr(),
        actionButton("analyze", "Diseñar carta", class = "btn-success btn-lg"),
        
        tags$small("El cálculo Bootstrap puede demorar dependiendo de las réplicas.", 
                   style = "color: #666; font-style: italic; display:block; margin-top:5px;"),
        
        br(), 
        actionButton("reset_puntos", "Restaurar puntos", class = "btn-warning btn-sm", icon = icon("refresh")),
        
        hr(style="border-top: 2px solid #333;"),
        h4("3. Finalizar diseño"),
        p("Al validar los límites, proceda a la fase de monitoreo."), 
        actionButton("pasar_fase2", "Fijar límites y monitorear", class = "btn-primary", icon = icon("arrow-right"))
      ),
      
      # --- Fase 2: Monitoreo ---
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
      tabsetPanel(id = "fases_principales",
                  
                  # Fase 0 pestaña
                  tabPanel("Análisis de distribución",
                           br(),
                           tabsetPanel(id = "tabs_fase0",
                                       tabPanel("Vista previa de datos",
                                                br(),
                                                uiOutput("titulo_tabla_0"),
                                                helpText("Verifique aquí que sus datos se leen correctamente antes de ajustar."),
                                                tableOutput("tabla_fase0")
                                       ),
                                       tabPanel("Ajuste y densidad",
                                                br(),
                                                h4("Comparación visual"),
                                                plotOutput("plot_densidad", height = "450px"),
                                                br(),
                                                h4("Comparación estadística (AIC)"),
                                                helpText("Menor AIC indica mejor ajuste. Se comparan 6 modelos."),
                                                tableOutput("tabla_aic")
                                       ),
                                       tabPanel("Diagnóstico (Worm Plot)",
                                                br(),
                                                helpText("Si los puntos siguen la línea horizontal, el ajuste es bueno."),
                                                plotOutput("plot_worm", height = "500px")
                                       )
                           )
                  ),
                  
                  # Fase 1 pestaña
                  tabPanel("Diseño y ajuste", 
                           br(),
                           tabsetPanel(id = "tabs_fase1",
                                       tabPanel("Carta de control", 
                                                br(),
                                                div(class = "alert alert-info", icon("mouse-pointer"), 
                                                    " Modo interactivo: Haga clic en un punto para eliminarlo/restaurarlo."),
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
                  
                  # Fase 2 pestaña
                  tabPanel("Monitoreo", 
                           br(),
                           uiOutput("estado_fase2"), 
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

##-----------------------------------------------------
# -------------------- Servidor -----------------------

server <- function(input, output, session) {
  
  valores <- reactiveValues(
    omitidos = numeric(0),        
    limites_fijos = NULL,         
    n_fase1 = NULL                
  )
  
  # Mensaje de error
  msg_error_lectura <- paste0(
    
    "\n\n⚠️ ERROR DE LECTURA DE DATOS\n\n",
    "No se encontraron datos numéricos válidos.\n\n",
    "Sugerencia: Verifique si esta incluyendo el titulo en los datos o si el separador de columnas o el decimal es el correcto en el panel izquierdo."
  )
  
  # ===========================================================================
  #                  --- Faso 0: Analisis ---
  
  
  observeEvent(input$check_dist, {
    req(datos_finales_0()) 
    updateTabsetPanel(session, "tabs_fase0", selected = "Ajuste y densidad")
  })
  
  observeEvent(input$ir_fase1, {
    updateTabsetPanel(session, "fases_principales", selected = "Diseño y ajuste")
  })
  
  output$opciones_csv_0 <- renderUI({
    req(input$file0); if(tools::file_ext(input$file0$name) != "csv") return(NULL)
    wellPanel(radioButtons("sep0", "Separador:", c(";"=";", ","=",", "Tab"="\t"), ";"),
              radioButtons("dec0", "Decimal:", c(","=",", "."="."), ","))
  })
  
  datos_crudos_0 <- reactive({
    req(input$file0); archivo <- input$file0; ext <- tools::file_ext(archivo$name); df <- NULL
    tryCatch({
      if (ext == "xlsx") df <- read_excel(archivo$datapath, col_names = input$header0)
      else if (ext == "csv") {
        mi_sep <- if(is.null(input$sep0)) ";" else input$sep0; mi_dec <- if(is.null(input$dec0)) "," else input$dec0 
        df <- read.csv(archivo$datapath, header = input$header0, sep = mi_sep, dec = mi_dec)
      }
    }, error = function(e) return(NULL))
    
    if(is.null(df) || ncol(df)==0) return(NULL)
    df_num <- df %>% select(where(is.numeric))
    if(ncol(df_num) == 0) return(NULL)
    return(df_num)
  })
  
  output$panel_agrupacion_0 <- renderUI({
    req(datos_crudos_0()); df <- datos_crudos_0()
    if (ncol(df) == 1) {
      wellPanel(checkboxInput("agrupar0", "Agrupar datos en subgrupos", FALSE),
                conditionalPanel("input.agrupar0 == true", numericInput("n0", "Tamaño (n):", 5, min=2)))
    } else return(NULL)
  })
  
  datos_finales_0 <- reactive({
    df <- datos_crudos_0()
    if(is.null(df)) return(NULL) 
    if (ncol(df) == 1 && isTRUE(input$agrupar0)) {
      vec <- df[[1]]; n <- input$n0; res <- length(vec) %% n
      if (res > 0) vec <- vec[1:(length(vec)-res)]
      mat <- matrix(vec, ncol = n, byrow = TRUE); dft <- as.data.frame(mat); colnames(dft) <- paste0("Muestra", 1:n)
      return(dft)
    } else return(df)
  })
  
  output$titulo_tabla_0 <- renderUI({
    if(is.null(datos_finales_0())) return(h4("⚠️ Esperando datos válidos o corrección de formato...", style="color:gray;"))
    h4(paste("Datos cargados - Dimensiones:", nrow(datos_finales_0()), "x", ncol(datos_finales_0())))
  })
  
  output$tabla_fase0 <- renderTable({
    validate(need(datos_finales_0(), msg_error_lectura))
    head(datos_finales_0(), 15)
  })
  
  # Calculo AIC
  ajuste_distribucion <- eventReactive(input$check_dist, {
    validate(need(datos_finales_0(), msg_error_lectura))
    
    df <- datos_finales_0()
    datos <- as.vector(as.matrix(df))
    datos <- as.numeric(datos[!is.na(datos)])
    validate(need(length(datos) >= 2, "⚠️ Datos insuficientes. Se necesitan al menos 2 datos."))
    
    m_no  <- try(gamlss(datos ~ 1, family = NO, trace=FALSE), silent=TRUE)
    m_ga  <- try(gamlss(datos ~ 1, family = GA, trace=FALSE), silent=TRUE)
    m_ln  <- try(gamlss(datos ~ 1, family = LOGNO, trace=FALSE), silent=TRUE)
    m_wei <- try(gamlss(datos ~ 1, family = WEI, trace=FALSE), silent=TRUE) 
    m_ig  <- try(gamlss(datos ~ 1, family = IG, trace=FALSE), silent=TRUE)  
    
    if(exists("BS6")) {
      m_bs <- try(gamlss(datos ~ 1, family = BS6, trace=FALSE), silent=TRUE)
    } else {
      m_bs <- try(stop("Falta BS6"), silent=TRUE)
    }
    
    sacar_aic <- function(mod) if(inherits(mod, "try-error")) Inf else AIC(mod)
    
    aics <- data.frame(
      Modelo = c("Birnbaum-Saunders (BS6)", "Weibull (WEI)", "Inv. Gaussiana (IG)", 
                 "Lognormal (LOGNO)", "Gamma (GA)", "Normal (NO)"),
      AIC = c(sacar_aic(m_bs), sacar_aic(m_wei), sacar_aic(m_ig), 
              sacar_aic(m_ln), sacar_aic(m_ga), sacar_aic(m_no))
    )
    aics <- aics[order(aics$AIC), ] 
    
    list(datos = datos, modelos = list(BS=m_bs), tabla = aics)
  })
  
  output$plot_densidad <- renderPlot({
    validate(need(input$check_dist > 0, "Presione 'Verificar ajuste y comparar' para iniciar."))
    req(ajuste_distribucion())
    
    res <- ajuste_distribucion()
    datos <- res$datos
    m_bs <- res$modelos$BS
    
    hist(datos, freq = FALSE, col = "lightblue", border = "white", 
         main = "Histograma y ajuste Birnbaum-Saunders basada en la media", xlab = "Datos", 
         ylim = c(0, max(density(datos)$y)*1.2))
    lines(density(datos), lty=2, col="black", lwd=2) 
    
    if(!inherits(m_bs, "try-error") && exists("dBS6")) {
      mu_hat <- fitted(m_bs, "mu")[1]
      sigma_hat <- fitted(m_bs, "sigma")[1]
      curve(dBS6(x, mu=mu_hat, sigma=sigma_hat), add=TRUE, col="red", lwd=3)
      legend("topright", legend=c("Densidad observada", "Densidad teórica"), col=c("black", "red"), lty=c(2,1), lwd=2)
    }
  })
  
  output$tabla_aic <- renderTable({
    req(ajuste_distribucion())
    ajuste_distribucion()$tabla
  }, digits = 2)
  
  output$plot_worm <- renderPlot({
    validate(need(input$check_dist > 0, "Presione 'Verificar ajuste y comparar'."))
    req(ajuste_distribucion())
    m_bs <- ajuste_distribucion()$modelos$BS
    
    validate(need(!inherits(m_bs, "try-error"), "No se pudo ajustar el modelo (posibles datos negativos)."))
    wp(m_bs, ylim.all = 1) 
    title("Worm Plot: Diagnóstico de residuos (BS6)")
  })
  
  outputOptions(output, "plot_worm", suspendWhenHidden = FALSE)
  outputOptions(output, "plot_densidad", suspendWhenHidden = FALSE)
  
  
  # ===========================================================================
  #                  --- Fase 1: Diseño y ajuste ---
  
  
  observeEvent(input$ayuda_f1, {
    showModal(modalDialog(
      title = "Ayuda en formato de carga de datos",
      h4("Opción A: Subgrupos (Estándar)"),
      p("El archivo tiene varias columnas (Muestra1, Muestra2...)."),
      h4("Opción B: Columna única"),
      p("El archivo tiene una sola columna."),
      easyClose = TRUE, footer = modalButton("Entendido")
    ))
  })
  
  observeEvent(input$file1, { 
    valores$omitidos <- numeric(0)
    valores$limites_fijos <- NULL
    updateTabsetPanel(session, "tabs_fase1", selected = "Vista previa de datos")
  })
  
  observeEvent(input$reset_puntos, { valores$omitidos <- numeric(0) })
  
  observeEvent(input$analyze, {
    updateTabsetPanel(session, "tabs_fase1", selected = "Carta de control")
  })
  
  #Reciclo codigo
  output$opciones_csv_1 <- renderUI({
    req(input$file1); if(tools::file_ext(input$file1$name) != "csv") return(NULL)
    wellPanel(radioButtons("sep1", "Separador de columnas:", c(";"=";", ","=",", "Tab"="\t"), ";"),
              radioButtons("dec1", "Separador decimal:", c(","=",", "."="."), ","))
  })
  
  datos_crudos_1 <- reactive({
    req(input$file1); archivo <- input$file1; ext <- tools::file_ext(archivo$name); df <- NULL
    tryCatch({
      if (ext == "xlsx") df <- read_excel(archivo$datapath, col_names = input$header1)
      else if (ext == "csv") {
        mi_sep <- if(is.null(input$sep1)) ";" else input$sep1; mi_dec <- if(is.null(input$dec1)) "," else input$dec1 
        df <- read.csv(archivo$datapath, header = input$header1, sep = mi_sep, dec = mi_dec)
      }
    }, error = function(e) return(NULL))
    if(is.null(df) || ncol(df)==0) return(NULL)
    df_num <- df %>% select(where(is.numeric))
    if(ncol(df_num) == 0) return(NULL) 
    return(df_num)
  })
  
  output$panel_agrupacion_1 <- renderUI({
    req(datos_crudos_1()); df <- datos_crudos_1()
    if (ncol(df) == 1) {
      wellPanel(checkboxInput("agrupar1", "Agrupar datos en subgrupos", FALSE),
                conditionalPanel("input.agrupar1 == true", numericInput("n1", "Tamaño (n):", 5, min=2)))
    } else return(NULL)
  })
  
  datos_finales_1 <- reactive({
    df <- datos_crudos_1()
    if(is.null(df)) return(NULL)
    if (ncol(df) == 1 && isTRUE(input$agrupar1)) {
      vec <- df[[1]]; n <- input$n1; res <- length(vec) %% n
      if (res > 0) vec <- vec[1:(length(vec)-res)]
      mat <- matrix(vec, ncol = n, byrow = TRUE); dft <- as.data.frame(mat); colnames(dft) <- paste0("Muestra", 1:n)
      return(dft)
    } else return(df)
  })
  
  observeEvent(input$click_fase1, {
    res <- resultados_fase1() 
    if(!is.null(res)) {
      punto <- nearPoints(res$df_grafica, input$click_fase1, xvar="ID", yvar="Media", threshold=10, maxpoints=1)
      if (nrow(punto) > 0) {
        id <- punto$ID
        if (id %in% valores$omitidos) valores$omitidos <- setdiff(valores$omitidos, id)
        else valores$omitidos <- c(valores$omitidos, id)
      }
    }
  })
  
  resultados_fase1 <- eventReactive(c(input$analyze, valores$omitidos), {
    validate(need(datos_finales_1(), msg_error_lectura))
    if (input$analyze == 0) return(NULL)
    
    df <- datos_finales_1(); n <- ncol(df); m <- nrow(df); B <- input$B_boot
    promedios <- rowMeans(df, na.rm=TRUE)
    validos <- setdiff(1:m, valores$omitidos)
    if(length(validos) < 2) return(NULL)
    
    datos_est <- as.vector(as.matrix(df[validos, ]))
    datos_est <- as.numeric(datos_est[!is.na(datos_est)]) 
    
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
  
  output$chart_fase1 <- renderPlot({
    validate(need(datos_finales_1(), msg_error_lectura),
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
      scale_y_continuous(n.breaks = 15) + 
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
    validate(need(datos_finales_1(), msg_error_lectura)); 
    head(datos_finales_1(), 15) 
  })
  
  observeEvent(input$pasar_fase2, {
    req(resultados_fase1())
    res <- resultados_fase1()
    valores$limites_fijos <- res$limites
    valores$n_fase1       <- res$n
    updateTabsetPanel(session, "fases_principales", selected = "Monitoreo") 
    showNotification("✅ Límites fijados. Fase de monitoreo activa.", type="message", duration=4)
  })
  
  
  # ===========================================================================
  #                  --- Fase 2: Monitoreo ---
  
  
  observeEvent(input$file2, {
    updateTabsetPanel(session, "tabs_fase2", selected = "Vista previa de datos")
  })
  
  #Reciclo codigo
  output$opciones_csv_2 <- renderUI({
    req(input$file2); if(tools::file_ext(input$file2$name) != "csv") return(NULL)
    wellPanel(radioButtons("sep2", "Separador de columnas:", c(";"=";", ","=",", "Tab"="\t"), ";"),
              radioButtons("dec2", "Separador decimal:", c(","=",", "."="."), ","))
  })
  
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
    if(ncol(df_num) == 0) return(NULL) 
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
    df <- datos_crudos_2()
    if(is.null(df)) return(NULL)
    
    if (ncol(df) == 1 && isTRUE(input$agrupar2)) {
      vec <- df[[1]]; n <- input$n2; res <- length(vec) %% n
      if (res > 0) vec <- vec[1:(length(vec)-res)]
      mat <- matrix(vec, ncol = n, byrow = TRUE); dft <- as.data.frame(mat); colnames(dft) <- paste0("Muestra", 1:n)
      return(dft)
    } else return(df)
  })
  
  output$estado_fase2 <- renderUI({
    if (is.null(valores$limites_fijos)) {
      return(div(class="alert alert-danger", "Primero debe completar 'Diseño y ajuste' y fijar los límites."))
    }
    if (is.null(input$file2)) {
      return(div(class="alert alert-warning", "⚠️ Esperando carga de nuevos datos..."))
    }
    if (is.null(datos_finales_2())) {
      return(div(class="alert alert-danger", "⚠️ Error de lectura: Verifique el separador de columnas."))
    }
    n_nuevo <- ncol(datos_finales_2()); n_viejo <- valores$n_fase1
    if (n_nuevo != n_viejo) {
      return(div(class="alert alert-warning", paste("⚠️ ADVERTENCIA: Diferente tamaño de subgrupo (n). Diseño:", n_viejo, "vs Monitoreo:", n_nuevo)))
    }
    div(class="alert alert-success", "✅ Monitoreando proceso con límites fijos.")
  })
  
  output$chart_fase2 <- renderPlot({
    validate(
      need(valores$limites_fijos, "Complete la Fase 1 primero."),
      need(datos_finales_2(), msg_error_lectura)
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
      scale_y_continuous(n.breaks = 15) + 
      theme_bw() + theme(legend.position="none") +
      labs(title="Carta de la parametrización Birnbaum Saunders basada en la media", subtitle=paste("Límites fijos (n =", valores$n_fase1, ")"), y="Media muestral", x="Subgrupos")
  })
  
  output$titulo_tabla_2 <- renderUI({
    if(is.null(datos_finales_2())) return(h4("⚠️ Error de datos", style="color:gray;"))
    h4(paste("Dimensiones:", nrow(datos_finales_2()), "x", ncol(datos_finales_2())))
  })
  
  output$tabla_fase2 <- renderTable({ 
    validate(need(datos_finales_2(), msg_error_lectura)); 
    head(datos_finales_2(), 15) 
  })
}

shinyApp(ui, server)