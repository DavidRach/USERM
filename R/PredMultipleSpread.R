#' Launch an interactive Shiny app to explore predicted spread for multiple populations
#'
#' This function launches a Shiny application that allows users to interactively explore the predicted spread
#' of unmixed signals for multiple populations based on a given userm object (\code{Userm}).
#' The app provides multiple tabs for visualizing the unmixing matrix, prediction plots, slope and intercept matrices,
#' and exporting results.
#'
#' @param Userm A userm object created by \code{\link{CreateUserm}}, containing unmixing matrix, detector and fluor information, and intensity data.
#' @param population_ids A character vector specifying one or more column names in \code{Userm$Intensity_mtx} representing the populations to be visualized.
#'
#' @return A Shiny UI object that can be passed to \code{shinyApp()} or used within a Shiny server function.
#'
#' @details
#' The app includes the following features:
#' \itemize{
#'   \item Fluor and detector selection panels
#'   \item Interactive prediction plot with scale and axis controls
#'   \item Visualization of unmixing, pseudoinverse, slope, and intercept matrices
#'   \item Export options for HTML reports and NxN/Nx1 plots
#' }
#'
#' @examples
#' \dontrun{
#' PredMultipleSpread(Userm, population_ids = c("P1", "P2"))
#' }
#' @export

PredMultipleSpread = function(Userm,population_ids){

  missing_cols <- setdiff(population_ids, colnames(Userm$Intensity_mtx))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns from 'population_ids' are missing in Userm$Intensity_mtx:",
               paste(missing_cols, collapse = ", ")))
  }


  ui <- fluidPage(
    titlePanel("USERM"),

    sidebarLayout(
      sidebarPanel(
        tabsetPanel(id = "input_tabs",
                    tabPanel("Fluors",
                             tagList(
                               lapply(Userm$fluors, function(fluor) {
                                 checkboxInput(inputId = fluor, label = fluor, value = TRUE)
                               })
                             )
                    ),
                    tabPanel("Detectors",
                             tagList(
                               lapply(Userm$detectors, function(detector) {
                                 checkboxInput(inputId = detector, label = detector, value = TRUE)
                               })
                             )
                    )
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("unmixing_matrix",
                   uiOutput("unmixing_matrix_html")
          ),
          tabPanel("prediction",
                   selectInput("prediction_plot_mode", "select display mode：", choices = NULL),

                   # X axis controls in one row
                   fluidRow(
                     column(4, selectInput("fluor_selector_x_prediction", "select x axis：", choices = NULL)),
                     column(4, numericInput("x_min_input", "Min", value = 0)),
                     column(4, numericInput("x_max_input", "Max", value = 1000))
                   ),
                   fluidRow(
                     column(4, selectInput("x_scale", "Scale", choices = NULL)),
                     column(8, numericInput("x_cofactor", "Cofactor(only for Arcsinh)", value = 100, min = 1))
                   ),

                   # Y axis controls in one row
                   fluidRow(
                     column(4, selectInput("fluor_selector_y_prediction", "select y axis：", choices = NULL)),
                     column(4, numericInput("y_min_input", "Min", value = 0)),
                     column(4, numericInput("y_max_input", "Max", value = 1000))
                   ),
                   fluidRow(
                     column(4, selectInput("y_scale", "Scale", choices = NULL)),
                     column(8, numericInput("y_cofactor", "Cofactor(only for Arcsinh)", value = 100, min = 1))
                   ),

                   plotOutput("prediction_plot")
          ),

          tabPanel("pinv_matrix",
                   uiOutput("pinv_matrix_html")
          ),
          tabPanel("intercept_matrix",
                   tabsetPanel(
                     tabPanel("raw",
                              selectInput("fluor_selector_intercept_matrix_raw", "select fluor SCC file：", choices = NULL),
                              uiOutput("intercept_matrix_raw_html")),
                     tabPanel("weighted",
                              uiOutput("intercept_matrix_weighted_html"))
                   )
          ),
          tabPanel("slop_matrix",
                   tabsetPanel(
                     tabPanel("raw",
                              selectInput("fluor_selector_slop_matrix_raw", "select fluor SCC file：", choices = NULL),
                              uiOutput("slop_matrix_raw_html"))
                   )
          ),
          tabPanel("Export",
                   selectInput("prediction_plot_mode_report", "select display mode：", choices = c("Pseudo-color","Contour line"), selected = "Pseudo-color"),
                   downloadButton("downloadHTML", "Export HTML report"),
                   downloadButton("downloadNxN", "Export N x N plot"),
                   p("Note: all parameters are from input Userm object."),
                   selectInput("prediction_Nx1_fluor", "select fluor for Nx1 plot：", choices = NULL),
                   downloadButton("downloadNx1", "Export N x 1 plot")


          )

        )
      )

    )
  )


  # Server
  server <- function(input, output, session) {

    # Initialize cache object

    #detector_cache
    detector_cache <- reactiveValues(selected = Userm$detectors)
    observe({
      selected <- Userm$detectors[unlist(lapply(Userm$detectors, function(det) input[[det]]))]
      detector_cache$selected <- selected
    })

    #fluor_cache
    fluor_cache <- reactiveValues(selected = Userm$fluors)
    observe({
      selected <- Userm$fluors[unlist(lapply(Userm$fluors, function(fluor) input[[fluor]]))]
      fluor_cache$selected <- selected
    })

    #intensity_cache
    intensity_cache <- reactiveValues()
    intensity_matrix = Userm$Intensity_mtx[,population_ids,drop = FALSE]
    intensity_cache$matrix = intensity_matrix

    #intercept_matrix_cache
    intercept_matrix_cache <- reactiveValues()

    #observe and calculate intercept_matrix
    observe({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)

      weighted_matrix <- array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)

      for (i in 1:ncol(A)) {
        fluor = colnames(A)[i]
        intercept_matrix = Userm$Res[[fluor]]$interceptMtx
        intercept_matrix = intercept_matrix[detector_cache$selected,detector_cache$selected]
        weighted_matrix[,i] = diag((A_pinv %*% intercept_matrix) %*% t(A_pinv))
      }

      intercept_matrix_cache$matrix <- weighted_matrix
    })


    # unmixing_matrix_table
    output$unmixing_matrix_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      A = Userm$A
      mat = A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      HTML(Userm_html_table(mat = mat, val_min = -1, val_mid = 0, val_max = 1,
                           colormin = "#2166ac", colormid = "#f7f7f7", colormax = "#b2182b"))
    })

    # pinv_matrix_table
    output$pinv_matrix_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      A = Userm$A
      mat = A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      mat_pinv = ginv(mat)
      colnames(mat_pinv) = detector_cache$selected
      rownames(mat_pinv) = fluor_cache$selected
      HTML(Userm_html_table(mat = mat_pinv, val_min = -1, val_mid = 0, val_max = 1,
                           colormin = "#2166ac", colormid = "#f7f7f7", colormax = "#b2182b"))
    })

    # intercept_matrix_table
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "fluor_selector_intercept_matrix_raw",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
    })
    output$intercept_matrix_raw_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      intercept_matrix = Userm$Res[[input$fluor_selector_intercept_matrix_raw]]$interceptMtx
      intercept_matrix = intercept_matrix[detector_cache$selected,detector_cache$selected]
      HTML(Userm_html_table(mat = intercept_matrix, val_min = -10, val_mid = 0, val_max = 10,
                           colormin = "#2166ac", colormid = "#f7f7f7", colormax = "#b2182b"))
    })
    output$intercept_matrix_weighted_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]

      req(intercept_matrix_cache$matrix)
      weighted_matrix = intercept_matrix_cache$matrix

      final_col = matrix(apply(weighted_matrix,1,median),ncol = 1)
      sqrt_col = sqrt(abs(final_col))
      final_matrix = cbind(weighted_matrix, final_col)
      final_matrix = cbind(final_matrix, sqrt_col)
      colnames(final_matrix) = c(colnames(A),"median(Sigma^2)","abs(Sigma)")
      rownames(final_matrix) = colnames(A)

      HTML(Userm_html_table(mat = final_matrix, val_min = -10, val_mid = 0, val_max = 10,
                           colormin = "#2166ac", colormid = "#f7f7f7", colormax = "#b2182b"))
    })

    # slop_matrix_table
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "fluor_selector_slop_matrix_raw",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
    })
    output$slop_matrix_raw_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      slop_matrix = Userm$Res[[input$fluor_selector_slop_matrix_raw]]$slopMtx
      slop_matrix = slop_matrix[detector_cache$selected,detector_cache$selected]
      HTML(Userm_html_table(mat = slop_matrix, val_min = -1, val_mid = 0, val_max = 1,
                           colormin = "#2166ac", colormid = "#f7f7f7", colormax = "#b2182b"))
    })

    #prediction
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "prediction_plot_mode",
                        choices = c("Pseudo-color","Contour line"),
                        selected = "Pseudo-color")
      updateSelectInput(session, "fluor_selector_x_prediction",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
      updateSelectInput(session, "fluor_selector_y_prediction",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
      updateSelectInput(session, "x_scale",
                        choices = c("Linear","Log10","Arcsinh"),
                        selected = "Linear")
      updateSelectInput(session, "y_scale",
                        choices = c("Linear","Log10","Arcsinh"),
                        selected = "Linear")

    })
    output$prediction_plot <- renderPlot({
      req(detector_cache$selected, fluor_cache$selected)
      req(intercept_matrix_cache$matrix)
      intercept_matrix = intercept_matrix_cache$matrix
      intercept_col = matrix(apply(intercept_matrix,1,median),ncol = 1)

      # A = Userm$A
      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)

      # create grid
      x <- create_seq(min = input$x_min_input,
                      max = input$x_max_input, len=100,
                      scale = input$x_scale, cofactor = input$x_cofactor)
      y <- create_seq(min = input$y_min_input,
                      max = input$y_max_input, len=100,
                      scale = input$y_scale, cofactor = input$y_cofactor)
      grid <- expand.grid(x = x, y = y)
      grid$z = 0

      for (popultion in population_ids) {

        grid_tmp = grid

        weighted_matrix = array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)
        intensity_col = array(0, dim = c(ncol(A), 1)) #(channel, 1)
        for (i in 1:ncol(A)) {
          fluor = colnames(A)[i]
          slop_matrix = Userm$Res[[fluor]]$slopMtx
          slop_matrix = slop_matrix[detector_cache$selected,detector_cache$selected]
          weighted_matrix[,i] = diag((A_pinv %*% slop_matrix) %*% t(A_pinv)) * intensity_cache$matrix[fluor,popultion]
          intensity_col[i,1] = intensity_cache$matrix[fluor,popultion]
        }
        spread_col = matrix(apply(weighted_matrix,1,sum),ncol = 1)

        final_col = spread_col + intercept_col

        sqrt_col = sqrt(abs(final_col))
        rownames(sqrt_col) = colnames(A)
        rownames(intensity_col) = colnames(A)
        #Calculate iso band data

        mu_x = intensity_col[input$fluor_selector_x_prediction,1]
        mu_y = intensity_col[input$fluor_selector_y_prediction,1]
        sigma_x = sqrt_col[input$fluor_selector_x_prediction,1]
        sigma_y = sqrt_col[input$fluor_selector_y_prediction,1]

        # calculate intensity for x and y
        fx <- dnorm(grid_tmp$x, mean = mu_x, sd = sigma_x)
        fy <- dnorm(grid_tmp$y, mean = mu_y, sd = sigma_y)

        # Joint probability density (due to independence, it equals the product of individual densities)
        grid_tmp$z <- fx * fy

        grid$z = grid$z + grid_tmp$z
      }


      #visualize iso band data
      p = PredIsobandPlot(grid = grid,
                          x_scale = input$x_scale,
                          y_scale = input$y_scale,
                          x_cofactor = input$x_cofactor,
                          y_cofactor = input$y_cofactor,
                          x_label = input$fluor_selector_x_prediction,
                          y_label = input$fluor_selector_y_prediction,
                          x_min = input$x_min_input,
                          x_max = input$x_max_input,
                          y_min = input$y_min_input,
                          y_max = input$y_max_input,
                          mode = input$prediction_plot_mode,
                          label_population = population_ids,
                          intensity_matrix = intensity_cache$matrix)
      p

    })

    #Export
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "prediction_Nx1_fluor",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
    })
    #downloadHTML
    output$downloadHTML <- downloadHandler(
      filename = function() {
        paste0("UsermMultiple_report-", format(Sys.time(), "%Y-%m-%d_%H_%M_%S"), ".html")
      },
      content = function(file) {
        tempReport <- system.file("report_tmp", "UsermMultiplereport.Rmd", package = "USERM")
        # tempReport <- file.path(report_tmp_dir, "UsermMultiplereport.Rmd")
        file.copy("UsermMultiplereport.Rmd", tempReport, overwrite = TRUE)

        params <- list(Userm = Userm,
                       detector_selected = detector_cache$selected,
                       fluor_selected = fluor_cache$selected,
                       prediction_plot_mode_report = input$prediction_plot_mode_report,
                       population_ids = population_ids
        )

        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          output_format = "html_document",
                          envir = new.env(parent = globalenv()))
      }
    )
    #downloadNxN
    output$downloadNxN <- downloadHandler(
      filename = function() {
        paste0("UsermMultiple_NxN-", format(Sys.time(), "%Y-%m-%d_%H_%M_%S"), ".html")
      },
      content = function(file) {
        tempReport <- system.file("report_tmp", "UsermMultipleNxN.Rmd", package = "USERM")
        # tempReport <- file.path(report_tmp_dir, "UsermMultipleNxN.Rmd")
        file.copy("UsermMultipleNxN.Rmd", tempReport, overwrite = TRUE)

        params <- list(Userm = Userm,
                       detector_selected = detector_cache$selected,
                       fluor_selected = fluor_cache$selected,
                       prediction_plot_mode_report = input$prediction_plot_mode_report,
                       population_ids = population_ids)

        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          output_format = "html_document",
                          envir = new.env(parent = globalenv()))
      }
    )
    #downloadNx1
    output$downloadNx1 <- downloadHandler(
      filename = function() {
        paste0("UsermMultiple_Nx1-", format(Sys.time(), "%Y-%m-%d_%H_%M_%S"), ".html")
      },
      content = function(file) {
        tempReport <- system.file("report_tmp", "UsermMultipleNx1.Rmd", package = "USERM")
        # tempReport <- file.path(report_tmp_dir, "UsermMultipleNx1.Rmd")
        file.copy("UsermMultipleNx1.Rmd", tempReport, overwrite = TRUE)

        params <- list(Userm = Userm,
                       detector_selected = detector_cache$selected,
                       fluor_selected = fluor_cache$selected,
                       prediction_plot_mode_report = input$prediction_plot_mode_report,
                       x_fluor = input$prediction_Nx1_fluor,
                       population_ids = population_ids)

        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          output_format = "html_document",
                          envir = new.env(parent = globalenv()))
      }
    )

  }

  # running App
  shinyApp(ui = ui, server = server)

}
