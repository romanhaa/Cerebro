##----------------------------------------------------------------------------##
## Cerebro
## version 1.0
##
## Author:    Roman Hillje
## Institute: IEO
## Lab:       PGP
## Date:      
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
if (grepl(tolower(Sys.info()['sysname']), pattern='^win')) {
    
  .libPaths(paste0(getwd(), "/R-Portable-Win/library"))

  plot_export_path <- paste0(Sys.getenv("USERPROFILE"), "\\Desktop\\")

} else {

  .libPaths(paste0(getwd(), "/R-Portable-Mac/library"))

  plot_export_path <- "~/Desktop/"

}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
if ( !is.element(el = "BiocManager", set = rownames(installed.packages())) ) {
  install.packages(
    "BiocManager",
    repos = "http://cran.us.r-project.org",
    dependencies = TRUE
  )
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
required_packages_CRAN <- c(
    "DT",
    "formattable",
    "ggplot2",
    "ggtree",
    "Matrix",
    "msigdbr",
    "plotly",
    "RColorBrewer",
    "reshape2",
    "scales",
    "scatterD3",
    "shiny",
    "shinydashboard",
    "shinyWidgets"
  )

for ( package in required_packages_CRAN ) {
  if ( !is.element(el = package, set = rownames(installed.packages())) ) {
    BiocManager::install(package)
  }
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
source("resources/descriptions.txt")

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
library("dplyr")
library("formattable")
library("shiny")
library("shinydashboard")

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
# system("type R")
# print(sessionInfo())
# Sys.setenv("R_LIBS_USER" = "")
# print(Sys.getenv())

##----------------------------------------------------------------------------##
## Check if data already exists.
##----------------------------------------------------------------------------##
mode <- ifelse(file.exists("resources/data.rds"), "boxed", "open")

##----------------------------------------------------------------------------##
## Allow upload of files up to 400 MB.
##----------------------------------------------------------------------------##
options(shiny.maxRequestSize = 400*1024^2) 

##----------------------------------------------------------------------------##
## Color management.
##----------------------------------------------------------------------------##
# Dutch palette from flatuicolors.com
colors_dutch <- c("#FFC312","#C4E538","#12CBC4","#FDA7DF","#ED4C67",
                  "#F79F1F","#A3CB38","#1289A7","#D980FA","#B53471",
                  "#EE5A24","#009432","#0652DD","#9980FA","#833471",
                  "#EA2027","#006266","#1B1464","#5758BB","#6F1E51")

# Spanish palette from flatuicolors.com
colors_spanish <- c("#40407a","#706fd3","#f7f1e3","#34ace0","#33d9b2",
                    "#2c2c54","#474787","#aaa69d","#227093","#218c74",
                    "#ff5252","#ff793f","#d1ccc0","#ffb142","#ffda79",
                    "#b33939","#cd6133","#84817a","#cc8e35","#ccae62")

colors <- c(colors_dutch, colors_spanish)

cell_cycle_colorset <- setNames(
    c("#45aaf2","#f1c40f","#e74c3c", "#7f8c8d"),
    c("G1",     "S",      "G2M",     "-")
  )

##----------------------------------------------------------------------------##
## Server.
##----------------------------------------------------------------------------##
server <- function(input, output, session) {

  ##--------------------------------------------------------------------------##
  ## Sidebar menu.
  ##--------------------------------------------------------------------------##
  output$sidebar_menu <- renderMenu({
    if ( mode == "open" ) {
      sidebarMenu(id = "sidebar",
        menuItem(
          "Load data", tabName = "loadData",
          icon = icon("spinner"), selected = TRUE
        ),
        menuItem(
          "Overview", tabName = "overview",
          icon = icon("binoculars")
        ),
        menuItem(
          "Samples", tabName = "samples",
          icon = icon("star")
        ),
        menuItem(
          "Clusters", tabName = "clusters",
          icon = icon("braille")
        ),
        menuItem(
          "Top expressed genes", tabName = "topExpressedGenes",
          icon = icon("bullhorn")
        ),
        menuItem(
          "Marker genes", tabName = "markerGenes",
          icon = icon("magnet")
        ),
        menuItem(
          "Enriched pathways", tabName = "enrichedPathways",
          icon = icon("sitemap")
        ),
        menuItem(
          "Gene expression", tabName = "geneExpression",
          icon = icon("signal")
        ),
        menuItem(
          "Gene set expression", tabName = "geneSetExpression",
          icon = icon("list")
        ),
        menuItem(
          "Gene ID conversion", tabName = "geneIdConversion",
          icon = icon("barcode")
        ),
        menuItem(
          "Sample info", tabName = "info",
          icon = icon("info")
        ),
        menuItem(
          "About", tabName = "about",
          icon = icon("at")
        )
      )
    } else {
      sidebarMenu(id = "sidebar",
        menuItem(
          "Load data", tabName = "loadData", icon = icon("spinner")
        ),
        menuItem(
          "Overview", tabName = "overview", icon = icon("binoculars"),
          selected = TRUE
        ),
        menuItem(
          "Samples", tabName = "samples", icon = icon("star")
        ),
        menuItem(
          "Clusters", tabName = "clusters", icon = icon("braille")
        ),
        menuItem(
          "Top expressed genes", tabName = "topExpressedGenes",
          icon = icon("bullhorn")
        ),
        menuItem(
          "Marker genes", tabName = "markerGenes", icon = icon("magnet")
        ),
        menuItem(
          "Enriched pathways", tabName = "enrichedPathways",
          icon = icon("sitemap")
        ),
        menuItem(
          "Gene expression", tabName = "geneExpression", icon = icon("signal")
        ),
        menuItem(
          "Gene set expression", tabName = "geneSetExpression",
          icon = icon("list")
        ),
        menuItem(
          "Gene ID conversion", tabName = "geneIdConversion",
          icon = icon("barcode")
        ),
        menuItem(
          "Sample info", tabName = "info", icon = icon("info")
        ),
        menuItem(
          "About", tabName = "about", icon = icon("at")
        )
      )
    }
  })

  ##--------------------------------------------------------------------------##
  ## Sample data.
  ##--------------------------------------------------------------------------##
  sample_data <- reactive({

    if ( mode == "boxed" ) {
      sample_data <- readRDS("resources/data.rds")
    } else if ( is.null(input$RDS_file) || is.na(input$RDS_file) ) {
      sample_data <- readRDS("resources/example.rds")
    } else {
      req(input$RDS_file)
      sample_data <- readRDS(input$RDS_file$datapath)
    }

    sample_data
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Overview.
  ##--------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##

  output$overview_UI <- renderUI({
    tagList(
      selectInput("overview_projection_to_display", label = "Projection",
        choices = names(sample_data()$projections)),
      selectInput("overview_projection_cell_color", label = "Color cells by",
        choices = names(sample_data()$cells)[! names(sample_data()$cells) %in% c("cell_barcode")]),
      checkboxInput("overview_projection_ellipses",
        label = "Confidence ellipses",
        value = FALSE),
      shinyWidgets::pickerInput("overview_samples_to_display",
        label = "Samples to display",
        choices = sample_data()$samples$overview$sample,
        selected = sample_data()$samples$overview$sample,
        options = list("actions-box" = TRUE),
        multiple = TRUE),
      shinyWidgets::pickerInput("overview_clusters_to_display",
        label = "Clusters to display",
        choices = sample_data()$clusters$overview$cluster,
        selected = sample_data()$clusters$overview$cluster,
        options = list("actions-box" = TRUE),
        multiple = TRUE),
      selectInput("overview_projection_cell_size_variable",
        label = "Change point size by",
        choices = c("None", "nUMI", "nGene"),
        selected = "None"),
      sliderInput("overview_projection_cell_size", label = "Point size",
        min = 0, max = 50, value = 15, step = 1),
      sliderInput("overview_projection_cell_opacity", label = "Point opacity",
        min = 0, max = 1, value = 1, step = 0.05)
    )
  })

  ##--------------------------------------------------------------------------##

  output$overview_scales <- renderUI({
    projection_to_display <- if ( is.null(input$overview_projection_to_display) || is.na(input$overview_projection_to_display) ) names(sample_data()$projections)[1] else input$overview_projection_to_display
    range_x_min <- round(min(sample_data()$projections[[ projection_to_display ]][,1])*1.1)
    range_x_max <- round(max(sample_data()$projections[[ projection_to_display ]][,1])*1.1)
    range_y_min <- round(min(sample_data()$projections[[ projection_to_display ]][,2])*1.1)
    range_y_max <- round(max(sample_data()$projections[[ projection_to_display ]][,2])*1.1)
    tagList(
      sliderInput("overview_projection_scale_x_manual_range",
        label = "X axis",
        min = range_x_min,
        max = range_x_max,
        value = c(range_x_min, range_x_max)
      ),
      sliderInput("overview_projection_scale_y_manual_range",
        label = "Y axis",
        min = range_y_min,
        max = range_y_max,
        value = c(range_y_min, range_y_max)
      )
    )
  })

  ##--------------------------------------------------------------------------##

  output$overview.projection <- scatterD3::renderScatterD3({

    # don't do anything before these inputs are selected
    req(input$overview_projection_to_display)
    req(input$overview_samples_to_display)
    req(input$overview_clusters_to_display)
    req(input$overview_projection_scale_x_manual_range)

    # define which projection should be plotted
    if ( is.null(input$overview_projection_to_display) || is.na(input$overview_projection_to_display) ) {
      projection_to_display <- names(sample_data()$projections)[1]
    } else {
      projection_to_display <- input$overview_projection_to_display
    }
    
    # define which samples should be plotted
    if ( is.null(input$overview_samples_to_display) || is.na(input$overview_samples_to_display) ) {
      samples_to_display <- sample_data()$samples$overview$sample
    } else {
      samples_to_display <- input$overview_samples_to_display
    }

    # define which clusters should be plotted
    if ( is.null(input$overview_clusters_to_display) || is.na(input$overview_clusters_to_display) ) {
      clusters_to_display <- sample_data()$clusters$overview$cluster
    } else {
      clusters_to_display <- input$overview_clusters_to_display
    }

    # define which cells should be plotted
    cells_to_display <- which(
        grepl(
          sample_data()$cells$sample,
          pattern = paste0("^", samples_to_display, "$", collapse="|")
        ) & 
        grepl(
          sample_data()$cells$cluster,
          pattern = paste0("^", clusters_to_display, "$", collapse="|")
        )
      )

    # extract cells to plot
    to_plot <- cbind(
        sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
        sample_data()$cells[ cells_to_display , ]
      )
    to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]

    # define variable used to color cells by
    col_var <- to_plot[ , input$overview_projection_cell_color ]

    # define colors
    if ( is.null(input$overview_projection_cell_color) || is.na(input$overview_projection_cell_color) ) {
      colors <- NULL
    } else if ( input$overview_projection_cell_color == "sample" ) {
      colors <- sample_data()$samples$colors
    } else if ( input$overview_projection_cell_color == "cluster" ) {
      colors <- sample_data()$clusters$colors
    } else if ( input$overview_projection_cell_color %in% c("cell_cycle_Regev","cell_cycle_Cyclone") ) {
      colors <- cell_cycle_colorset
    } else if ( is.factor(to_plot[,input$overview_projection_cell_color]) ) {
      colors <- setNames(colors[1:length(levels(to_plot[,input$overview_projection_cell_color]))], levels(to_plot[,input$overview_projection_cell_color]))
    } else if ( is.character(to_plot[,input$overview_projection_cell_color]) ) {
      colors <- colors
    } else {
      colors <- NULL
    }

    # define variable used for cell size
    size_var <- if ( input$overview_projection_cell_size_variable == "None" ) NULL else to_plot[ , input$overview_projection_cell_size_variable ]

    # plot
    scatterD3::scatterD3(
      x = to_plot[ , 1 ],
      y = to_plot[ , 2 ],
      xlab = colnames(to_plot)[ 1 ],
      ylab = colnames(to_plot)[ 2 ],
      xlim = c(
          input$overview_projection_scale_x_manual_range[1],
          input$overview_projection_scale_x_manual_range[2]
        ),
      ylim = c(
          input$overview_projection_scale_y_manual_range[1],
          input$overview_projection_scale_y_manual_range[2]
        ),
      point_size = input$overview_projection_cell_size,
      col_var = col_var,
      col_lab = input$overview_projection_cell_color,
      colors = colors,
      ellipses = input$overview_projection_ellipses,
      size_var = size_var,
      point_opacity = input$overview_projection_cell_opacity,
      transitions = FALSE,
      menu = FALSE,
      tooltip_text  = paste0(
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br/>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br/>",
        "<b>nUMI</b>: ", to_plot[ , "nUMI" ], "<br/>",
        "<b>nGene</b>: ", to_plot[ , "nGene" ], "<br/>",
        "<b>Expr. MT</b>: ", format(to_plot[ , "percent_mt" ]*100, digits=1), "%<br/>",
        "<b>Expr. ribo</b>: ", format(to_plot[ , "percent_ribo" ]*100, digits=1), "%<br/>"))
  })

  observeEvent(input$overview_projection_info, {
    showModal(
      modalDialog(
        overview_projection_info_text,
        title = overview_projection_info_title, easyClose = TRUE, footer = NULL
      )
    )
  })

  observeEvent(input$overview_projection_export, {
    library("ggplot2")

    projection_to_display <- input$overview_projection_to_display
    samples_to_display <- input$overview_samples_to_display
    clusters_to_display <- input$overview_clusters_to_display
    cells_to_display <- which(
        grepl(
          sample_data()$cells$sample,
          pattern = paste0("^", samples_to_display, "$", collapse = "|")
        ) &
        grepl(
          sample_data()$cells$cluster,
          pattern = paste0("^", clusters_to_display, "$", collapse = "|")
        )
      )
    to_plot <- cbind(
        sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
        sample_data()$cells[ cells_to_display , ]
      )
    to_plot <- to_plot[ sample(1:nrow(to_plot)) , ]

    xlim <- c(
        input$overview_projection_scale_x_manual_range[1],
        input$overview_projection_scale_x_manual_range[2]
      )
    ylim <- c(
        input$overview_projection_scale_y_manual_range[1],
        input$overview_projection_scale_y_manual_range[2]
      )

    if ( is.factor(to_plot[,input$overview_projection_cell_color]) | is.character(to_plot[,input$overview_projection_cell_color]) ) {
      if ( input$overview_projection_cell_color == "sample" ) {
        cols <- sample_data()$samples$colors
      } else if ( input$overview_projection_cell_color == "cluster" ) {
        cols <- sample_data()$clusters$colors
      } else if ( input$overview_projection_cell_color %in% c("cell_cycle_Regev","cell_cycle_Cyclone") ) {
        cols <- cell_cycle_colorset
      } else if ( is.factor(to_plot[,input$overview_projection_cell_color]) ) {
        cols <- setNames(colors[1:length(levels(to_plot[,input$overview_projection_cell_color]))], levels(to_plot[,input$overview_projection_cell_color]))
      } else {
        cols <- colors
      }
      p <- ggplot(
          to_plot,
          aes_q(
            x = as.name(colnames(to_plot)[1]),
            y = as.name(colnames(to_plot)[2]),
            colour = as.name(input$overview_projection_cell_color)
          )
        ) +
        geom_point() +
        scale_colour_manual(values = cols) +
        lims(x = xlim, y = ylim) +
        theme_bw()
    } else {
      p <- ggplot(
          to_plot,
          aes_q(
            x = as.name(colnames(to_plot)[1]),
            y = as.name(colnames(to_plot)[2]),
            colour = as.name(input$overview_projection_cell_color)
          )
        ) +
        geom_point() +
        viridis::scale_colour_viridis(
          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
        ) +
        lims(x = xlim, y = ylim) +
        theme_bw()
    }

    out_filename <- paste0(
        plot_export_path, "Cerebro_",
        gsub(
          sample_data()$experiment$experiment_name,
          pattern = " ", replacement = "_"
        ),
        "_overview_", input$overview_projection_to_display, "_by_",
        gsub(
          input$overview_projection_cell_color,
          pattern = "\\.", replacement = "_"
        ),
        ".pdf"
      )

    pdf(NULL)
    ggsave(out_filename, p, height = 8, width = 11)

    if ( file.exists(out_filename) ) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Success!",
        text = paste0("Plot saved successfully as: ", out_filename),
        type = "success"
      )
    } else {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Error!",
        text = "Sorry, it seems something went wrong...",
        type = "error"
      )
    }
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Samples.
  ##--------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##

  output$samples_overview <- DT::renderDataTable({
    sample_data()$samples$overview %>%
    rename(
      Sample = sample,
      Color = color,
      "# of cells" = number_of_cells
    ) %>%
    mutate(
      "Mean number of UMI" = round(mean_nUMI, digits = 1),
      "Mean number of expressed genes" = round(mean_nGene, digits = 1)
    ) %>%
    select(-c(path, data_type, mean_nUMI, mean_nGene)) %>%    
    DT::datatable(
      filter = "none",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      class = "cell-border stripe",
      options = list(
        scrollX = TRUE,
        sDom = '<"top">lrt<"bottom">ip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15
      )
    )
  })

  observeEvent(input$samples_overview_info, {
    showModal(
      modalDialog(samples_overview_info_text,
        title = samples_overview_info_title, easyClose = TRUE, footer = NULL
      )
    )
  })

  ##--------------------------------------------------------------------------##

  output$samples_by_cluster_table <- DT::renderDataTable({
    sample_data()$samples$by_cluster %>%
    rename(
      Sample = sample,
      "# of cells" = total_cell_count
    ) %>%
    DT::datatable(
      filter = "none",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      class = "cell-border stripe",
      options = list(
        scrollX = TRUE,
        sDom = '<"top">lrt<"bottom">ip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15
      )
    )
  })

  observeEvent(input$samples_by_cluster_table_info, {
    showModal(
      modalDialog(samples_by_cluster_table_info_text,
        title = samples_by_cluster_table_info_title, easyClose = TRUE,
        footer = NULL
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of samples by cluster
  output$samples_by_cluster_UI <- renderUI({
    if ( nrow(sample_data()$clusters$overview) > 1 ) {
      plotly::plotlyOutput("samples_by_cluster_plot")
    } else {
      textOutput("samples_by_cluster_text")
    }
  })

  # bar plot of samples by cluster
  output$samples_by_cluster_plot <- plotly::renderPlotly({
    sample_data()$samples$by_cluster %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "sample") %>%
    rename(cluster = variable, cells = value) %>%
    left_join(
      .,
      sample_data()$samples$by_cluster[ , c("sample", "total_cell_count") ],
      by = "sample"
    ) %>%
    mutate(pct = cells / total_cell_count) %>%
    plotly::plot_ly(
      x = ~sample,
      y = ~pct*100,
      type = "bar",
      color = ~cluster,
      colors = sample_data()$clusters$colors,
      text = ~pct*100,
      hoverinfo = "name+y"
    ) %>%
    plotly::layout(
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage (%)", hoverformat = ".2f"),
      barmode = "stack",
      hovermode = "compare"
    ) 
  })

  # alternative text output for bar plot of samples by cluster
  output$samples_by_cluster_text <- renderText({
      "Only 1 cluster in this data set."
    })

  observeEvent(input$samples_by_cluster_info, {
    showModal(
      modalDialog(
        title = "Samples by cluster", easyClose = TRUE, footer = NULL,
        p("Percentage bar plot representation of the table shown above. Allows to see which clusters contribute most strongly to each sample. Clusters can be removed from the plot by clicking on them in the legend.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of transcripts per sample
  output$samples_box_nUMI <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~sample,
      y = ~nUMI,
      type = "box",
      color = ~sample,
      colors = sample_data()$samples$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Number of UMIs", type = "log", hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$samples_box_nUMI_info, {
    showModal(
      modalDialog(
        title = "Number of transcripts", easyClose = TRUE, footer = NULL,
        p("Box plot of the number of transcripts (UMIs) found in each sample.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of expressed genes per sample
  output$samples_box_nGene <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~sample,
      y = ~nGene,
      type = "box",
      color = ~sample,
      colors = sample_data()$samples$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Number of expressed genes", type = "log",
        hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$samples_box_nGene_info, {
    showModal(
      modalDialog(
        title = "Number of expressed genes", easyClose = TRUE, footer = NULL,
        p("Box plot of the number of expressed genes found in each sample.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of mitochondrial gene expression per sample
  output$samples_box_percent_mt <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~sample,
      y = ~percent_mt,
      type = "box",
      color = ~sample,
      colors = sample_data()$samples$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage of mitochondrial gene expression",
        range = c(0,1), hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$samples_box_percent_mt_info, {
    showModal(
      modalDialog(
        title = "Mitochondrial gene expression", easyClose = TRUE,
        footer = NULL,
        p("Box plot of the percentage of mitochondrial gene expression found in each sample. This reflects the contribution of mitochondrial transcripts to the entire transcriptome in each cell. A list of all genes considered to be mitochondrial can be found in the 'Sample info' tab on the left.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of ribosomal gene expression per sample
  output$samples_box_percent_ribo <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~sample,
      y = ~percent_ribo,
      type = "box",
      color = ~sample,
      colors = sample_data()$samples$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage of ribosomal gene expression",
        range = c(0,1), hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$samples_box_percent_ribo_info, {
    showModal(
      modalDialog(
        title = "Ribosomal gene expression", easyClose = TRUE, footer = NULL,
        p("Box plot of the percentage of ribosomal gene expression found in each sample. This reflects the contribution of ribosomal transcripts to the entire transcriptome in each cell. A list of all genes considered to be ribosomal can be found in the 'Sample info' tab on the left.")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  
  # UI element for bar plot of samples by cell cycle (Regev)
  output$samples_by_cell_cycle_Regev_UI <- renderUI({
    if ( !is.null(sample_data()$samples$by_cell_cycle_Regev) ) {
      plotly::plotlyOutput("samples_by_cell_cycle_Regev_plot")
    } else {
      textOutput("samples_by_cell_cycle_Regev_text")
    }
  })

  # bar plot of samples by cell cycle (Regev)  
  output$samples_by_cell_cycle_Regev_plot <- plotly::renderPlotly({
    sample_data()$samples$by_cell_cycle_Regev %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "sample") %>%
    rename(phase = variable, cells = value) %>%
    mutate(
      phase = factor(phase, levels = c("G1", "S", "G2M")),
    ) %>%
    left_join(
      .,
      sample_data()$samples$by_cell_cycle_Regev[ , c("sample", "total_cell_count") ],
      by = "sample"
    ) %>%
    mutate(pct = cells / total_cell_count) %>%
    plotly::plot_ly(
      x = ~sample,
      y = ~pct*100,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      text = ~pct*100,
      hoverinfo = "name+y"
    ) %>%
    plotly::layout(
      xaxis = list(title =""),
      yaxis = list(title = "Percentage (%)", hoverformat = ".2f"),
      barmode = "stack",
      hovermode = "compare"
    ) 
  })

  # alternative text for bar plot of samples by cell cycle (Regev)  
  output$samples_by_cell_cycle_Regev_text <- renderText({
      "Data not available."
    })

  observeEvent(input$samples_by_cell_cycle_Regev_info, {
    showModal(
      modalDialog(
        title = "Cell cycle analysis (Regev)", easyClose = TRUE, footer = NULL,
        p("Cell cycle distribution by sample using the method embedded in the Seurat framework. For each cell, it calculates scores for both G2M and S phase based on lists of genes (see 'Sample info' tab on the left) and assigns the cell cycle phase on the basis of these scores.")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  
  # UI element for bar plot of samples by cell cycle (Cyclone)
  output$samples_by_cell_cycle_Cyclone_UI <- renderUI({
    if ( !is.null(sample_data()$samples$by_cell_cycle_Cyclone) ) {
      plotly::plotlyOutput("samples_by_cell_cycle_Cyclone_plot")
    } else {
      textOutput("samples_by_cell_cycle_Cyclone_text")
    }
  })

  # bar plot of samples by cell cycle (Cyclone)
  output$samples_by_cell_cycle_Cyclone_plot <- plotly::renderPlotly({
    sample_data()$samples$by_cell_cycle_Cyclone %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "sample") %>%
    rename(phase = variable, cells = value) %>%
    mutate(
      phase = factor(phase, levels = c("G1", "S", "G2M", "-")),
    ) %>%
    left_join(
      .,
      sample_data()$samples$by_cell_cycle_Cyclone[ , c("sample", "total_cell_count") ],
      by = "sample"
    ) %>%
    mutate(pct = cells / total_cell_count) %>%
    plotly::plot_ly(
      x = ~sample,
      y = ~pct*100,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      text = ~pct*100,
      hoverinfo = "name+y"
    ) %>%
    plotly::layout(
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage (%)", hoverformat = ".2f"),
      barmode = "stack",
      hovermode = "compare"
    ) 
  })

  # alternative text for bar plot of samples by cell cycle (Cyclone)
  output$samples_by_cell_cycle_Cyclone_text <- renderText({
      "Data not available."
    })

  observeEvent(input$samples_by_cell_cycle_Cyclone_info, {
    showModal(
      modalDialog(
        title = "Cell cycle analysis (Cyclone)", easyClose = TRUE,
        footer = NULL,
        p("Cell cycle distribution by sample using the machine learning-based Cyclone method published by Scialdone et al (2015). It assigns the cell cycle phase based on scores calculated using relative expression of lists of gene pairs. In contrast to the Seurat/Regev method, scores are calculated for G1 and G2M phase. Cells with a low score for both are assigned S phase. Inability to predict the cell cycle phase for a given cell with this method is most likely a result of very few expressed genes in the respective cell.")
      )
    )
  }) 


  ##--------------------------------------------------------------------------##
  ## Panel: Clusters.
  ##--------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$samples$count
  ## - sample_data()$clusters$overview$cluster
  ## - sample_data()$clusters$tree
  ## - sample_data()$cells$cluster
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent_mt
  ## - sample_data()$cells$percent_ribo
  ## - sample_data()$cells$cell_cycle_Regev (optional)
  ## - sample_data()$cells$cell_cycle_Cyclone (optional)
  ##--------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##

  # UI element for cluster tree
  output$clusters_tree_UI <- renderUI({
    if ( !is.null(sample_data()$clusters$tree) ) {
      plotOutput("clusters_tree_plot")
    } else {
      textOutput("clusters_tree_text")
    }
  })

  output$clusters_tree_plot <- renderPlot({
    library("ggtree")
    tree <- sample_data()$clusters$tree
    tree$tip.label <- paste0("Cluster ", tree$tip.label)
    colors_tree <- colors[1:length(tree$tip.label)]
    ggplot(tree, aes(x, y)) + 
      ggplot2::scale_y_reverse() +
      xlim(0, max(tree$edge.length * 1.1)) +
      geom_tree() +
      theme_tree() +
      geom_tiplab(size = 5, hjust = -0.2) +
      geom_tippoint(color = colors_tree, shape = 16, size = 6)
  })

  output$clusters_tree_text <- renderText({ "Data not available." })

  observeEvent(input$clusters_tree_info, {
    showModal(
      modalDialog(
        title = "Cluster tree", easyClose = TRUE, footer = NULL,
        p("Cluster tree reflecting the similarity of clusters based on their expression profiles. Instead of using the expression values, the correlation is calculated using the user-specified number of principal components (see 'Sample info' tab on the left).")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # cell counts by cluster and sample
  output$clusters_by_sample <- DT::renderDataTable({
    sample_data()$clusters$by_sample %>%
    rename(
      Cluster = cluster,
      "# of cells" = total_cell_count
    ) %>%
    DT::datatable(
      filter = "none",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      class = "cell-border stripe",
      options = list(
        scrollX = TRUE,
        sDom = '<"top">lrt<"bottom">ip',
        lengthMenu = c(20, 30, 50, 100),
        pageLength = 20
      )
    )
  })

  observeEvent(input$clusters_by_sample_info, {
    showModal(
      modalDialog(
        title = "Overview of clusters", easyClose = TRUE, footer = NULL,
        p("Table of clusters (by row) stratified by sample (columns).")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of clusters by sample
  output$clusters_by_sample_UI <- renderUI({
    if ( nrow(sample_data()$samples$overview) > 1 ) {
      plotly::plotlyOutput("clusters_by_sample_plot")
    } else {
      textOutput("clusters_by_sample_text")
    }
  })

  # bar plot of clusters by sample
  output$clusters_by_sample_plot <- plotly::renderPlotly({
    sample_data()$clusters$by_sample %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "cluster") %>%
    rename(sample = variable, cells = value) %>%
    left_join(
      .,
      sample_data()$clusters$by_sample[ , c("cluster", "total_cell_count") ],
      by = "cluster"
    ) %>%
    mutate(pct = cells / total_cell_count) %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~pct*100,
      type = "bar",
      color = ~sample,
      colors = sample_data()$samples$colors,
      text = ~pct*100,
      hoverinfo = "name+y"
    ) %>%
    plotly::layout(
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage (%)", hoverformat = ".2f"),
      barmode = "stack",
      hovermode = "compare"
    )
  })

  # alternative text for bar plot of clusters by sample
  output$clusters_by_sample_text <- renderText({
      "Only 1 sample in this data set."
    })

  observeEvent(input$clusters_by_sample_info, {
    showModal(
      modalDialog(
        title = "Clusters by samples", easyClose = TRUE, footer = NULL,
        p("Percentage bar plot representation of the table shown above. Allows to see which samples contribute most strongly to each cluster. Samples can be removed from the plot by clicking on them in the legend.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of transcripts per cluster
  output$clusters_box_nUMI <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~cluster,
      y = ~nUMI,
      type = "box",
      color = ~cluster,
      colors = sample_data()$clusters$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""), 
      yaxis = list(title = "Number of UMIs", type = "log", hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$clusters_box_nUMI_info, {
    showModal(
      modalDialog(
        title = "Number of transcripts", easyClose = TRUE, footer = NULL,
        p("Box plot of the number of transcripts (UMIs) found in each cluster.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of expressed genes per cluster
  output$clusters_box_nGene <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~cluster, 
      y = ~nGene,
      type = "box",
      color = ~cluster,
      colors = sample_data()$clusters$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title =""),
      yaxis = list(title = "Number of expressed genes", type = "log",
          hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$clusters_box_nGene_info, {
    showModal(
      modalDialog(
        title = "Number of expressed genes", easyClose = TRUE, footer = NULL,
        p("Box plot of the number of expressed genes found in each cluster.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of mitochondrial gene expression per cluster
  output$clusters_box_percent_mt <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~cluster,
      y = ~percent_mt,
      type = "box",
      color = ~cluster,
      colors = sample_data()$clusters$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage of mitochondrial gene expression",
        range = c(0, 1), hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$clusters_box_percent_mt_info, {
    showModal(
      modalDialog(
        title = "Mitochondrial gene expression", easyClose = TRUE,
        footer = NULL,
        p("Box plot of the percentage of mitochondrial gene expression found in each cluster. This reflects the contribution of mitochondrial transcripts to the entire transcriptome in each cell. A list of all genes considered to be mitochondrial can be found in the 'Sample info' tab on the left.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of ribosomal gene expression per cluster
  output$clusters_box_percent_ribo <- plotly::renderPlotly({
    plotly::plot_ly(
      sample_data()$cells,
      x = ~cluster,
      y = ~percent_ribo,
      type = "box",
      color = ~cluster,
      colors = sample_data()$clusters$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage of ribosomal gene expression",
        range = c(0, 1), hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$clusters_box_percent_ribo_info, {
    showModal(
      modalDialog(
        title = "Ribosomal gene expression", easyClose = TRUE, footer = NULL,
        p("Box plot of the percentage of ribosomal gene expression found in each cluster. This reflects the contribution of ribosomal transcripts to the entire transcriptome in each cell. A list of all genes considered to be ribosomal can be found in the 'Sample info' tab on the left.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of clusters by cell cycle (Regev)
  output$clusters_by_cell_cycle_Regev_UI <- renderUI({
    if ( !is.null(sample_data()$clusters$by_cell_cycle_Regev) ) {
      plotly::plotlyOutput("clusters_by_cell_cycle_Regev_plot")
    } else {
      textOutput("clusters_by_cell_cycle_Regev_text")
    }
  })

  # bar plot of clusters by cell cycle (Regev)
  output$clusters_by_cell_cycle_Regev_plot <- plotly::renderPlotly({
    sample_data()$clusters$by_cell_cycle_Regev %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "cluster") %>%
    rename(phase = variable, cells = value) %>%
    mutate(
      phase = factor(phase, levels = c("G1", "S", "G2M")),
    ) %>%
    left_join(
      .,
      sample_data()$clusters$by_cell_cycle_Regev[ , c("cluster", "total_cell_count") ],
      by = "cluster"
    ) %>%
    mutate(pct = cells / total_cell_count) %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~pct*100,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      text = ~pct*100,
      hoverinfo = "name+y"
    ) %>%
    plotly::layout(
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage (%)", hoverformat = ".2f"),
      barmode = "stack",
      hovermode = "compare"
    ) 
  })

  # alternative text for bar plot of clusters by cell cycle (Regev)
  output$clusters_by_cell_cycle_Regev_text <- renderText({
      "Data not available."
    })

  observeEvent(input$clusters_by_cell_cycle_Regev_info, {
    showModal(
      modalDialog(
        title = "Cell cycle analysis (Regev)", easyClose = TRUE, footer = NULL,
        p("Cell cycle distribution by cluster using the method embedded in the Seurat framework. For each cell, it calculates scores for both G2M and S phase based on lists of genes (see 'Sample info' tab on the left) and assigns the cell cycle phase on the basis of these scores.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of clusters by cell cycle (Cyclone)
  output$clusters_by_cell_cycle_Cyclone_UI <- renderUI(
    if ( !is.null(sample_data()$cells$cell_cycle_Cyclone) ) {
      plotly::plotlyOutput("clusters_by_cell_cycle_Cyclone_plot")
    } else {
      textOutput("clusters_by_cell_cycle_Cyclone_text")
    }
  )

  # bar plot of clusters by cell cycle (Cyclone)
  output$clusters_by_cell_cycle_Cyclone_plot <- plotly::renderPlotly({
    sample_data()$clusters$by_cell_cycle_Cyclone %>%
    select(-total_cell_count) %>%
    reshape2::melt(id.vars = "cluster") %>%
    rename(phase = variable, cells = value) %>%
    mutate(
      phase = factor(phase, levels = c("G1", "S", "G2M", "-")),
    ) %>%
    left_join(
      .,
      sample_data()$clusters$by_cell_cycle_Cyclone[ , c("cluster", "total_cell_count") ],
      by = "cluster"
    ) %>%
    mutate(pct = cells / total_cell_count) %>%
    plotly::plot_ly( 
      x = ~cluster,
      y = ~pct*100,
      type = "bar",
      color = ~phase,
      colors = cell_cycle_colorset,
      text = ~pct*100,
      hoverinfo = "name+y"
    ) %>%
    plotly::layout(
      xaxis = list(title = ""),
      yaxis = list(title = "Percentage (%)", hoverformat = ".2f"), 
      barmode = "stack",
      hovermode = "compare"
    )
  })

  # alternative text for bar plot of clusters by cell cycle (Cyclone)
  output$clusters_by_cell_cycle_Cyclone_text <- renderText({
      "Data not available."
    })

  observeEvent(input$clusters_by_cell_cycle_Cyclone_info, {
    showModal(
      modalDialog(
        title = "Cell cycle analysis (Cyclone)", easyClose = TRUE,
        footer = NULL,
        p("Cell cycle distribution by cluster using the machine learning-based Cyclone method published by Scialdone et al (2015). It assigns the cell cycle phase based on scores calculated using relative expression of lists of gene pairs. In contrast to the Seurat/Regev method, scores are calculated for G1 and G2M phase. Cells with a low score for both are assigned S phase.")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Top expressed genes.
  ##--------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$samples$overview$sample
  ## - sample_data()$most_expressed_genes$by_sample
  ## - sample_data()$clusters$overview$cluster
  ## - sample_data()$most_expressed_genes$by_cluster
  ##--------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##

  # by sample
  output$top_expressed_genes_by_sample_UI <- renderUI({
    if ( !is.null(sample_data()$most_expressed_genes$by_sample) ) {
      fluidRow(
        column(12,
          selectInput("top_expressed_genes_by_sample_input", label = NULL,
            choices = sample_data()$samples$overview$sample),
          DT::dataTableOutput("top_expressed_genes_by_sample_table_present")
        )
      )
    } else {
      textOutput("top_expressed_genes_by_sample_table_missing")
    }
  })

  output$top_expressed_genes_by_sample_table_present <- DT::renderDataTable(server = FALSE, {
    req(input$top_expressed_genes_by_sample_input)
    sample_data()$most_expressed_genes$by_sample %>%
    filter(sample == input$top_expressed_genes_by_sample_input) %>%
    mutate(pct = formattable::percent(round(pct/100, digits = 4))) %>%
    rename(
      Sample = sample,
      Gene = gene,
      "% of total expression" = pct
    ) %>%
    formattable::formattable(
      list(
        "Sample" = formattable::color_tile(
            colors[ which(sample_data()$samples$overview$sample == input$top_expressed_genes_by_sample_input) ],
            colors[ which(sample_data()$samples$overview$sample == input$top_expressed_genes_by_sample_input) ]
          ),
        "% of total expression" = formattable::color_bar("pink")
      )
    ) %>%
    formattable::as.datatable(
      filter = "none",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        dom = "Bfrtip",
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons = list(
          "colvis", 
           list(
             extend = "collection",
             text = "Download",
             buttons = list(
              list(
                extend = "csv",
                filename = "top_expressed_genes_per_sample",
                title = "Top expressed genes per sample"
              ),
              list(
                extend = "excel",
                filename = "top_expressed_genes_per_sample",
                title = "Top expressed genes per sample"
              ),
              list(
                extend = "pdf",
                filename = "top_expressed_genes_per_sample",
                title = "Top expressed genes per sample"
              )
            )
          )
        )
      )
    ) %>%
    DT::formatStyle("% of total expression", textAlign = "right")
  })

  output$top_expressed_genes_by_sample_table_missing <- renderText({
      "Data not available."
    })

  observeEvent(input$top_expressed_genes_by_sample_info, {
    showModal(
      modalDialog(
        title = "Top expressed genes per sample", easyClose = TRUE,
        footer = NULL,
        p("Table of top 100 most expressed genes in each sample. For example, if gene XY contributes with 5% to the total expression, that means 5% of all transcripts found in all cells of this sample come from that respective gene. These lists can help to identify/verify the dominant cell types.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # by cluster
  output$top_expressed_genes_by_cluster_UI <- renderUI({
    if ( !is.null(sample_data()$most_expressed_genes$by_cluster) ) {
      fluidRow(
        column(12,
          selectInput("top_expressed_genes_by_cluster_input", label = NULL,
          choices = sample_data()$clusters$overview$cluster),
          DT::dataTableOutput("top_expressed_genes_by_cluster_table_present")
        )
      )
    } else {
      textOutput("top_expressed_genes_by_cluster_table_missing")
    }
  })

  output$top_expressed_genes_by_cluster_table_present <- DT::renderDataTable(server = FALSE, {
    req(input$top_expressed_genes_by_cluster_input)
    sample_data()$most_expressed_genes$by_cluster %>%
    filter(cluster == input$top_expressed_genes_by_cluster_input) %>%
    mutate(pct = formattable::percent(round(pct/100, digits = 4))) %>%
    rename(
      Cluster = cluster,
      Gene = gene,
      "% of total expression" = pct
    ) %>%
    formattable::formattable(
      list(
        "Cluster" = formattable::color_tile(
            colors[ which(sample_data()$clusters$overview$cluster == input$top_expressed_genes_by_cluster_input) ],
            colors[ which(sample_data()$clusters$overview$cluster == input$top_expressed_genes_by_cluster_input) ]
          ),
        "% of total expression" = formattable::color_bar("pink")
      )
    ) %>%
    formattable::as.datatable(
      filter = "none",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        dom = "Bfrtip",
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons = list(
          "colvis", 
           list(
             extend = "collection",
             text = "Download",
             buttons = list(
              list(
                extend = "csv",
                filename = "top_expressed_genes_per_cluster",
                title = "Top expressed genes per cluster"
              ),
              list(
                extend = "excel",
                filename = "top_expressed_genes_per_cluster",
                title = "Top expressed genes per cluster"
              ),
              list(
                extend = "pdf",
                filename = "top_expressed_genes_per_cluster",
                title = "Top expressed genes per cluster"
              )
            )
          )
        )
      )
    ) %>%
    DT::formatStyle("% of total expression", textAlign = "right")
  })

  output$top_expressed_genes_by_cluster_table_missing <- renderText({
      "Data not available."
    })

  observeEvent(input$top_expressed_genes_by_cluster_info, {
    showModal(
      modalDialog(
        title = "Top expressed genes per cluster", easyClose = TRUE,
        footer = NULL,
        p("Table of top 100 most expressed genes in each cluster. For example, if gene XY contributes with 5% to the total expression, that means 5% of all transcripts found in all cells of this cluster come from that respective gene. These lists can help to identify/verify the dominant cell types.")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Marker genes.
  ##--------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$samples$overview$sample
  ## - sample_data()$samples$count
  ## - sample_data()$marker_genes$by_sample (optional)
  ## - sample_data()$clusters$overview$cluster
  ## - sample_data()$clusters$count
  ## - sample_data()$marker_genes$by_cluster (optional)
  ##--------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##

  # by sample
  output$marker_genes_by_sample_UI <- renderUI({
    if ( nrow(sample_data()$samples$overview) > 1 & !is.null(sample_data()$marker_genes$by_sample) ) {
      fluidRow(
        column(12,
          selectInput(
            "marker_genes_by_sample_input", label = NULL,
            choices = sample_data()$samples$overview$sample
          ),
          DT::dataTableOutput("marker_genes_by_sample_table_present")
        )
      )
    } else {
      textOutput("marker_genes_by_sample_table_missing")
    }
  })

  output$marker_genes_by_sample_table_present <- DT::renderDataTable(server = FALSE, {
    req(input$marker_genes_by_sample_input)
    if ( "on_cell_surface" %in% colnames(sample_data()$marker_genes$by_sample) ) {
      table <- sample_data()$marker_genes$by_sample[ which(sample_data()$marker_genes$by_sample$sample == input$marker_genes_by_sample_input) , c(2,4,5,6,7,8) ]
      colnames(table) <- c("Gene", "avg. logFC", "% cells in this sample", "% cells in other samples", "adj. p-value", "present on cell surface")
      table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
      table$"% cells in this sample" <- formattable::percent(table$"% cells in this sample")
      table$"% cells in other samples" <- formattable::percent(table$"% cells in other samples")
      table <- table %>%
      mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
      formattable::formattable(list(
        "avg. logFC" = formattable::color_tile("white", "orange"),
        "% cells in this sample" = formattable::color_bar("pink"),
        "% cells in other samples" = formattable::color_bar("pink"),
        "present on cell surface" = formattable::formatter("span", style=x~style(color=ifelse(x, "green", "red")))
      ))
    } else if ( tolower(sample_data()$experiment$organism) %in% c("hg","mm") ) {
      if ( !exists("genes_surface") ) {
        genes_surface <- read.table(paste0("resources/genes_surface_", tolower(sample_data()$experiment$organism), ".txt"), sep="\t", header=FALSE, stringsAsFactors=FALSE)[,1]
      }
      table <- sample_data()$marker_genes$by_sample[ which(sample_data()$marker_genes$by_sample$sample == input$marker_genes_by_sample_input) , c(2,4,5,6,7) ]
      table$surface <- table$gene %in% genes_surface
      colnames(table) <- c("Gene", "avg. logFC", "% cells in this sample", "% cells in other samples", "adj. p-value", "present on cell surface")
      table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
      table$"% cells in this sample" <- formattable::percent(table$"% cells in this sample")
      table$"% cells in other samples" <- formattable::percent(table$"% cells in other samples")
      table <- table %>%
      mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
      formattable::formattable(list(
        "avg. logFC" = formattable::color_tile("white", "orange"),
        "% cells in this sample" = formattable::color_bar("pink"),
        "% cells in other samples" = formattable::color_bar("pink"),
        "present on cell surface" = formattable::formatter("span", style=x~style(color=ifelse(x, "green", "red")))
      ))
    } else {
      table <- sample_data()$marker_genes$by_sample[ which(sample_data()$marker_genes$by_sample$sample == input$marker_genes_by_sample_input) , c(2,4,5,6,7) ]
      colnames(table) <- c("Gene", "avg. logFC", "% cells in this sample", "% cells in other samples", "adj. p-value")
      table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
      table$"% cells in this sample" <- formattable::percent(table$"% cells in this sample")
      table$"% cells in other samples" <- formattable::percent(table$"% cells in other samples")
      table <- table %>%
      mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
      formattable::formattable(list(
        "avg. logFC" = formattable::color_tile("white", "orange"),
        "% cells in this sample" = formattable::color_bar("pink"),
        "% cells in other samples" = formattable::color_bar("pink")
      ))
    }
    formattable::as.datatable(
      table,
      filter = "top",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        dom = "Bfrtip",
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons = list(
          "colvis", 
          list(
            extend = "collection",
            text = "Download",
            buttons = list(
              list(
                extend = "csv",
                filename = "marker_genes_by_sample",
                title = "Marker genes by sample"
              ),
              list(
                extend = "excel",
                filename = "marker_genes_by_sample",
                title = "Marker genes by sample"
              ),
              list(
                extend = "pdf",
                filename = "marker_genes_by_sample",
                title = "Marker genes by sample"
              )
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(
      columns = c("avg. logFC", "% cells in this sample", "% cells in other samples", "adj. p-value"),
      textAlign = "right"
    )
  })

  output$marker_genes_by_sample_table_missing <- renderText({
      "Only 1 sample in this data set or data not available."
    })

  observeEvent(input$marker_genes_by_sample_info, {
    showModal(
      modalDialog(
        title = "Marker genes per sample", easyClose = TRUE, footer = NULL,
        p("Shown here are the marker genes identified for each sample - resembling bulk RNA-seq. These genes should help to identify the cell type in this sample or find new markers to purify it. In this analysis, each sample is compared to all other samples combined. Only genes with a positive average log-fold change of at least 0.25 are reported - meaning only over-expressed genes are shown. Also, marker genes must be expressed in at least 70% of the cells of the respective sample. Statistical analysis is performed using a classical t-test as it has been shown to be very accurate in single cell RNA-seq. Finally, if data is available, the last column reports for each gene if it is associated with gene ontology term GO:0009986 which is an indicator that the respective gene is present on the cell surface (which could make it more interesting to purify a given population).")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # by cluster
  output$marker_genes_by_cluster_UI <- renderUI({
    if ( nrow(sample_data()$clusters$overview) > 1 & !is.null(sample_data()$marker_genes$by_cluster) ) {
      fluidRow(
        column(12,
          selectInput(
            "marker_genes_by_cluster_input", label = NULL,
            choices = sample_data()$clusters$overview$cluster
          ),
          DT::dataTableOutput("marker_genes_by_cluster_table_present")
        )
      )
    } else {
      textOutput("marker_genes_by_cluster_table_missing")
    }
  })

  output$marker_genes_by_cluster_table_present <- DT::renderDataTable(server = FALSE, {
    req(input$marker_genes_by_cluster_input)
    if ("on_cell_surface" %in% colnames(sample_data()$marker_genes$by_cluster)) {
      table <- sample_data()$marker_genes$by_cluster[ which(sample_data()$marker_genes$by_cluster$cluster == input$marker_genes_by_cluster_input) , c(2,4,5,6,7,8) ]
      colnames(table) <- c("Gene", "avg. logFC", "% cells in this cluster", "% cells in other clusters", "adj. p-value", "present on cell surface")
      table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
      table$"% cells in this cluster" <- formattable::percent(table$"% cells in this cluster")
      table$"% cells in other clusters" <- formattable::percent(table$"% cells in other clusters")
      table <- table %>% mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
      formattable::formattable(list(
        "avg. logFC" = formattable::color_tile("white", "orange"),
        "% cells in this cluster" = formattable::color_bar("pink"),
        "% cells in other clusters" = formattable::color_bar("pink"),
        "present on cell surface" = formattable::formatter("span", style=x~style(color=ifelse(x, "green", "red")))
      ))
    } else if ( tolower(sample_data()$experiment$organism) %in% c("hg","mm") ) {
      if ( !exists("genes_surface") ) {
        genes_surface <- read.table(paste0("resources/genes_surface_", tolower(sample_data()$experiment$organism), ".txt"), sep="\t", header=FALSE, stringsAsFactors=FALSE)[,1]
      }
      table <- sample_data()$marker_genes$by_cluster[ which(sample_data()$marker_genes$by_cluster$cluster == input$marker_genes_by_cluster_input) , c(2,4,5,6,7) ]
      table$surface <- table$gene %in% genes_surface
      colnames(table) <- c("Gene", "avg. logFC", "% cells in this cluster", "% cells in other clusters", "adj. p-value", "present on cell surface")
      table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
      table$"% cells in this cluster" <- formattable::percent(table$"% cells in this cluster")
      table$"% cells in other clusters" <- formattable::percent(table$"% cells in other clusters")
      table <- table %>%
      mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
      formattable::formattable(list(
        "avg. logFC" = formattable::color_tile("white", "orange"),
        "% cells in this cluster" = formattable::color_bar("pink"),
        "% cells in other clusters" = formattable::color_bar("pink"),
        "present on cell surface" = formattable::formatter("span", style=x~style(color=ifelse(x, "green", "red")))
      ))
    } else {
      table <- sample_data()$marker_genes$by_cluster[ which(sample_data()$marker_genes$by_cluster$cluster == input$marker_genes_by_cluster_input) , c(2,4,5,6,7) ]
      colnames(table) <- c("Gene", "avg. logFC", "% cells in this cluster", "% cells in other clusters", "adj. p-value")
      table$"avg. logFC" <- round(table$"avg. logFC", digits=3)
      table$"% cells in this cluster" <- formattable::percent(table$"% cells in this cluster")
      table$"% cells in other clusters" <- formattable::percent(table$"% cells in other clusters")
      table <- table %>% mutate("adj. p-value"=formatC(table$"adj. p-value", format="e", digits=3)) %>%
      formattable::formattable(list(
        "avg. logFC" = formattable::color_tile("white", "orange"),
        "% cells in this cluster" = formattable::color_bar("pink"),
        "% cells in other clusters" = formattable::color_bar("pink")
      ))
    }
    formattable::as.datatable(
      table,
      filter = "top",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        dom = "Bfrtip",
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons = list(
          "colvis", 
          list(
            extend = "collection",
            text = "Download",
            buttons = list(
              list(
                extend = "csv",
                filename = "marker_genes_by_cluster",
                title = "Marker genes by cluster"
              ),
              list(
                extend = "excel",
                filename = "marker_genes_by_cluster",
                title = "Marker genes by cluster"
              ),
              list(
                extend = "pdf",
                filename = "marker_genes_by_cluster",
                title = "Marker genes by cluster"
              )
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(
      columns = c("avg. logFC", "% cells in this cluster", "% cells in other clusters", "adj. p-value"),
      textAlign = "right"
    )
  })

  output$marker_genes_by_cluster_table_missing <- renderText({
      "Only 1 cluster in this data set or data not available."
    })

  observeEvent(input$marker_genes_by_cluster_info, {
    showModal(
      modalDialog(
        title = "Marker genes per cluster", easyClose = TRUE, footer = NULL,
        p("Shown here are the marker genes identified for each cluster. These genes should help to identify the cell type in this cluster or find new markers to purify it. In this analysis, each cluster is compared to all other clusters combined. Only genes with a positive average log-fold change of at least 0.25 are reported - meaning only over-expressed genes are shown. Also, marker genes must be expressed in at least 70% of the cells of the respective cluster. Statistical analysis is performed using a classical t-test as it has been shown to be very accurate in single cell RNA-seq. Finally, if data is available, the last column reports for each gene if it is associated with gene ontology term GO:0009986 which is an indicator that the respective gene is present on the cell surface (which could make it more interesting to purify a given population).")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Enriched pathways
  ##--------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$parameters$enrichr_dbs
  ## - sample_data()$samples$count
  ## - sample_data()$samples$overview$sample
  ## - sample_data()$clusters$count
  ## - sample_data()$clusters$overview$cluster
  ## - sample_data()$marker_genes$by_sample_annotation (optional)
  ## - sample_data()$marker_genes$by_cluster_annotation (optional)
  ##--------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##

  # by sample
  output$enriched_pathways_by_sample_UI <- renderUI({
    if ( nrow(sample_data()$samples$overview) > 1 & !is.null(sample_data()$marker_genes$by_sample_annotation) ) {    
      tagList(
        fluidRow(
          column(4,
            selectInput("enriched_pathways_select_sample", label = NULL,
              choices = sample_data()$samples$overview$sample)
          ),
          column(8,
            selectInput("enriched_pathways_select_db_for_sample", label = NULL,
              choices = sample_data()$parameters$enrichr_dbs)
          )
        ),
        fluidRow(
          column(12,
            DT::dataTableOutput("enriched_pathways_by_sample_table_present")
          )
        )
      )
    } else {
      textOutput("enriched_pathways_by_sample_table_missing")
    }
  })

  output$enriched_pathways_by_sample_table_present <- DT::renderDataTable(server = FALSE, {
    req(input$enriched_pathways_select_sample)
    req(input$enriched_pathways_select_db_for_sample)
    sample_data()$marker_genes$by_sample_annotation[[ input$enriched_pathways_select_sample ]][[ input$enriched_pathways_select_db_for_sample ]] %>%
    select(c(1,2,3,4,8,9)) %>%
    mutate(
      P.value = formatC(P.value, format = "e", digits = 3),
      Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 3),
      Combined.Score = formatC(Combined.Score, format = "f", digits = 2)
    ) %>%
    rename(
      "p-value" = P.value,
      "adj. p-value" = Adjusted.P.value,
      "combined score" = Combined.Score,
    ) %>%
    formattable::formattable(
      list("combined score" = formattable::color_bar("pink"))
    ) %>%
    formattable::as.datatable(
      filter = "top",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        columnDefs = list(list(visible = FALSE, targets = c(2,5))),
        scrollX = TRUE,
        dom = "Bfrtip",
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons = list(
          "colvis", 
          list(
            extend = "collection",
            text = "Download",
            buttons = list(
              list(
                extend = "csv",
                filename = "enriched_pathways_by_sample",
                title = "Enriched pathways by sample"
              ),
              list(
                extend = "excel",
                filename = "enriched_pathways_by_sample",
                title = "Enriched pathways by sample"
              ),
              list(
                extend = "pdf",
                filename = "enriched_pathways_by_sample",
                title = "Enriched pathways by sample"
              )
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(columns = c("combined score"), textAlign = "right")
  })

  output$enriched_pathways_by_sample_table_missing <- renderText({
      "Only 1 sample in this data set or data not available."
    })

  observeEvent(input$enriched_pathways_by_sample_info, {
    showModal(
      modalDialog(
        title = "Enriched pathways by sample", easyClose = TRUE, footer = NULL,
        p("Using all marker genes identified for a respective sample, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the 'Column visibility' button. For the details on the combined score is calculated, please refer to the Enrichr website and publication: http://amp.pharm.mssm.edu/Enrichr/.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # by cluster
  output$enriched_pathways_by_cluster_UI <- renderUI({
    if ( nrow(sample_data()$clusters$overview) > 1 & !is.null(sample_data()$marker_genes$by_cluster_annotation) ) {    
      tagList(
        fluidRow(
          column(4,
            selectInput("enriched_pathways_select_cluster", label = NULL,
              choices = sample_data()$clusters$overview$cluster)
          ),
          column(8,
            selectInput("enriched_pathways_select_db_for_cluster", label = NULL,
              choices = sample_data()$parameters$enrichr_dbs)
          )
        ),
        fluidRow(
          column(12,
            DT::dataTableOutput("enriched_pathways_by_cluster_table_present")
          )
        )
      )
    } else {
      textOutput("enriched_pathways_by_cluster_table_missing")
    }
  })

  output$enriched_pathways_by_cluster_table_present <- DT::renderDataTable(server = FALSE, {
    req(input$enriched_pathways_select_cluster)
    req(input$enriched_pathways_select_db_for_cluster)
    sample_data()$marker_genes$by_cluster_annotation[[ input$enriched_pathways_select_cluster ]][[ input$enriched_pathways_select_db_for_cluster ]][ , c(1,2,3,4,8,9) ] %>%
    mutate(
      P.value = formatC(P.value, format = "e", digits = 3),
      Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 3),
      Combined.Score = formatC(Combined.Score, format = "f", digits = 2)
    ) %>%
    rename(
      "p-value" = P.value,
      "adj. p-value" = Adjusted.P.value,
      "combined score" = Combined.Score,
    ) %>%
    formattable::formattable(
      list("combined score" = formattable::color_bar("pink"))
    ) %>%
    formattable::as.datatable(
      filter = "top",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      extensions = c("Buttons"),
      class = "cell-border stripe",
      options = list(
        columnDefs = list(list(visible = FALSE, targets = c(2,5))),
        scrollX = TRUE,
        dom = "Bfrtip", 
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons = list(
          "colvis",
          list(
            extend = "collection",
            text = "Download",
            buttons = list(
              list(
                extend = "csv",
                filename = "enriched_pathways_by_cluster",
                title = "Enriched pathways by cluster"
              ),
              list(
                extend = "excel",
                filename = "enriched_pathways_by_cluster",
                title = "Enriched pathways by cluster"
              ),
              list(
                extend = "pdf",
                filename = "enriched_pathways_by_cluster",
                title = "Enriched pathways by cluster"
              )
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(columns = c("combined score"), textAlign = "right")
  })

  output$enriched_pathways_by_cluster_table_missing <- renderText({
      "Only 1 cluster in this data set or data not available."
    })

  observeEvent(input$enriched_pathways_by_cluster_info, {
    showModal(
      modalDialog(
        title = "Enriched pathways by cluster", easyClose = TRUE, footer = NULL,
        p("Using all marker genes identified for a respective cluster, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the 'Column visibility' button. For the details on the combined score is calculated, please refer to the Enrichr website and publication: http://amp.pharm.mssm.edu/Enrichr/")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Gene expression
  ##--------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$projections
  ## - sample_data()$samples$overview$sample
  ## - sample_data()$clusters$overview$cluster
  ## - sample_data()$expression
  ## - sample_data()$cells$sample
  ## - sample_data()$cells$cluster
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent_mt
  ## - sample_data()$cells$percent_ribo
  ##--------------------------------------------------------------------------##
 
  # reactive data
  genesToPlot <- reactive({
    genesToPlot <- list()
    if ( is.null(input$expression_genes_input) ) {
      genesToPlot$genes_to_display <- ""
    } else {
      genesToPlot$genes_to_display <- input$expression_genes_input %>%
        strsplit(",| |;|\n") %>%
        unlist() %>%
        gsub(pattern = " ", replacement = "", fixed = TRUE) %>%
        unique()
    }
    genesToPlot$genes_to_display_here <- rownames(sample_data()$expression)[ match(tolower(genesToPlot$genes_to_display), tolower(rownames(sample_data()$expression))) ]
    genesToPlot$genes_to_display_present <- genesToPlot$genes_to_display_here[ which(!is.na(genesToPlot$genes_to_display_here)) ]
    genesToPlot$genes_to_display_missing <- genesToPlot$genes_to_display[ which(is.na(genesToPlot$genes_to_display_here)) ]
    genesToPlot
  })

  # select genes to be displayed
  output$expression_genes_displayed <- renderText({
    paste0(
      "<b>Showing expression for ",
      length(genesToPlot()$genes_to_display_present), " gene(s):</b><br>",
      paste0(genesToPlot()$genes_to_display_present, collapse = ", "),
      "<br><b>",
      length(genesToPlot()$genes_to_display_missing),
      " gene(s) are not in data set: </b><br>",
      paste0(genesToPlot()$genes_to_display_missing, collapse = ", ")
    )
  })

  ##--------------------------------------------------------------------------##
 # UI
  output$expression_UI <- renderUI({
    tagList(
      selectInput("expression_projection_to_display", label = "Projection:",
        choices = names(sample_data()$projections)
      ),
      textAreaInput("expression_genes_input", label = "Gene(s):",
        value = "",
        placeholder = "Insert genes here."
      ),
      shinyWidgets::pickerInput(
        "expression_samples_to_display",
        label = "Samples to display:",
        choices = sample_data()$samples$overview$sample,
        selected = sample_data()$samples$overview$sample,
        options = list("actions-box" = TRUE),
        multiple = TRUE
      ),
      shinyWidgets::pickerInput(
        "expression_clusters_to_display",
        label = "Clusters to display:",
        choices = sample_data()$clusters$overview$cluster,
        selected = sample_data()$clusters$overview$cluster,
        options = list("actions-box" = TRUE),
        multiple = TRUE
      ),
      selectInput("expression_plotting_order", label = "Plotting order:",
        choices = c("Random", "Highest expression on top")
      ),
      sliderInput("expression_projection_dot.size", label = "Point size:",
        min = 0, max = 50, value = 25, step = 1
      ),
      sliderInput("expression_projection_opacity", label = "Point opacity:",
        min = 0, max = 1, value = 1, step = 0.05
      )
    )
  })

  output$expression_scales <- renderUI({
    projection_to_display <- if ( is.null(input$expression_projection_to_display) || is.na(input$expression_projection_to_display) ) names(sample_data()$projections)[1] else input$expression_projection_to_display
    range_x_min <- round(
        min(sample_data()$projections[[ projection_to_display ]][,1]) * 1.1)
    range_x_max <- round(
        max(sample_data()$projections[[ projection_to_display ]][,1]) * 1.1)
    range_y_min <- round(
        min(sample_data()$projections[[ projection_to_display ]][,2]) * 1.1)
    range_y_max <- round(
        max(sample_data()$projections[[ projection_to_display ]][,2]) * 1.1)
    tagList(
      sliderInput(
        "expression_projection_scale_x_manual_range",
        label = "X axis",
        min = range_x_min, max = range_x_max,
        value = c(range_x_min, range_x_max)
      ),
      sliderInput(
        "expression_projection_scale_y_manual_range",
        label = "Y axis",
        min = range_y_min, max = range_y_max,
        value = c(range_y_min, range_y_max)
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # projection with scatterD3
  output$expression_projection <- scatterD3::renderScatterD3({
    req(input$expression_projection_to_display)
    req(input$expression_samples_to_display)
    req(input$expression_clusters_to_display)
    req(input$expression_plotting_order)
    projection_to_display <- input$expression_projection_to_display
    samples_to_display <- input$expression_samples_to_display
    clusters_to_display <- input$expression_clusters_to_display
    cells_to_display <- which(
        grepl(sample_data()$cells$sample, pattern = paste0("^", samples_to_display, "$", collapse = "|")) == TRUE & 
        grepl(sample_data()$cells$cluster, pattern = paste0("^", clusters_to_display, "$", collapse = "|")) == TRUE
      )
    to_plot <- cbind(
        sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
        sample_data()$cells[ cells_to_display , ]
      )
    if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
      to_plot$level <- 0
    } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
      to_plot$level <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ]
    } else {
      to_plot$level <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ] %>%
        colMeans() %>%
        as.vector()
    }
    if ( input$expression_plotting_order == "Random" ) {
      to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
    } else if ( input$expression_plotting_order == "Highest expression on top" ) {
      to_plot <- to_plot[ order(to_plot$level, decreasing = FALSE) , ]
    }
    scatterD3::scatterD3(
      x = to_plot[ , 1 ],
      y = to_plot[ , 2 ],
      xlab = colnames(to_plot)[ 1 ],
      ylab = colnames(to_plot)[ 2 ],
      xlim = c(
          input$expression_projection_scale_x_manual_range[1],
          input$expression_projection_scale_x_manual_range[2]
        ),
      ylim = c(
          input$expression_projection_scale_y_manual_range[1],
          input$expression_projection_scale_y_manual_range[2]
        ),
      point_size = input$expression_projection_dot.size,
      col_var = to_plot$level,
      col_lab = "Gene expression",
      col_continuous = TRUE,
      point_opacity = input$expression_projection_opacity,
      transitions = FALSE,
      legend_width = 0,
      menu = FALSE,
      tooltip_text = paste0(
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br/>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br/>",
        "<b>nUMI</b>: ", to_plot[ , "nUMI" ], "<br/>",
        "<b>nGene</b>: ", to_plot[ , "nGene" ], "<br/>",
        "<b>Expr. MT</b>: ", format(to_plot[ , "percent_mt"   ]*100, digits = 1), "%<br/>",
        "<b>Expr. ribo</b>: ", format(to_plot[ , "percent_ribo" ]*100, digits = 1), "%<br/>"))
  })

  observeEvent(input$expression_projection_info, {
    showModal(
      modalDialog(
        title = "Dimensional reduction", easyClose = TRUE, footer = NULL,
        p("Interactive projection of cells into 2-dimensional space based on their expression profile.", 
          tags$ul(tags$li("Both tSNE and UMAP are frequently used algorithms for dimensional reduction in single cell transcriptomics. While they generally allow to make similar conclusions, some differences exist between the two (please refer to Google)."),
                  tags$li("Cell color reflects the log-normalised expression of entered genes. If more than 1 gene is entered, the color reflects the average expression of all genes. Genes must be in separate lines or separated by a space, comma, or semicolon. Reported below the projection are the genes that are present and absent in this data set. Absent genes could either have been annotated with a different name or were not expressed in any of the cells. Matching of gene names is case-insensitive, that means Myc/MYC/myc are treated equally."),
                  tags$li("Samples and clusters can be removed from the plot individually to highlight a contrast of interest."),
                  tags$li("Cells can be plotted either randomly (which a more unbiased image) or in the order of expression (with highest expression plotted last), sometimes resulting in a more appealing figure."),
                  tags$li("By default, the dot size is set to 15 without any transparency but both these attributes can be changed using the sliders on the left."),
                  tags$li("The last 2 slider elements on the left can be used to resize the projection axes. This can be particularly useful when a projection contains a population of cell that is very far away from the rest and therefore creates a big empty space (which is not uncommon for UMAPs).")
          ),
          "The plot is interactive (drag and zoom) but depending on the computer of the user and the number of cells displayed it can become very slow."
        )
      )
    )
  })

  observeEvent(input$expression_projection_export, {
    library("ggplot2")

    projection_to_display <- input$expression_projection_to_display
    samples_to_display <- input$expression_samples_to_display
    clusters_to_display <- input$expression_clusters_to_display
    cells_to_display <- which(
        grepl(
          sample_data()$cells$sample,
          pattern = paste0("^", samples_to_display, "$", collapse = "|")
        ) &
        grepl(
          sample_data()$cells$cluster,
          pattern = paste0("^", clusters_to_display, "$", collapse = "|")
        )
      )
    to_plot <- cbind(
        sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
        sample_data()$cells[ cells_to_display , ]
      )

    xlim <- c(
        input$expression_projection_scale_x_manual_range[1],
        input$expression_projection_scale_x_manual_range[2]
      )
    ylim <- c(
        input$expression_projection_scale_y_manual_range[1],
        input$expression_projection_scale_y_manual_range[2]
      )

    if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
      to_plot$level <- 0
      out_filename <- paste0(
          plot_export_path, "Cerebro_",
          sample_data()$experiment$experiment_name, "_gene_expression_none"
        )
    } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
      to_plot$level <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ]
      out_filename <- paste0(
          plot_export_path, "Cerebro_",
          sample_data()$experiment$experiment_name, "_gene_expression_",
          genesToPlot()$genes_to_display_present, "_",
          input$expression_projection_to_display
        )
    } else {
      to_plot$level <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ] %>%
        colMeans() %>%
        as.vector()
      out_filename <- paste0(
          plot_export_path, "Cerebro_",
          sample_data()$experiment$experiment_name, "_gene_expression_",
          genesToPlot()$genes_to_display_present[1],
          "_and_others_", input$expression_projection_to_display
        )
    }

    if ( input$expression_plotting_order == "Random" ) {
      to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
      out_filename <- paste0(out_filename, "_random_order.pdf")
    } else if ( input$expression_plotting_order == "Highest expression on top" ) {
      to_plot <- to_plot[ order(to_plot$level, decreasing = FALSE) , ]
      out_filename <- paste0(out_filename, "_highest_expression_on_top.pdf")
    }

    p <- ggplot(
        to_plot,
        aes_q(
          x = as.name(colnames(to_plot)[1]),
          y = as.name(colnames(to_plot)[2]),
          colour = as.name("level")
        )
      ) +
      geom_point() +
      viridis::scale_colour_viridis(
        name = "Log-normalised\nexpression",
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      ) +
      lims(x = xlim, y = ylim) +
      theme_bw()

    pdf(NULL)
    ggsave(out_filename, p, height = 8, width = 11)

    if ( file.exists(out_filename) ) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Success!",
        text = paste0("Plot saved successfully as: ", out_filename),
        type = "success"
      )
    } else {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Error!",
        text = "Sorry, it seems something went wrong...",
        type = "error"
      )
    }
  })

  ##--------------------------------------------------------------------------##

  # box plot by sample
  output$expression_by_sample <- plotly::renderPlotly({
    if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
      expression_levels <- 0
    } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
      expression_levels <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , ]
    } else {
      expression_levels <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , ] %>%
        colMeans() %>%
        as.vector()
    }
    data.frame(
      "sample" = sample_data()$cells[ , "sample" ],
      "expression" = expression_levels
    ) %>%
    plotly::plot_ly(
      x = ~sample,
      y = ~expression,
      type = "box",
      color = ~sample,
      colors = sample_data()$samples$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Expression level", hoverformat = ".4f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$expression_by_sample_info, {
    showModal(
      modalDialog(
        title = "Expression levels by sample", easyClose = TRUE, footer = NULL,
        p("Log-normalised expression of genes inserted above by sample. If more than 1 gene was given, this reflects the average across all cells of each sample.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot by cluster
  output$expression_by_cluster <- plotly::renderPlotly({
    if ( length(genesToPlot()$genes_to_display_present) == 0 ) {
      expression_levels <- 0
    } else if ( length(genesToPlot()$genes_to_display_present) == 1 ) {
      expression_levels <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , ]
    } else {
      expression_levels <- genesToPlot()$genes_to_display_present %>%
        sample_data()$expression[ . , ] %>%
        colMeans() %>%
        as.vector()
    }
    data.frame(
      "cluster" = sample_data()$cells[ , "cluster" ],
      "expression" = expression_levels
    ) %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~expression,
      type = "box",
      color = ~cluster,
      colors = sample_data()$clusters$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Expression level", hoverformat = ".4f"),
      dragmode =  "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$expression_by_cluster_info, {
    showModal(
      modalDialog(
        title = "Expression levels by cluster", easyClose = TRUE, footer = NULL,
        p("Log-normalised expression of genes inserted above by cluster. If more than 1 gene was given, this reflects the average across all cells of each cluster.")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Gene set expression
  ##--------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$parameters$organism
  ## - sample_data()$expression
  ## - sample_data()$projections
  ## - sample_data()$samples$overview$sample
  ## - sample_data()$clusters$overview$cluster
  ## - sample_data()$cells$sample
  ## - sample_data()$cells$cluster
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent_mt
  ## - sample_data()$cells$percent_ribo
  ##--------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##

  # reactive data
  geneSets <- reactive({
    if ( sample_data()$experiment$organism == "mm" ) {
      msigdbr::msigdbr(species = "Mus musculus")
    } else if ( sample_data()$experiment$organism == "hg" ) {
      msigdbr::msigdbr(species = "Homo sapiens")
    } else {
      msigdbr::msigdbr(species = "Mus musculus")
    }
  })

  # reactive data
  geneSetData <- reactive({
    geneSetData <- list()
    if ( is.null(input$geneSetexpression_select_geneSet) || is.na(input$geneSetexpression_select_geneSet) || input$geneSetexpression_select_geneSet == "-" ) {
      geneSetData$genes_to_display_present <- NULL
    } else {
      geneSetData$genes_to_display <- geneSets()[ which(geneSets()$gs_name == input$geneSetexpression_select_geneSet) , "gene_symbol" ]$gene_symbol
      geneSetData$genes_to_display <- unique(geneSetData$genes_to_display)
      geneSetData$genes_to_display_here <- rownames(sample_data()$expression)[ match(tolower(geneSetData$genes_to_display), tolower(rownames(sample_data()$expression))) ]
      geneSetData$genes_to_display_present <- geneSetData$genes_to_display_here[ which(!is.na(geneSetData$genes_to_display_here)) ]
      geneSetData$genes_to_display_missing <- geneSetData$genes_to_display[ which(is.na(geneSetData$genes_to_display_here)) ]
    }
    geneSetData
  })

  # show which genes are in the data set and which aren"t
  output$geneSetexpression_genes_displayed <- renderText({
    genes_to_display_text <- paste0(
        "<br><b>Total unique genes in gene set:</b><br>",
        length(geneSetData()$genes_to_display),
        "<br><b>Showing expression for ",
        length(geneSetData()$genes_to_display_present),
        " genes:</b><br>",
        paste0(geneSetData()$genes_to_display_present, collapse = ", "),
        "<br><b>",
        length(geneSetData()$genes_to_display_missing),
        " gene(s) are not in data set: </b><br>",
        paste0(geneSetData()$genes_to_display_missing, collapse = ", ")
      )
    if ( sample_data()$experiment$organism != "mm" & sample_data()$experiment$organism != "hg" ) {
      paste0(
        '<br><b><font color="red">Note:</b> Gene sets are available for human and mouse only. Organism for loaded samples is either not set or none of the two. Mouse gene sets are loaded and can be used.</font><br>',
        genes_to_display_text
      )
    } else {
      genes_to_display_text  
    }
  })

  ##--------------------------------------------------------------------------##

  # UI
  output$geneSetexpression_UI <- renderUI({
    tagList(
      selectInput(
        "geneSetexpression_projection_to_display", label = "Projection:",
        choices = names(sample_data()$projections)
      ),
      selectInput(
        "geneSetexpression_select_geneSet", label = "Gene set:",
        choices = c("-", unique(geneSets()$gs_name)), selected = "-"
      ),
      shinyWidgets::pickerInput(
        "geneSetexpression_samples_to_display", label = "Samples to display:",
        choices = sample_data()$samples$overview$sample,
        selected = sample_data()$samples$overview$sample,
        options = list("actions-box"=TRUE), multiple = TRUE
      ),
      shinyWidgets::pickerInput(
        "geneSetexpression_clusters_to_display", label = "Clusters to display:",
        choices = sample_data()$clusters$overview$cluster,
        selected = sample_data()$clusters$overview$cluster,
        options = list("actions-box"=TRUE), multiple = TRUE
      ),
      selectInput(
        "geneSetexpression_plotting_order", label = "Plotting order:",
        choices = c("Random", "Highest expression on top")
      ),
      sliderInput(
        "geneSetexpression_projection_dot.size", label = "Point size:",
        min = 0, max = 50, value = 25, step = 1
      ),
      sliderInput(
        "geneSetexpression_projection_opacity", label = "Point opacity:",
        min = 0, max = 1, value = 1, step = 0.05
      )
    )
  })

  output$geneSetexpression_scales <- renderUI({
    projection_to_display <- if ( is.null(input$geneSetexpression_projection_to_display) || is.na(input$geneSetexpression_projection_to_display) ) names(sample_data()$projections)[1] else input$geneSetexpression_projection_to_display
    range_x_min <- round(
      min(sample_data()$projections[[ projection_to_display ]][,1]) * 1.1)
    range_x_max <- round(
      max(sample_data()$projections[[ projection_to_display ]][,1]) * 1.1)
    range_y_min <- round(
      min(sample_data()$projections[[ projection_to_display ]][,2]) * 1.1)
    range_y_max <- round(
      max(sample_data()$projections[[ projection_to_display ]][,2]) * 1.1)
    tagList(
      sliderInput(
        "geneSetexpression_projection_scale_x_manual_range", label = "X axis",
        min = range_x_min, max = range_x_max,
        value = c(range_x_min, range_x_max)
      ),
      sliderInput(
        "geneSetexpression_projection_scale_y_manual_range", label = "Y axis",
        min = range_y_min, max = range_y_max,
        value = c(range_y_min, range_y_max)
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # projection with scatterD3
  output$geneSetexpression_projection <- scatterD3::renderScatterD3({
    req(input$geneSetexpression_projection_to_display)
    req(input$geneSetexpression_samples_to_display)
    req(input$geneSetexpression_clusters_to_display)
    req(input$geneSetexpression_plotting_order)
    projection_to_display <- input$geneSetexpression_projection_to_display
    samples_to_display <- input$geneSetexpression_samples_to_display
    clusters_to_display <- input$geneSetexpression_clusters_to_display
    cells_to_display <- which(
        grepl(
          sample_data()$cells$sample,
          pattern = paste0("^", samples_to_display, "$", collapse = "|")
        ) &
        grepl(
          sample_data()$cells$cluster,
          pattern = paste0("^", clusters_to_display, "$", collapse = "|")
        )
      )
    to_plot <- cbind(
        sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
        sample_data()$cells[ cells_to_display , ]
      )

    if ( length(geneSetData()$genes_to_display_present) == 0 ) {
      to_plot$level <- 0
    } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
      to_plot$level <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ]
    } else {
      to_plot$level <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ] %>%
        colMeans() %>%
        as.vector()
    }

    if ( input$geneSetexpression_plotting_order == "Random" ) {
      to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
    } else if ( input$geneSetexpression_plotting_order == "Highest expression on top" ) {
      to_plot <- to_plot[ order(to_plot$level, decreasing=FALSE) , ]
    }
    scatterD3::scatterD3(
      x = to_plot[ , 1 ],
      y = to_plot[ , 2 ],
      xlab = colnames(to_plot)[ 1 ],
      ylab = colnames(to_plot)[ 2 ],
      xlim = c(
          input$geneSetexpression_projection_scale_x_manual_range[1],
          input$geneSetexpression_projection_scale_x_manual_range[2]
        ),
      ylim = c(
          input$geneSetexpression_projection_scale_y_manual_range[1],
          input$geneSetexpression_projection_scale_y_manual_range[2]
        ),
      point_size = input$geneSetexpression_projection_dot.size,
      col_var = to_plot$level,
      col_lab = "Gene expression",
      col_continuous = TRUE,
      point_opacity = input$geneSetexpression_projection_opacity,
      transitions = FALSE,
      legend_width = 0,
      menu = FALSE,
      tooltip_text = paste0(
        "<b>Sample</b>: ", to_plot[ , "sample" ], "<br/>",
        "<b>Cluster</b>: ", to_plot[ , "cluster" ], "<br/>",
        "<b>nUMI</b>: ", to_plot[ , "nUMI" ], "<br/>",
        "<b>nGene</b>: ", to_plot[ , "nGene" ], "<br/>",
        "<b>Expr. MT</b>: ", format(to_plot[ , "percent_mt" ]*100, digits = 1), "%<br/>",
        "<b>Expr. ribo</b>: ", format(to_plot[ , "percent_ribo" ]*100, digits = 1), "%<br/>"))
  })

  observeEvent(input$geneSetexpression_projection_info, {
    showModal(
      modalDialog(
        title = "Dimensional reduction", easyClose = TRUE, footer = NULL,
        p("Interactive projection of cells into 2-dimensional space based on their expression profile.", 
          tags$ul(tags$li("Both tSNE and UMAP are frequently used algorithms for dimensional reduction in single cell transcriptomics. While they generally allow to make similar conclusions, some differences exist between the two (please refer to Google)."),
                  tags$li("For human and murine data sets, all organism-specific gene sets from the MSigDB can be selected. If the experiment was performed in another organism, the murine gene sets will be available."),
                  tags$li("Cell color reflects the average log-normalised expression of the genes in the selected gene set. Reported below the projection are the genes that are present and absent in this data set. Absent genes could either have been annotated with a different name or were not expressed in any of the cells. Matching of gene names is case-insensitive, that means Myc/MYC/myc are treated equally."),
                  tags$li("Samples and clusters can be removed from the plot individually to highlight a contrast of interest."),
                  tags$li("Cells can be plotted either randomly (which a more unbiased image) or in the order of expression (with highest expression plotted last), sometimes resulting in a more appealing figure."),
                  tags$li("By default, the dot size is set to 15 without any transparency but both these attributes can be changed using the sliders on the left."),
                  tags$li("The last 2 slider elements on the left can be used to resize the projection axes. This can be particularly useful when a projection contains a population of cell that is very far away from the rest and therefore creates a big empty space (which is not uncommon for UMAPs).")
          ),
          "The plot is interactive (drag and zoom) but depending on the computer of the user and the number of cells displayed it can become very slow."
        )
      )
    )
  })

  observeEvent(input$geneSetexpression_projection_export, {
    library("ggplot2")

    projection_to_display <- input$geneSetexpression_projection_to_display
    samples_to_display <- input$geneSetexpression_samples_to_display
    clusters_to_display <- input$geneSetexpression_clusters_to_display
    cells_to_display <- which(
        grepl(
          sample_data()$cells$sample,
          pattern = paste0("^", samples_to_display, "$", collapse = "|")
        ) &
        grepl(
          sample_data()$cells$cluster,
          pattern = paste0("^", clusters_to_display, "$", collapse = "|")
        )
      )
    to_plot <- cbind(
        sample_data()$projections[[ projection_to_display ]][ cells_to_display , ],
        sample_data()$cells[ cells_to_display , ]
      )

    xlim <- c(
        input$geneSetexpression_projection_scale_x_manual_range[1],
        input$geneSetexpression_projection_scale_x_manual_range[2]
      )
    ylim <- c(
        input$geneSetexpression_projection_scale_y_manual_range[1],
        input$geneSetexpression_projection_scale_y_manual_range[2]
      )

    if ( length(geneSetData()$genes_to_display_present) == 0 ) {
      to_plot$level <- 0
    } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
      to_plot$level <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ]
    } else {
      to_plot$level <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , cells_to_display ] %>%
        colMeans() %>%
        as.vector()
    }

    out_filename <- paste0(
        plot_export_path, "Cerebro_", sample_data()$experiment$experiment_name,
        "_gene_set_expression_", input$geneSetexpression_select_geneSet, "_",
        input$geneSetexpression_projection_to_display
      )

    if ( input$geneSetexpression_plotting_order == "Random" ) {
      to_plot <- sample(1:nrow(to_plot), nrow(to_plot)) %>% to_plot[ . , ]
      out_filename <- paste0(out_filename, "_random_order.pdf")
    } else if ( input$geneSetexpression_plotting_order == "Highest expression on top" ) {
      to_plot <- to_plot[ order(to_plot$level, decreasing = FALSE) , ]
      out_filename <- paste0(out_filename, "_highest_expression_on_top.pdf")
    }

    p <- ggplot(
        to_plot,
        aes_q(
          x = as.name(colnames(to_plot)[1]),
          y = as.name(colnames(to_plot)[2]),
          colour = as.name("level")
        )
      ) +
      geom_point() +
      viridis::scale_colour_viridis(
        name = "Average\nlog-normalised\nexpression",
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      lims(x = xlim, y = ylim) +
      theme_bw()

    pdf(NULL)
    ggsave(out_filename, p, height = 8, width = 11)

    if (file.exists(out_filename)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Success!",
        text = paste0("Plot saved successfully as: ", out_filename),
        type = "success"
      )
    } else {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = "Error!",
        text = "Sorry, it seems something went wrong...",
        type = "error"
      )
    }
  })

  ##--------------------------------------------------------------------------##

  # box plot by sample
  output$geneSetexpression_by_sample <- plotly::renderPlotly({
    if ( length(geneSetData()$genes_to_display_present) == 0 ) {
      expression_levels <- 0
    } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
      expression_levels <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , ]
    } else {
      expression_levels <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , ] %>%
        colMeans() %>%
        as.vector()
    }
    data.frame(
      "sample" = sample_data()$cells[ , "sample" ],
      "expression" = expression_levels
    ) %>%
    plotly::plot_ly(
      x = ~sample,
      y = ~expression,
      type = "box",
      color = ~sample,
      colors = sample_data()$samples$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Expression level", hoverformat=".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$geneSetexpression_by_sample_info, {
    showModal(
      modalDialog(
        title = "Expression levels by sample", easyClose = TRUE, footer = NULL,
        p("Average log-normalised expression of genes in selected gene set by sample.")
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot by cluster
  output$geneSetexpression_by_cluster <- plotly::renderPlotly({
    if ( length(geneSetData()$genes_to_display_present) == 0 ) {
      expression_levels <- 0
    } else if ( length(geneSetData()$genes_to_display_present) == 1 ) {
      expression_levels <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , ]
    } else {
      expression_levels <- geneSetData()$genes_to_display_present %>%
        sample_data()$expression[ . , ] %>%
        colMeans() %>%
        as.vector()
    }
    data.frame(
      "cluster" = sample_data()$cells[ , "cluster" ],
      "expression" = expression_levels
    ) %>%
    plotly::plot_ly(
      x = ~cluster,
      y = ~expression,
      type = "box",
      color = ~cluster,
      colors = sample_data()$clusters$colors,
      source = "subset",
      showlegend = FALSE,
      hoverinfo = "y",
      marker = list(size = 5)
    ) %>%
    plotly::layout(
      title = "",
      xaxis = list(title = ""),
      yaxis = list(title = "Expression level", hoverformat = ".2f"),
      dragmode = "select",
      hovermode = "compare"
    )
  })

  observeEvent(input$geneSetexpression_by_cluster_info, {
    showModal(
      modalDialog(
        title = "Expression levels by cluster", easyClose = TRUE, footer = NULL,
        p("Average log-normalised expression of genes in selected gene set by cluster.")
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Gene id/symbol conversion.
  ##--------------------------------------------------------------------------##
  output$gene_info <- DT::renderDataTable({
    if ( input$geneIdConversion_organism == "mouse" ) {
      conversion_table <- read.table("resources/mm10_gene_ID_name.txt",
        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else if ( input$geneIdConversion_organism == "human" ) {
      conversion_table <- read.table("resources/hg38_gene_ID_name.txt",
        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    }
    DT::datatable(
      conversion_table,
      filter = "none",
      selection = "multiple",
      escape = FALSE,
      autoHideNavigation = TRUE,
      rownames = FALSE,
      options = list(
        scrollX = FALSE,
        dom = "Bfrtip",
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 50
      )
    )
  })

  ##--------------------------------------------------------------------------##
  ## Panel: Sample info.
  ##--------------------------------------------------------------------------##
  ## Expected data:
  ## - see below
  ##--------------------------------------------------------------------------##
  output$sample_info_general <- renderText({
    info <- paste0(
        "<ul>",
          "<li><b>Experiment name:</b> ",
          sample_data()$experiment$experiment_name,
          "<li><b>Organism:</b> ",
          sample_data()$experiment$organism,
          "<li><b>Samples:</b> ",
          paste0(sample_data()$samples$overview$sample, collapse = ", "),
          "<li><b># of clusters found:</b> ",
          nrow(sample_data()$clusters$overview),
          "<li><b>Min. cells:</b> ",
          sample_data()$parameters$discard_genes_expressed_in_fewer_cells_than,
          "<li><b>Keep mitochondrial genes:</b> ",
          sample_data()$parameters$keep_mitochondrial_genes,
          "<li><b>Variables to regress:</b> ",
          sample_data()$parameters$variables_to_regress_out,
          "<li><b>Number of principal components:</b> ",
          sample_data()$parameters$number_PCs,
          "<li><b>tSNE perplexity:</b> ",
          sample_data()$parameters$tSNE_perplexity,
          "<li><b>Cluster resolution:</b> ",
          sample_data()$parameters$cluster_resolution,
          "<li><b>enrichR databases:</b> ",
          paste0(sample_data()$parameters$enrichr_dbs, collapse = ", "),
          "<li><b>Mitochondrial genes:</b> ",
          paste0(sample_data()$gene_lists$mitochondrial_genes, collapse = ", "),
          "<li><b>Ribosomal genes:</b> ",
          paste0(sample_data()$gene_lists$ribosomal_genes, collapse = ", "),
          "<li><b>Genes used for apoptotic score:</b> ",
          paste0(sample_data()$gene_lists$apoptosis_genes, collapse = ", "),
          "<li><b>S phase genes:</b> ",
          paste0(sample_data()$gene_lists$S_phase_genes, collapse = ", "),
          "<li><b>G2M phase genes:</b> ",
          paste0(sample_data()$gene_lists$G2M_phase_genes, collapse = ", "),
          "<li><b>Min/max # of UMI:</b> ",
          paste0(
            sample_data()$parameters$filtering$UMI_min, " / ",
            sample_data()$parameters$filtering$UMI_max),
          "<li><b>Min/max # of expressed genes:</b> ",
          paste0(
            sample_data()$parameters$filtering$genes_min, " / ",
            sample_data()$parameters$filtering$genes_max),
          "<li><b>Min/max % of mitochonrial genes:</b> ",
          paste0(
            sample_data()$parameters$filtering$percent_mt_min, " / ",
            sample_data()$parameters$filtering$percent_mt_min),
          "<li><b>Min/max % of ribosomal genes:</b> ",
          paste0(
            sample_data()$parameters$filtering$percent_ribo_min, " / ",
            sample_data()$parameters$filtering$percent_ribo_min),
        "</ul>"
      )
    info_R_raw <- sample_data()$technical_info$R
    info_R <- c()
    for ( i in 1:length(info_R_raw) ) {
      info_R <- paste(info_R, "<br>", info_R_raw[i])
    }
    paste0(
      info,
      "<br><b>R environment and packages use in analysis:</b><br><pre>",
      info_R,
      "</pre>"
    )
  })

  output$sample_info_R <- renderPrint({
    if ( !is.null(sample_data()$technical_info$R) ) {
      capture.output(sample_data()$technical_info$R)
    } else {
      print("Not available")
    }
  })


  ##--------------------------------------------------------------------------##
  ## Panel: About.
  ##--------------------------------------------------------------------------##
  output$about <- renderText({
    '<b>Version:</b><br>
     1.0 (January 2019)<br>
     <br>
     <b>Author:</b><br>
     Roman Hillje<br>
     Department of Experimental Oncology<br>
     IEO, European Institute of Oncology IRCCS, Milan<br>
     <br>
     <b>Contact:</b><br>
     <a href="mailto:roman.hillje@ieo.it?subject=Cerebro">roman.hillje@ieo.it</a><br>
     <br>
     <u>Do not share this application without permission.</u><br>
     <br>
     <b>Credit where credit is due:</b><br>
     <ul>
      <li>App icon made by <a href="https://www.flaticon.com/authors/kiranshastry" title="Kiranshastry" target="_blank">Kiranshastry</a> from <a href="https://www.flaticon.com/" title="Flaticon" target="_blank">www.flaticon.com</a> is licensed by <a href="http://creativecommons.org/licenses/by/3.0/" title="Creative Commons BY 3.0" target="_blank">CC 3.0 BY</a></li>
      <li>Sample and cluster color palettes taken from <a href="https://flatuicolors.com/" title="Flat UI Colors 2" target="_blank">https://flatuicolors.com/</a></li>
     </ul>'
  })
}


##----------------------------------------------------------------------------##
## UI.
##----------------------------------------------------------------------------##
ui <- dashboardPage(
  dashboardHeader(title = "Cerebro"),
  dashboardSidebar(
    tags$head(tags$style(HTML(".content-wrapper {overflow-x: scroll;}"))),
    sidebarMenu(
      sidebarMenuOutput("sidebar_menu")
    )
  ),
  dashboardBody(
    tags$script(HTML('$("body").addClass("fixed");')),
    tabItems(
      tabItem(tabName = "loadData",
        fluidRow(
          column(12,
            titlePanel("Load data"),
            fileInput(
              inputId = "RDS_file", label = "Choose RDS file...",
              multiple = FALSE, accept = c(".rds"), width = NULL,
              buttonLabel = "Browse...", placeholder = "No file selected")
          )
        )
      ),
      tabItem(tabName = "overview",
        tagList(
          fluidRow(
            column(width = 3, offset = 0, style = "padding: 0px;",
              box(
                title = "Input parameters", status = "primary",
                solidHeader = TRUE, width = 12, collapsible = TRUE,
                  tagList(
                    uiOutput("overview_UI"),
                    uiOutput("overview_scales")
                  )
              )
            ),
            column(width = 9, offset = 0, style = "padding: 0px;",
              box(
                title = tagList(
                  p("Dimensional reduction",
                    style = "padding-right: 5px; display: inline"
                  ),
                  actionButton(
                    inputId = "overview_projection_info", label = "info",
                    icon = NULL, class = "btn-xs",
                    title = "Show additional information for this panel.",
                    style = "margin-right: 5px"
                  ),
                  actionButton(
                    inputId = "overview_projection_export",
                    label = "export to PDF", icon = NULL, class = "btn-xs",
                    title = "Export dimensional reduction to PDF file."
                  )
                ),
                status = "primary", solidHeader = TRUE, width = 12,
                collapsible = TRUE,
                scatterD3::scatterD3Output(
                  "overview.projection", height = "720px"
                )
              )
            )
          )
        )
      ),
      tabItem(tabName = "samples",
        box(
          title = tagList(
            p("Overview of samples",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "samples_table_info", label = "info", icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          DT::dataTableOutput("samples_overview")
        ),
        box(
          title = tagList(
            p("Samples by cluster",
              style = "padding-right: 5px; display: inline"),
            actionButton(
              inputId = "samples_by_cluster_info",
              label = "info", icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          DT::dataTableOutput("samples_by_cluster_table"),
          uiOutput("samples_by_cluster_UI")
        ),
        box(
          title = tagList(
            p("Number of transcripts",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "samples_box_nUMI_info", label = "info", icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("samples_box_nUMI")
        ),
        box(
          title = tagList(
            p("Number of expressed genes",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "samples_box_nGene_info", label = "info", icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("samples_box_nGene")
        ),
        box(
          title = tagList(
            p("Mitochondrial gene expression",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "samples_box_percent_mt_info",
              label = "info", icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("samples_box_percent_mt")
        ),
        box(
          title = tagList(
            p("Ribosomal gene expression",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "samples_box_percent_ribo_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("samples_box_percent_ribo")
        ),
        box(
          title = tagList(
            p("Cell cycle analysis (Regev)",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "samples_by_cell_cycle_Regev_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("samples_by_cell_cycle_Regev_UI")
        ),
        box(
          title = tagList(
            p("Cell cycle analysis (Cyclone)",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "samples_by_cell_cycle_Cyclone_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("samples_by_cell_cycle_Cyclone_UI")
        )
      ),
      tabItem(tabName = "clusters",
        box(
          title = tagList(
            p("Cluster tree", style = "padding-right: 5px; display: inline"),
            actionButton(
              inputId = "clusters_tree_info", label = "info", icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("clusters_tree_UI")
        ),
        box(
          title = tagList(
            p(
              "Clusters by samples",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "clusters_by_sample_info", label = "info", icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          DT::dataTableOutput("clusters_by_sample"),
          uiOutput("clusters_by_sample_UI")
        ),
        box(
          title = tagList(
            p(
              "Number of transcripts",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "clusters_box_nUMI_info", label = "info", icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("clusters_box_nUMI")
        ),
        box(
          title = tagList(
            p(
              "Number of expressed genes",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "clusters_box_nGene_info", label = "info", icon = NULL,
              class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("clusters_box_nGene")
        ),
        box(
          title = tagList(
            p(
              "Mitochondrial gene expression",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "clusters_box_percent_mt_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("clusters_box_percent_mt")
        ),
        box(
          title = tagList(
            p(
              "Ribosomal gene expression",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "clusters_box_percent_ribo_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          plotly::plotlyOutput("clusters_box_percent_ribo")
        ),
        box(
          title = tagList(
            p(
              "Cell cycle analysis (Regev)",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "clusters_by_cell_cycle_Regev_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("clusters_by_cell_cycle_Regev_UI")
        ),
        box(
          title = tagList(
            p(
              "Cell cycle analysis (Cyclone)",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "clusters_by_cell_cycle_Cyclone_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("clusters_by_cell_cycle_Cyclone_UI")
        )
      ),
      tabItem(tabName = "topExpressedGenes",
        box(
          title = tagList(
            p(
              "Top expressed genes per sample",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "top_expressed_genes_by_sample_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("top_expressed_genes_by_sample_UI")
        ),
        box(
          title = tagList(
            p(
              "Top expressed genes per cluster",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "top_expressed_genes_by_cluster_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("top_expressed_genes_by_cluster_UI")
        )
      ),
      tabItem(tabName = "markerGenes",
        box(
          title = tagList(
            p("Marker genes per sample",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "marker_genes_by_sample_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("marker_genes_by_sample_UI")
        ),
        box(
          title = tagList(
            p("Marker genes per cluster",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "marker_genes_by_cluster_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("marker_genes_by_cluster_UI")
        )
      ),
      tabItem(tabName = "enrichedPathways",
        box(
          title = tagList(
            p("Enriched pathways by sample",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "enriched_pathways_by_sample_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("enriched_pathways_by_sample_UI")
        ),
        box(
          title = tagList(
            p("Enriched pathways by cluster",
              style = "padding-right: 5px; display: inline"
            ),
            actionButton(
              inputId = "enriched_pathways_by_cluster_info", label = "info",
              icon = NULL, class = "btn-xs",
              title = "Show additional information for this panel."
            )
          ),
          status = "primary", solidHeader = TRUE, width = 12,
          collapsible = TRUE,
          uiOutput("enriched_pathways_by_cluster_UI")
        )
      ),
      tabItem(tabName = "geneExpression",
        tagList(
          fluidRow(
            column(width = 3, offset = 0, style = "padding: 0px;",
              box(
                title = "Input parameters", status = "primary",
                solidHeader = TRUE, width = 12, collapsible = TRUE,
                tagList(
                  uiOutput("expression_UI"),
                  uiOutput("expression_scales")
                )
              )
            ),
            column(width = 9, offset = 0, style = "padding: 0px;",
              box(
                title = tagList(
                  p("Dimensional reduction",
                    style = "padding-right: 5px; display: inline"
                  ),
                  actionButton(
                    inputId = "expression_projection_info", label = "info",
                    icon = NULL, class = "btn-xs",
                    title = "Show additional information for this panel.",
                    style = "margin-right: 5px"
                  ),
                  actionButton(
                    inputId = "expression_projection_export",
                    label = "export to PDF", icon = NULL, class = "btn-xs",
                    title = "Export dimensional reduction to PDF file."
                  )
                ),
                status = "primary", solidHeader = TRUE, width = 12,
                collapsible = TRUE,
                tagList(
                  scatterD3::scatterD3Output(
                    "expression_projection", height = "720px"
                  ),
                  tags$br(),
                  htmlOutput("expression_genes_displayed")
                )
              )
            )
          ),
          fluidRow(
            box(
              title = tagList(
                p("Expression levels by sample",
                  style = "padding-right: 5px; display: inline"
                ),
                actionButton(
                  inputId = "expression_by_sample_info", label = "info",
                  icon = NULL, class = "btn-xs",
                  title = "Show additional information for this panel."
                )
              ),
              status = "primary", solidHeader = TRUE, width = 12,
              collapsible = TRUE,
              plotly::plotlyOutput("expression_by_sample")
            )
          ),
          fluidRow(
            box(
              title = tagList(
                p("Expression levels by cluster",
                  style = "padding-right: 5px; display: inline"
                ),
                actionButton(
                  inputId = "expression_by_cluster_info", label = "info",
                  icon = NULL, class = "btn-xs",
                  title = "Show additional information for this panel."
                )
              ),
              status = "primary", solidHeader = TRUE, width = 12,
              collapsible = TRUE,
              plotly::plotlyOutput("expression_by_cluster")
            )
          )
        )
      ),
      tabItem(tabName = "geneSetExpression",
        tagList(
          fluidRow(
            column(width = 3, offset = 0, style = "padding:0px;",
              box(
                title = "Input parameters", status = "primary",
                solidHeader = TRUE, width = 12, collapsible = TRUE,
                tagList(
                  uiOutput("geneSetexpression_UI"),
                  uiOutput("geneSetexpression_scales")
                )
              )
            ),
            column(width = 9, offset = 0, style = "padding:0px;",
              box(
                title = tagList(
                  p("Dimensional reduction",
                    style = "padding-right: 5px; display: inline"
                  ),
                  actionButton(
                    inputId = "geneSetexpression_projection_info",
                    label = "info", icon = NULL, class = "btn-xs",
                    title = "Show additional information for this panel.",
                    style = "margin-right: 5px"
                  ),
                  actionButton(
                    inputId = "geneSetexpression_projection_export",
                    label = "export to PDF", icon = NULL, class = "btn-xs",
                    title = "Export dimensional reduction to PDF file."
                  )
                ),
                status = "primary", solidHeader = TRUE, width = 12,
                collapsible = TRUE,
                tagList(
                  scatterD3::scatterD3Output(
                    "geneSetexpression_projection", height = "720px"
                  ),
                  tags$br(),
                  htmlOutput("geneSetexpression_genes_displayed")
                )
              )
            )
          ),
          fluidRow(
            box(
              title = tagList(
                p("Expression levels by sample",
                  style = "padding-right: 5px; display: inline"
                ),
                actionButton(
                  inputId = "geneSetexpression_by_sample_info", label = "info",
                  icon = NULL, class = "btn-xs",
                  title = "Show additional information for this panel."
                )
              ),
              status = "primary", solidHeader = TRUE, width = 12,
              collapsible = TRUE,
              plotly::plotlyOutput("geneSetexpression_by_sample")
            )
          ),
          fluidRow(
            box(
              title = tagList(
                p("Expression levels by cluster",
                style = "padding-right: 5px; display: inline"
                ),
                actionButton(
                  inputId = "geneSetexpression_by_cluster_info", label = "info",
                  icon = NULL, class = "btn-xs",
                  title = "Show additional information for this panel."
                )
              ),
              status = "primary", solidHeader = TRUE, width = 12,
              collapsible = TRUE,
              plotly::plotlyOutput("geneSetexpression_by_cluster")
            )
          )
        )
      ),
      tabItem(tabName = "geneIdConversion",
        box(
          title = "Convert gene ID <-> gene symbol", status = "primary",
          solidHeader = TRUE, width = 12, collapsible = FALSE,
          tagList(
            selectInput(
              "geneIdConversion_organism", "Organism:",
              choices = c("mouse", "human")
            ),
            DT::dataTableOutput("gene_info")
          )
        )
      ),
      tabItem(tabName = "info",
        fluidPage(
          fluidRow(
            column(12,
              titlePanel("Information about samples and analysis"),
              htmlOutput("sample_info_general")
            )
          )
        )
      ),
      tabItem(tabName = "about",
        fluidPage(
          fluidRow(
            column(12,
              titlePanel("About this application"),
              htmlOutput("about")
            )
          )
        )
      )
    )
  )
)


##----------------------------------------------------------------------------##
## App.
##----------------------------------------------------------------------------##
shinyApp(ui, server)


