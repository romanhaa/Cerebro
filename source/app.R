##----------------------------------------------------------------------------##
## single cell browser
## version 1.0.3
##
## Author:    Roman Hillje
## Institute: IEO
## Lab:       PGP
## Date:      2018-08-18
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Data to be loaded is stored in this structure.
##----------------------------------------------------------------------------##
## Header
## - version                            e.g. '1.1'
## - type                               'transcriptomics'
##
## Specific info
## - parameters                         list
##          - project.name              character
##          - organism                  e.g. 'mm'
##          - reference.version         e.g. 'mm10'
##          - annotation                e.g. 'GENCODE'
##          - annotation.type           name or ID
##          - enrichr.dbs               list of databases queried using enrichR
##          - min.genes                 integer
##          - vars.to.regress           vector
##          - numberPCs                 integer
##          - tsne.perplexity           integer
## - gene.lists                         list
##          - genes.mt                  vector
##          - genes.ribo                vector
##          - genes.apoptosis           vector
##          - genes.S                   vector
##          - genes.G2M                 vector
## - samples                            list
##          - names                     named vector
##          - count                     integer
##          - names.long                vector
##          - table.before.filtering    data frame
##          - table.after.filtering     data frame
## - clusters                           list
##          - names                     vector
##          - count                     integer
## - cells:                             data frame
##          - coordinates of tSNE       data frame (2 columns x # of cells)
##          - coordinates of UMAP       data frame (2 columns x # of cells)
##          - cell indices              vector
##          - sample                    vector
##          - cluster                   vector
##          - nUMI                      vector
##          - nGene                     vector
##          - percent.mt                vector
##          - percent.ribo              vector
##          - apoptotic score           vector
##          - cell cycle (Regev)        vector
##          - cell cycle (Cyclone)      vector
## - most.expressed.genes:              list
##          - by.sample:                data frame
##          - by.cluster:               data frame
## - marker.genes:                      list
##          - by.sample:                data frame
##          - by.sample.annotation:     list of data frames
##          - by.cluster:               data frame
##          - by.cluster.annotation:    list of data frames
## - expression:                        sparse expression matrix (seurat@data)
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
required_packages_CRAN <- c('shiny', 'shinydashboard', 'DT', 'plotly', 'reshape',
                            'scales', 'RColorBrewer', 'Matrix', 'scatterD3',
                            'ggplot2', 'msigdbr', 'shinyWidgets', 'formattable')
for ( package in required_packages_CRAN ) {
  if ( is.element(el=package, set=rownames(installed.packages())) == FALSE ) {
    install.packages(package, repos='http://cran.us.r-project.org', dependencies=TRUE)
  }# else {
  #  library(package, character.only=TRUE)
  #}
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
required_packages_bioconductor <- c('ggtree')
for ( package in required_packages_bioconductor ) {
  if ( is.element(el=package, set=rownames(installed.packages())) == FALSE ) {
    source('https://bioconductor.org/biocLite.R')
    biocLite(package, suppressAutoUpdate=TRUE, suppressUpdates=TRUE)
  }
}

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
source('resources/descriptions.txt')

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
library('shiny')
library('shinydashboard')
library('formattable')
#library('Cairo')
#library('rlang')
library('dplyr')

##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
system('type R')
print(.libPaths())
print(sessionInfo())
print(Sys.getenv())
Sys.setenv('R_LIBS_USER' = '')
print(Sys.getenv())


##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
if (grepl(tolower(Sys.info()['sysname']), pattern='^win')) {
  plot.export.path <- paste0(Sys.getenv('USERPROFILE'), '\\Desktop\\')
} else {
  plot.export.path <- '~/Desktop/'
}


##----------------------------------------------------------------------------##
## Check if data already exists.
##----------------------------------------------------------------------------##
if ( file.exists('resources/data.rds') ) {
  mode <- 'boxed'
} else {
  mode <- 'open'
}


##----------------------------------------------------------------------------##
## Allow upload of files up to 400 MB.
##----------------------------------------------------------------------------##
options(shiny.maxRequestSize=400*1024^2) 


##----------------------------------------------------------------------------##
## Color management.
##----------------------------------------------------------------------------##
# Dutch palette from flatuicolors.com
colors.dutch <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
                  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
                  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
                  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51')

# Spanish palette from flatuicolors.com
colors.spanish <- c('#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
                    '#2c2c54','#474787','#aaa69d','#227093','#218c74',
                    '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
                    '#b33939','#cd6133','#84817a','#cc8e35','#ccae62')

colors <- c(colors.dutch, colors.spanish)

# colors <- c('#BEAED4','#FDC086','#BF5B17','#386CB0','#F0027F','#666666','#1B9E77','#7FC97F','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#FFFF99','#666666','#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928','#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#F2F2F2','#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC','#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3','#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F')
cell.cycle.colorset <- setNames(c('#45aaf2','#f1c40f','#e74c3c', '#7f8c8d'), c('G1','S','G2M','-'))


##----------------------------------------------------------------------------##
## Server.
##----------------------------------------------------------------------------##
server <- function(input, output, session) {

  ##----------------------------------------------------------------------------##
  ## Sidebar menu.
  ##----------------------------------------------------------------------------##
  output$sidebar.menu <- renderMenu({
    if ( mode == 'open' ) {
      sidebarMenu(id='sidebar',
        menuItem('Load data', tabName='loadData', icon=icon('spinner'), selected=TRUE),
        menuItem('Overview', tabName='overview', icon=icon('binoculars')),
        menuItem('Samples', tabName='samples', icon=icon('star')),
        menuItem('Clusters', tabName='clusters', icon=icon('braille')),
        menuItem('Top expressed genes', tabName='topExpressedGenes', icon=icon('bullhorn')),
        menuItem('Marker genes', tabName='markerGenes', icon=icon('magnet')),
        menuItem('Enriched pathways', tabName='enrichedPathways', icon=icon('sitemap')),
        menuItem('Gene expression', tabName='geneExpression', icon=icon('signal')),
        menuItem('Gene set expression', tabName='geneSetExpression', icon=icon('list')),
        menuItem('Gene ID conversion', tabName='geneIdConversion', icon=icon('barcode')),
        menuItem('Sample info', tabName='info', icon=icon('info')),
        menuItem('About', tabName='about', icon=icon('at'))
      )
    } else {
      sidebarMenu(id='sidebar',
        menuItem('Load data', tabName='loadData', icon=icon('spinner')),
        menuItem('Overview', tabName='overview', icon=icon('binoculars'), selected=TRUE),
        menuItem('Samples', tabName='samples', icon=icon('star')),
        menuItem('Clusters', tabName='clusters', icon=icon('braille')),
        menuItem('Top expressed genes', tabName='topExpressedGenes', icon=icon('bullhorn')),
        menuItem('Marker genes', tabName='markerGenes', icon=icon('magnet')),
        menuItem('Enriched pathways', tabName='enrichedPathways', icon=icon('sitemap')),
        menuItem('Gene expression', tabName='geneExpression', icon=icon('signal')),
        menuItem('Gene set expression', tabName='geneSetExpression', icon=icon('list')),
        menuItem('Gene ID conversion', tabName='geneIdConversion', icon=icon('barcode')),
        menuItem('Sample info', tabName='info', icon=icon('info')),
        menuItem('About', tabName='about', icon=icon('at'))
      )
    }
  })


  ##----------------------------------------------------------------------------##
  ## Sample data.
  ##----------------------------------------------------------------------------##
  sample_data <- reactive({

    if ( mode == 'boxed' ) {
      sample_data <- readRDS('resources/data.rds')
    } else if ( is.null(input$RDS_file) || is.na(input$RDS_file) ) {
      sample_data <- readRDS('resources/example.rds')
    } else {
      req(input$RDS_file)
      input_file  <- input$RDS_file
      sample_data <- readRDS(input_file$datapath)
      # updateTabsetPanel(session, inputId='sidebar', selected='overview')
    }


    ##----------------------------------------------------------------------------##
    ## Check if data set is old and, if so, convert it to new format.
    ##----------------------------------------------------------------------------##
    if ( is.null(sample_data$cells) ) {

      #
      sample_data$projections <- list('tSNE' = sample_data$tSNE[ , c(1,2) ])
      sample_data$cells <- sample_data$tSNE[ , c('sampleID','cluster','nUMI','nGene','pct.mt','pct.ribo') ]
      colnames(sample_data$cells) <- c('sample','cluster','nUMI','nGene','percent.mt','percent.ribo')
      sample_data$tSNE <- NULL

      #
      sample_data$samples <- list('names' = levels(sample_data$cells$sample),
                                  'count' = length(levels(sample_data$cells$sample)))

      #
      sample_data$clusters <- list('names' = levels(sample_data$cells$cluster),
                                   'count' = length(levels(sample_data$cells$cluster)))

      # #
      # temp.table           <- as.data.frame(table(sample_data$cells$sample))
      # colnames(temp.table) <- c('sample', 'total_cell_count')
      # temp.table$sample    <- as.character(temp.table$sample)

      # for ( i in sample_data$clusters$names ) {
      #   temp.counts           <- as.data.frame(table(sample_data$cells$sample[ which(sample_data$cells$cluster == i) ]))
      #   temp.counts[,1]       <- NULL
      #   colnames(temp.counts) <- i
      #   temp.table            <- cbind(temp.table, temp.counts)
      # }

      # sample_data$samples$by.cluster <- temp.table

      # #
      # temp.table           <- as.data.frame(table(sample_data$cells$cluster))
      # colnames(temp.table) <- c('cluster', 'total_cell_count')
      # temp.table$cluster   <- as.character(temp.table$cluster)

      # for ( i in sample_data$samples$names ) {
      #   temp.counts           <- as.data.frame(table(sample_data$cells$cluster[ which(sample_data$cells$sample == i) ]))
      #   temp.counts[,1]       <- NULL
      #   colnames(temp.counts) <- c(i)
      #   temp.table            <- cbind(temp.table, temp.counts)
      # }

      # sample_data$clusters$by.sample <- temp.table

      sample_data$table.samples <- NULL
      sample_data$table.clusters <- NULL

      #
      sample_data$most.expressed.genes <- list(
        'by.sample' = sample_data$most.expressed.genes.samples,
        'by.cluster' = sample_data$most.expressed.genes.clusters)

      sample_data$most.expressed.genes.samples <- NULL
      sample_data$most.expressed.genes.clusters <- NULL

      #
      sample_data$marker.genes <- list(
        'by.sample'             = sample_data$marker.genes.samples,
        'by.cluster'            = sample_data$marker.genes.clusters,
        'by.sample.annotation'  = NULL,
        'by.cluster.annotation' = NULL)

      sample_data$marker.genes.samples <- NULL
      sample_data$marker.genes.clusters <- NULL

      #
      sample_data$expression <- sample_data$expression_table
      sample_data$expression_table <- NULL

      #
      sample_data$data.version <- '1.0'
      sample_data$data.type <- 'transcriptomics'
      sample_data$parameters <- list(
        'project.name'      = 'not available',
        'min.genes'         = 'not available',
        'organism'          = 'not available',
        'reference.version' = 'not available',
        'annotation'        = 'not available',
        'annotation.type'   = 'not available',
        'enrichr.dbs'       = 'not available',
        'vars.to.regress'   = 'not available',
        'number.PCs'        = 'not available',
        'tSNE.perplexity'   = 'not available')

      #
      sample_data$gene.lists <- list(
        'genes.mt'        = 'not available',
        'genes.ribo'      = 'not available',
        'genes.apoptosis' = 'not available',
        'genes.S'         = 'not available',
        'genes.G2M'       = 'not available')

      }

    ##----------------------------------------------------------------------------##
    ## Replace NA values in Cyclone cell cycle analysis with '-' if there are any.
    ##----------------------------------------------------------------------------##
    if ( 'cell.cycle.Cyclone' %in% colnames(sample_data$cells) ) {
      if ( length(which(is.na(sample_data$cells$cell.cycle.Cyclone))) > 0 ) {
        sample_data$cells$cell.cycle.Cyclone <- factor(sample_data$cells$cell.cycle.Cyclone, levels=c('G1','S','G2M','-'))
        sample_data$cells$cell.cycle.Cyclone[ which(is.na(sample_data$cells$cell.cycle.Cyclone)) ] <- '-'
      }
    }

    return(sample_data)
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Overview.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$projections
  ## - sample_data()$cells$sample
  ## - sample_data()$cells$cluster
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent.mt
  ## - sample_data()$cells$percent.ribo
  ## - sample_data()$samples$names
  ## - sample_data()$samples$count
  ## - sample_data()$clusters$names
  ## - sample_data()$clusters$count
  ##----------------------------------------------------------------------------##
  output$overview.UI <- renderUI({
    tagList(
      selectInput('overview.projection.to.display', label='Projection',
        choices = names(sample_data()$projections)),
      selectInput('overview.projection.cell.color', label='Color cells by',
        choices = names(sample_data()$cells)[! names(sample_data()$cells) %in% c('cell.barcode')]),
      checkboxInput('overview.projection.ellipses', label='Confidence ellipses',
        value = FALSE),
      shinyWidgets::pickerInput('overview.samples.to.display', label='Samples to display',
        choices  = unname(sample_data()$samples$names),
        selected = unname(sample_data()$samples$names),
        options  = list('actions-box'=TRUE),
        multiple = TRUE),
      shinyWidgets::pickerInput('overview.clusters.to.display', label='Clusters to display',
        choices  = sample_data()$clusters$names,
        selected = sample_data()$clusters$names,
        options  = list('actions-box'=TRUE),
        multiple = TRUE),
      selectInput('overview.projection.cell.size.variable', label='Change point size by',
        choices  = c('None', 'nUMI', 'nGene'),
        selected = 'None'),
      sliderInput('overview.projection.cell.size', label='Point size',
        min=0, max=50, value=15, step=1),
      sliderInput('overview.projection.cell.opacity', label='Point opacity',
        min=0, max=1, value=1, step=0.05)
    )
  })

  output$overview.scales <- renderUI({
    projection.to.display <- if ( is.null(input$overview.projection.to.display) || is.na(input$overview.projection.to.display) ) names(sample_data()$projections)[1] else input$overview.projection.to.display
    range.x.min <- round(min(sample_data()$projections[[ projection.to.display ]][,1])*1.1)
    range.x.max <- round(max(sample_data()$projections[[ projection.to.display ]][,1])*1.1)
    range.y.min <- round(min(sample_data()$projections[[ projection.to.display ]][,2])*1.1)
    range.y.max <- round(max(sample_data()$projections[[ projection.to.display ]][,2])*1.1)
    tagList(
      sliderInput('overview.projection.scale.x.manual.range', label='X axis',
        min=range.x.min, max=range.x.max, value=c(range.x.min, range.x.max)),
      sliderInput('overview.projection.scale.y.manual.range', label='Y axis',
        min=range.y.min, max=range.y.max, value=c(range.y.min, range.y.max))
    )
  })

  output$overview.projection <- scatterD3::renderScatterD3({
    req(input$overview.projection.to.display)
    req(input$overview.samples.to.display)
    req(input$overview.clusters.to.display)
    req(input$overview.projection.scale.x.manual.range)
    samples.cols          <- colors[ 1:sample_data()$samples$count ]
    samples.colorset      <- setNames(samples.cols, unname(sample_data()$samples$names))
    clusters.cols         <- colors[ 1:sample_data()$clusters$count ]
    clusters.colorset     <- setNames(clusters.cols, unname(sample_data()$clusters$names))
    projection.to.display <- if ( is.null(input$overview.projection.to.display) || is.na(input$overview.projection.to.display) ) names(sample_data()$projections)[1] else input$overview.projection.to.display
    samples.to.display    <- if ( is.null(input$overview.samples.to.display) || is.na(input$overview.samples.to.display) ) sample_data()$samples$names else input$overview.samples.to.display
    clusters.to.display   <- if ( is.null(input$overview.clusters.to.display) || is.na(input$overview.clusters.to.display) ) sample_data()$clusters$names else input$overview.clusters.to.display
    cells.to.display      <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', samples.to.display, '$', collapse='|')) == TRUE & 
                                   grepl(sample_data()$cells$cluster, pattern=paste0('^', clusters.to.display, '$', collapse='|')) == TRUE)
    to.plot               <- cbind(sample_data()$projections[[ projection.to.display ]][ cells.to.display , ],
                                   sample_data()$cells[ cells.to.display , ])
    to.plot               <- to.plot[ sample(1:nrow(to.plot)) , ]
    col.var               <- to.plot[ , input$overview.projection.cell.color ]
    colors                <- if ( is.null(input$overview.projection.cell.color) || is.na(input$overview.projection.cell.color) ) {
                               NULL
                             } else if ( input$overview.projection.cell.color == 'sample' ) {
                               samples.colorset
                             } else if ( input$overview.projection.cell.color == 'cluster' ) {
                               clusters.colorset
                             } else if ( input$overview.projection.cell.color %in% c('cell.cycle.Regev','cell.cycle.Cyclone') ) {
                               cell.cycle.colorset
                             } else if ( is.factor(to.plot[,input$overview.projection.cell.color]) ) {
                               setNames(colors[1:length(levels(to.plot[,input$overview.projection.cell.color]))], levels(to.plot[,input$overview.projection.cell.color]))
                             } else if ( is.character(to.plot[,input$overview.projection.cell.color]) ) {
                               colors
                             } else {
                               NULL
                             }
    size.var              <- if ( input$overview.projection.cell.size.variable == 'None' ) NULL else to.plot[ , input$overview.projection.cell.size.variable ]
    scatterD3::scatterD3(
      x             = to.plot[ , 1 ],
      y             = to.plot[ , 2 ],
      xlab          = colnames(to.plot)[ 1 ],
      ylab          = colnames(to.plot)[ 2 ],
      xlim          = c(input$overview.projection.scale.x.manual.range[1], input$overview.projection.scale.x.manual.range[2]),
      ylim          = c(input$overview.projection.scale.y.manual.range[1], input$overview.projection.scale.y.manual.range[2]),
      point_size    = input$overview.projection.cell.size,
      col_var       = col.var,
      col_lab       = input$overview.projection.cell.color,
      colors        = colors,
      ellipses      = input$overview.projection.ellipses,
      size_var      = size.var,
      point_opacity = input$overview.projection.cell.opacity,
      transitions   = FALSE,
      menu          = FALSE,
      tooltip_text  = paste0(
        '<b>Sample</b>: ', to.plot[ , 'sample'       ], '<br/>',
        '<b>Cluster</b>: ',   to.plot[ , 'cluster'      ], '<br/>',
        '<b>nUMI</b>: ',      to.plot[ , 'nUMI'         ], '<br/>',
        '<b>nGene</b>: ',     to.plot[ , 'nGene'        ], '<br/>',
        '<b>Expr. MT</b>: ',   format(to.plot[ , 'percent.mt'   ]*100, digits=1), '%<br/>',
        '<b>Expr. ribo</b>: ', format(to.plot[ , 'percent.ribo' ]*100, digits=1), '%<br/>'))
  })

  observeEvent(input$overview.projection.info, {
    showModal(modalDialog(overview.projection.info.text, title=overview.projection.info.title, easyClose=TRUE, footer=NULL))
  })

  observeEvent(input$overview.projection.export, {
    #volumes <- c(Home = fs::path_home())
    #print(volumes)
    #shinyFiles::shinyDirChoose(input, overview.projection.export.directory, updateFreq = 0, defaultPath = "", defaultRoot = NULL)
    #shinyFiles::shinyDirChoose(input, overview.projection.export.directory)
    #out.path <- shinyFiles::shinyDirChoose(input, 'overview.projection.export', roots=c(home='~'), session=session, restrictions=system.file(package='base'))

    library('ggplot2')

    projection.to.display <- input$overview.projection.to.display
    samples.to.display    <- input$overview.samples.to.display
    clusters.to.display   <- input$overview.clusters.to.display
    cells.to.display      <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', samples.to.display, '$', collapse='|')) == TRUE & grepl(sample_data()$cells$cluster, pattern=paste0('^', clusters.to.display, '$', collapse='|')) == TRUE)
    to.plot               <- cbind(sample_data()$projections[[ projection.to.display ]][ cells.to.display , ],
                                   sample_data()$cells[ cells.to.display , ])
    to.plot               <- to.plot[ sample(1:nrow(to.plot)) , ]

    xlim <- c(input$overview.projection.scale.x.manual.range[1], input$overview.projection.scale.x.manual.range[2])
    ylim <- c(input$overview.projection.scale.y.manual.range[1], input$overview.projection.scale.y.manual.range[2])

    if ( is.factor(to.plot[,input$overview.projection.cell.color]) | is.character(to.plot[,input$overview.projection.cell.color]) ) {
      if ( input$overview.projection.cell.color == 'sample' ) {
        cols <- setNames(colors[ 1:sample_data()$samples$count ], unname(sample_data()$samples$names))
      } else if ( input$overview.projection.cell.color == 'cluster' ) {
        cols <- setNames(colors[ 1:sample_data()$clusters$count ], unname(sample_data()$clusters$names))
      } else if ( input$overview.projection.cell.color %in% c('cell.cycle.Regev','cell.cycle.Cyclone') ) {
        cols <- cell.cycle.colorset
      } else if ( is.factor(to.plot[,input$overview.projection.cell.color]) ) {
        cols <- setNames(colors[1:length(levels(to.plot[,input$overview.projection.cell.color]))], levels(to.plot[,input$overview.projection.cell.color]))
      } else {
        cols <- colors
      }
      p <- ggplot(to.plot, aes_q(x=as.name(colnames(to.plot)[1]), y=as.name(colnames(to.plot)[2]), colour=as.name(input$overview.projection.cell.color))) +
           geom_point() +
           scale_colour_manual(values=cols) +
           lims(x=xlim, y=ylim) +
           theme_bw()
    } else {
      p <- ggplot(to.plot, aes_q(x=as.name(colnames(to.plot)[1]), y=as.name(colnames(to.plot)[2]), colour=as.name(input$overview.projection.cell.color))) +
           geom_point() +
           viridis::scale_colour_viridis(guide=guide_colorbar(frame.colour='black', ticks.colour='black')) +
           lims(x=xlim, y=ylim) +
           theme_bw()
    }

    out.filename <- paste0(plot.export.path, 'scBrowser_', gsub(sample_data()$parameters$project.name, pattern=' ', replacement='_'), '_overview_', input$overview.projection.to.display, '_by_', gsub(input$overview.projection.cell.color, pattern='\\.', replacement='_'), '.pdf')

    pdf(NULL)
    ggsave(out.filename, p, height=8, width=11)

    if (file.exists(out.filename)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = 'Success!',
        text = paste0('Plot saved successfully as: ', out.filename),
        type = 'success'
      )
    } else {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = 'Error!',
        text = "Sorry, it seems something went wrong...",
        type = 'error'
      )
    }
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Samples.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$samples$names
  ## - sample_data()$cells$sample
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent.mt
  ## - sample_data()$cells$percent.ribo
  ## - sample_data()$cells$cell.cycle.Regev (optional)
  ## - sample_data()$cells$cell.cycle.Cyclone (optional)
  ##----------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##
  ## cell counts by sample

  output$samples.table <- DT::renderDataTable({
    table           <- as.data.frame(table(sample_data()$cells$sample))
    colnames(table) <- c('sample', 'total_cell_count')
    table$sample    <- as.character(table$sample)
    for ( i in sample_data()$clusters$names ) {
      counts           <- as.data.frame(table(sample_data()$cells$sample[ which(sample_data()$cells$cluster == i) ]))
      counts[,1]       <- NULL
      colnames(counts) <- c(i)
      table            <- cbind(table, counts)
    }
    DT::datatable(
      data               = table,
      filter             = 'none',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      colnames           = c('Sample', '# of cells', colnames(table)[ 3:length(table) ]),
      rownames           = FALSE,
      class              = 'cell-border stripe',
      options            = list(
        scrollX    = TRUE,
        sDom       = '<"top">lrt<"bottom">ip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15
      )
    )
  })

  observeEvent(input$samples.table.info, {
    showModal(modalDialog(samples.table.info.text, title=samples.table.info.title, easyClose=TRUE, footer=NULL))
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of samples by cluster
  output$samples.by.cluster.UI <- renderUI(
    if ( sample_data()$clusters$count > 1 ) {
      plotly::plotlyOutput('samples.by.cluster.plot')
    } else {
      textOutput('samples.by.cluster.text')
    }
  )

  # bar plot of samples by cluster
  output$samples.by.cluster.plot <- plotly::renderPlotly({
    table           <- as.data.frame(table(sample_data()$cells$sample))
    colnames(table) <- c('sample', 'total_cell_count')
    table$sample    <- as.character(table$sample)
    for ( i in sample_data()$clusters$names ) {
      counts           <- as.data.frame(table(sample_data()$cells$sample[ which(sample_data()$cells$cluster == i) ]))
      counts[,1]       <- NULL
      colnames(counts) <- c(i)
      table            <- cbind(table, counts)
    }
    table <- table[,c('sample', 'total_cell_count', sample_data()$clusters$names)]
    temp <- reshape::melt(table[ , c(1,3:length(table)) ], id.vars='sample')
    colnames(temp) <- c('sample','cluster','cells')
    temp$sample <- factor(temp$sample, levels=sample_data()$samples$names)
    temp$total <- 0
    for ( i in sample_data()$samples$names ) {
      temp$total[ which(temp$sample == i) ] <- table$total_cell_count[ which(table$sample == i) ] #sample_data()$samples$by.cluster$total_cell_count[ which(sample_data()$samples$by.cluster$sample == i) ]
    }
    temp$pct <- temp$cells / temp$total
    plotly::plot_ly(temp,
      x         = ~sample,
      y         = ~pct*100,
      type      = 'bar',
      color     = ~cluster,
      colors    = colors[ 1:sample_data()$clusters$count ],
      text      = ~pct*100,
      hoverinfo = 'name+y') %>%
    plotly::layout(
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage (%)', hoverformat='.2f'),
      barmode   = 'stack',
      hovermode = 'compare') 
  })

  # # alternative text output for bar plot of samples by cluster
  output$samples.by.cluster.text <- renderText({ 'Only 1 cluster in this data set.' })

  observeEvent(input$samples.by.cluster.info, {
    showModal(
      modalDialog(
        title='Samples by cluster', easyClose=TRUE, footer=NULL,
        p('Percentage bar plot representation of the table shown above. Allows to see which clusters contribute most strongly to each sample. Clusters can be removed from the plot by clicking on them in the legend.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of transcripts per sample
  output$samples.box.nUMI <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[ , c('sample', 'nUMI') ],
      x          = ~sample,
      y          = ~nUMI,
      type       = 'box',
      color      = ~sample,
      colors     = colors[ 1:sample_data()$samples$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size = 5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Number of UMIs', type='log', hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$samples.box.nUMI.info, {
    showModal(
      modalDialog(
        title='Number of transcripts', easyClose=TRUE, footer=NULL,
        p('Box plot of the number of transcripts (UMIs) found in each sample.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of expressed genes per sample
  output$samples.box.nGene <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[ , c('sample', 'nGene') ],
      x          = ~sample,
      y          = ~nGene,
      type       = 'box',
      color      = ~sample,
      colors     = colors[ 1:sample_data()$samples$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Number of expressed genes', type='log', hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$samples.box.nGene.info, {
    showModal(
      modalDialog(
        title='Number of expressed genes', easyClose=TRUE, footer=NULL,
        p('Box plot of the number of expressed genes found in each sample.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of mitochondrial gene expression per sample
  output$samples.box.percent.mt <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[,c('sample', 'percent.mt')],
      x          = ~sample,
      y          = ~percent.mt,
      type       = 'box',
      color      = ~sample,
      colors     = colors[ 1:sample_data()$samples$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage of mitochondrial gene expression', range=c(0,1), hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$samples.box.percent.mt.info, {
    showModal(
      modalDialog(
        title='Mitochondrial gene expression', easyClose=TRUE, footer=NULL,
        p('Box plot of the percentage of mitochondrial gene expression found in each sample. This reflects the contribution of mitochondrial transcripts to the entire transcriptome in each cell. A list of all genes considered to be mitochondrial can be found in the "Sample info" tab on the left.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of ribosomal gene expression per sample
  output$samples.box.percent.ribo <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[,c('sample', 'percent.ribo')],
      x          = ~sample,
      y          = ~percent.ribo,
      type       = 'box',
      color      = ~sample,
      colors     = colors[ 1:sample_data()$samples$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage of ribosomal gene expression', range=c(0,1), hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$samples.box.percent.ribo.info, {
    showModal(
      modalDialog(
        title='Ribosomal gene expression', easyClose=TRUE, footer=NULL,
        p('Box plot of the percentage of ribosomal gene expression found in each sample. This reflects the contribution of ribosomal transcripts to the entire transcriptome in each cell. A list of all genes considered to be ribosomal can be found in the "Sample info" tab on the left.')
      )
    )
  })

  ##--------------------------------------------------------------------------##
  
  # UI element for bar plot of samples by cell cycle (Regev)
  output$samples.by.cell.cycle.Regev.UI <- renderUI(
    if ( !is.null(sample_data()$cells$cell.cycle.Regev) ) {
      plotly::plotlyOutput('samples.by.cell.cycle.Regev.plot')
    } else {
      textOutput('samples.by.cell.cycle.Regev.text')
    }
  )

  # bar plot of samples by cell cycle (Regev)  
  output$samples.by.cell.cycle.Regev.plot <- plotly::renderPlotly({
    table           <- as.data.frame(table(sample_data()$cells$sample))
    colnames(table) <- c('sample', 'total_cell_count')
    table$sample    <- as.character(table$sample)
    for ( i in sample_data()$samples$names ) {
      table[ which(table$sample == i) , 'G1'  ] <- length(sample_data()$cells$sample[ which(sample_data()$cells$sample == i & sample_data()$cells$cell.cycle.Regev == 'G1')  ])
      table[ which(table$sample == i) , 'S'   ] <- length(sample_data()$cells$sample[ which(sample_data()$cells$sample == i & sample_data()$cells$cell.cycle.Regev == 'S')   ])
      table[ which(table$sample == i) , 'G2M' ] <- length(sample_data()$cells$sample[ which(sample_data()$cells$sample == i & sample_data()$cells$cell.cycle.Regev == 'G2M') ])
    }
    temp <- reshape::melt(table[ , c(1,3:length(table)) ], id.vars='sample')
    colnames(temp) <- c('sample','phase','cells')
    temp$sample <- factor(temp$sample, levels=sample_data()$samples$names)
    temp$phase <- factor(temp$phase, levels=c('G1','S','G2M'))
    temp$total <- 0
    for ( i in table$sample ) {
      temp$total[ which(temp$sample == i) ] <- table$total_cell_count[ which(table$sample == i) ]
    }
    temp$pct <- temp$cells / temp$total
    plotly::plot_ly(temp,
      x         = ~sample,
      y         = ~pct*100,
      type      = 'bar',
      color     = ~phase,
      colors    = c('#45aaf2','#f1c40f','#e74c3c'),
      text      = ~pct*100,
      hoverinfo = 'name+y') %>%
    plotly::layout(
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage (%)', hoverformat='.2f'),
      barmode   = 'stack',
      hovermode = 'compare') 
  })

  # alternative text for bar plot of samples by cell cycle (Regev)  
  output$samples.by.cell.cycle.Regev.text <- renderText({ 'Data not available.' })

  observeEvent(input$samples.by.cell.cycle.Regev.info, {
    showModal(
      modalDialog(
        title='Cell cycle analysis (Regev)', easyClose=TRUE, footer=NULL,
        p('Cell cycle distribution by sample using the method embedded in the Seurat framework. For each cell, it calculates scores for both G2M and S phase based on lists of genes (see "Sample info" tab on the left) and assigns the cell cycle phase on the basis of these scores.')
      )
    )
  })

  ##--------------------------------------------------------------------------##
  
  # UI element for bar plot of samples by cell cycle (Cyclone)
  output$samples.by.cell.cycle.Cyclone.UI <- renderUI(
    if ( !is.null(sample_data()$cells$cell.cycle.Cyclone) ) {
      plotly::plotlyOutput('samples.by.cell.cycle.Cyclone.plot')
    } else {
      textOutput('samples.by.cell.cycle.Cyclone.text')
    }
  )

  # bar plot of samples by cell cycle (Cyclone)
  output$samples.by.cell.cycle.Cyclone.plot <- plotly::renderPlotly({
    table           <- as.data.frame(table(sample_data()$cells$sample))
    colnames(table) <- c('sample', 'total_cell_count')
    table$sample    <- as.character(table$sample)
    for ( i in sample_data()$samples$names ) {
      table[ which(table$sample == i) , 'G1'  ] <- length(sample_data()$cells$sample[ which(sample_data()$cells$sample == i & sample_data()$cells$cell.cycle.Cyclone == 'G1')  ])
      table[ which(table$sample == i) , 'S'   ] <- length(sample_data()$cells$sample[ which(sample_data()$cells$sample == i & sample_data()$cells$cell.cycle.Cyclone == 'S')   ])
      table[ which(table$sample == i) , 'G2M' ] <- length(sample_data()$cells$sample[ which(sample_data()$cells$sample == i & sample_data()$cells$cell.cycle.Cyclone == 'G2M') ])
      table[ which(table$sample == i) , '-'   ] <- length(sample_data()$cells$sample[ which(sample_data()$cells$sample == i & ( is.na(sample_data()$cells$cell.cycle.Cyclone) | sample_data()$cells$cell.cycle.Cyclone == '-') ) ])
    }
    temp           <- reshape::melt(table[ , c(1,3:length(table)) ], id.vars='sample')
    colnames(temp) <- c('sample','phase','cells')
    temp$sample    <- factor(temp$sample, levels=sample_data()$samples$names)
    temp$phase     <- factor(temp$phase, levels=c('G1','S','G2M','-'))
    temp$total     <- 0
    for ( i in table$sample ) {
      temp$total[ which(temp$sample == i) ] <- table$total_cell_count[ which(table$sample == i) ]
    }
    temp$pct <- temp$cells / temp$total
    plotly::plot_ly(temp,
      x         = ~sample,
      y         = ~pct*100,
      type      = 'bar',
      color     = ~phase,
      colors    = c('#45aaf2','#f1c40f','#e74c3c', '#7f8c8d'),
      text      = ~pct*100,
      hoverinfo = 'name+y') %>%
    plotly::layout(
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage (%)', hoverformat='.2f'),
      barmode   = 'stack',
      hovermode = 'compare') 
  })

  # alternative text for bar plot of samples by cell cycle (Cyclone)
  output$samples.by.cell.cycle.Cyclone.text <- renderText({ 'Data not available.' })

  observeEvent(input$samples.by.cell.cycle.Cyclone.info, {
    showModal(
      modalDialog(
        title='Cell cycle analysis (Cyclone)', easyClose=TRUE, footer=NULL,
        p('Cell cycle distribution by sample using the machine learning-based Cyclone method published by Scialdone et al (2015). It assigns the cell cycle phase based on scores calculated using relative expression of lists of gene pairs. In contrast to the Seurat/Regev method, scores are calculated for G1 and G2M phase. Cells with a low score for both are assigned S phase. Inability to predict the cell cycle phase for a given cell with this method is most likely a result of very few expressed genes in the respective cell.')
      )
    )
  }) 


  ##----------------------------------------------------------------------------##
  ## Panel: Clusters.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$samples$count
  ## - sample_data()$clusters$names
  ## - sample_data()$clusters$tree
  ## - sample_data()$cells$cluster
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent.mt
  ## - sample_data()$cells$percent.ribo
  ## - sample_data()$cells$cell.cycle.Regev (optional)
  ## - sample_data()$cells$cell.cycle.Cyclone (optional)
  ##----------------------------------------------------------------------------##

  ##--------------------------------------------------------------------------##
  # cell counts by cluster and sample
  output$clusters.table <- DT::renderDataTable({
    table           <- as.data.frame(table(sample_data()$cells$cluster))
    colnames(table) <- c('cluster', 'total_cell_count')
    table$cluster   <- as.character(table$cluster)
    for ( i in sample_data()$samples$names ) {
      counts           <- as.data.frame(table(sample_data()$cells$cluster[ which(sample_data()$cells$sample == i) ]))
      counts[,1]       <- NULL
      colnames(counts) <- c(i)
      table            <- cbind(table, counts)
    }
    if ( sample_data()$samples$count > 1 ) {
      DT::datatable(
        data               = table,
        filter             = 'none',
        selection          = 'multiple',
        escape             = FALSE,
        autoHideNavigation = TRUE,
        colnames           = c('Cluster', '# of cells', colnames(table)[ 3:length(table) ]),
        rownames           = FALSE,
        class              = 'cell-border stripe',
        options            = list(
          scrollX    = TRUE,
          sDom       = '<"top">lrt<"bottom">ip',
          lengthMenu = c(20, 30, 50, 100),
          pageLength = 20)
      )
    } else {
      DT::datatable(
        data               = table,
        filter             = 'none',
        selection          = 'multiple',
        escape             = FALSE,
        autoHideNavigation = TRUE,
        colnames           = c('Cluster', '# of cells'),
        rownames           = FALSE,
        class              = 'cell-border stripe',
        options            = list(
          scrollX    = TRUE,
          sDom       = '<"top">lrt<"bottom">ip',
          lengthMenu = c(20, 30, 50, 100),
          pageLength = 20)
      )
    }
  })

  observeEvent(input$clusters.table.info, {
    showModal(
      modalDialog(
        title='Overview of clusters', easyClose=TRUE, footer=NULL,
        p('Table of clusters (by row) stratified by sample (columns).')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for cluster tree
  output$clusters.tree.UI <- renderUI(
    if ( !is.null(sample_data()$clusters$tree) ) {
      plotOutput('clusters.tree.plot')
    } else {
      textOutput('clusters.tree.text')
    }
  )

  output$clusters.tree.plot <- renderPlot({
    library('ggtree')
    tree <- sample_data()$clusters$tree
    tree$tip.label <- paste0('Cluster ', tree$tip.label)
    colors.tree <- colors[1:length(tree$tip.label)]
    ggplot(tree, aes(x, y)) + 
      ggplot2::scale_y_reverse() +
      xlim(0, max(tree$edge.length*1.1)) +
      geom_tree() +
      theme_tree() +
      geom_tiplab(size=5, hjust=-0.2) +
      geom_tippoint(color=colors.tree, shape=16, size=6)
  })

  output$clusters.tree.text <- renderText({ 'Data not available.' })

  observeEvent(input$clusters.tree.info, {
    showModal(
      modalDialog(
        title='Cluster tree', easyClose=TRUE, footer=NULL,
        p('Cluster tree reflecting the similarity of clusters based on their expression profiles. Instead of using the expression values, the correlation is calculated using the user-specified number of principal components (see "Sample info" tab on the left).')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of clusters by sample
  output$clusters.by.sample.UI <- renderUI(
    if ( sample_data()$samples$count > 1 ) {
      plotly::plotlyOutput('clusters.by.sample.plot')
    } else {
      textOutput('clusters.by.sample.text')
    }
  )

  # bar plot of clusters by sample
  output$clusters.by.sample.plot <- plotly::renderPlotly({
    table           <- as.data.frame(table(sample_data()$cells$cluster))
    colnames(table) <- c('cluster', 'total_cell_count')
    table$cluster   <- as.character(table$cluster)
    for ( i in sample_data()$samples$names ) {
      counts           <- as.data.frame(table(sample_data()$cells$cluster[ which(sample_data()$cells$sample == i) ]))
      counts[,1]       <- NULL
      colnames(counts) <- c(i)
      table            <- cbind(table, counts)
    }
    table <- table[,c('cluster', 'total_cell_count', sample_data()$samples$names)]
    temp <- reshape::melt(table[ , c(1,3:length(table)) ], id.vars='cluster')
    colnames(temp) <- c('cluster','sample','cells')
    temp$cluster <- factor(temp$cluster, levels=sample_data()$clusters$names)
    temp$total <- 0
    for ( i in sample_data()$clusters$names ) {
      temp$total[ which(temp$cluster == i) ] <- table$total_cell_count[ which(table$cluster == i) ] #sample_data()$clusters$by.sample$total_cell_count[ which(sample_data()$clusters$by.sample$cluster == i) ]
    }
    temp$pct <- temp$cells / temp$total
    plotly::plot_ly(temp,
      x         = ~cluster,
      y         = ~pct*100,
      type      = 'bar',
      color     = ~sample,
      colors    = colors[ 1:sample_data()$samples$count ],
      text      = ~pct*100,
      hoverinfo = 'name+y') %>%
    plotly::layout(
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage (%)', hoverformat='.2f'),
      barmode   = 'stack',
      hovermode = 'compare') 
  })

  # alternative text for bar plot of clusters by sample
  output$clusters.by.sample.text <- renderText({ 'Only 1 sample in this data set.' })

  observeEvent(input$clusters.by.sample.info, {
    showModal(
      modalDialog(
        title='Clusters by samples', easyClose=TRUE, footer=NULL,
        p('Percentage bar plot representation of the table shown above. Allows to see which samples contribute most strongly to each cluster. Samples can be removed from the plot by clicking on them in the legend.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of transcripts per cluster
  output$clusters.box.nUMI <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[ , c('cluster', 'nUMI') ],
      x          = ~cluster,
      y          = ~nUMI,
      type       = 'box',
      color      = ~cluster,
      colors     = colors[ 1:sample_data()$clusters$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''), 
      yaxis     = list(title='Number of UMIs', type='log', hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$clusters.box.nUMI.info, {
    showModal(
      modalDialog(
        title='Number of transcripts', easyClose=TRUE, footer=NULL,
        p('Box plot of the number of transcripts (UMIs) found in each cluster.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of number of expressed genes per cluster
  output$clusters.box.nGene <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[ , c('cluster', 'nGene') ],
      x          = ~cluster, 
      y          = ~nGene,
      type       = 'box',
      color      = ~cluster,
      colors     = colors[ 1:sample_data()$clusters$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Number of expressed genes', type='log', hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$clusters.box.nGene.info, {
    showModal(
      modalDialog(
        title='Number of expressed genes', easyClose=TRUE, footer=NULL,
        p('Box plot of the number of expressed genes found in each cluster.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of mitochondrial gene expression per cluster
  output$clusters.box.percent.mt <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[ , c('cluster', 'percent.mt') ],
      x          = ~cluster,
      y          = ~percent.mt,
      type       = 'box',
      color      = ~cluster,
      colors     = colors[ 1:sample_data()$clusters$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size = 5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage of mitochondrial gene expression', range=c(0, 1), hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$clusters.box.percent.mt.info, {
    showModal(
      modalDialog(
        title='Mitochondrial gene expression', easyClose=TRUE, footer=NULL,
        p('Box plot of the percentage of mitochondrial gene expression found in each cluster. This reflects the contribution of mitochondrial transcripts to the entire transcriptome in each cell. A list of all genes considered to be mitochondrial can be found in the "Sample info" tab on the left.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot of percentage of ribosomal gene expression per cluster
  output$clusters.box.percent.ribo <- plotly::renderPlotly({
    plotly::plot_ly(sample_data()$cells[ , c('cluster', 'percent.ribo') ],
      x          = ~cluster,
      y          = ~percent.ribo,
      type       = 'box',
      color      = ~cluster,
      colors     = colors[ 1:sample_data()$clusters$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage of ribosomal gene expression', range=c(0, 1), hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$clusters.box.percent.ribo.info, {
    showModal(
      modalDialog(
        title='Ribosomal gene expression', easyClose=TRUE, footer=NULL,
        p('Box plot of the percentage of ribosomal gene expression found in each cluster. This reflects the contribution of ribosomal transcripts to the entire transcriptome in each cell. A list of all genes considered to be ribosomal can be found in the "Sample info" tab on the left.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of clusters by cell cycle (Regev)
  output$clusters.by.cell.cycle.Regev.UI <- renderUI(
    if ( !is.null(sample_data()$cells$cell.cycle.Regev) ) {
      plotly::plotlyOutput('clusters.by.cell.cycle.Regev.plot')
    } else {
      textOutput('clusters.by.cell.cycle.Regev.text')
    }
  )

  # bar plot of clusters by cell cycle (Regev)
  output$clusters.by.cell.cycle.Regev.plot <- plotly::renderPlotly({
    table           <- as.data.frame(table(sample_data()$cells$cluster))
    colnames(table) <- c('cluster', 'total_cell_count')
    table$cluster   <- as.character(table$cluster)
    for ( i in sample_data()$clusters$names ) {
      table[ which(table$cluster == i) , 'G1'  ] <- length(sample_data()$cells$cluster[ which(sample_data()$cells$cluster == i & sample_data()$cells$cell.cycle.Regev == 'G1')  ])
      table[ which(table$cluster == i) , 'S'   ] <- length(sample_data()$cells$cluster[ which(sample_data()$cells$cluster == i & sample_data()$cells$cell.cycle.Regev == 'S')   ])
      table[ which(table$cluster == i) , 'G2M' ] <- length(sample_data()$cells$cluster[ which(sample_data()$cells$cluster == i & sample_data()$cells$cell.cycle.Regev == 'G2M') ])
    }
    temp <- reshape::melt(table[ , c(1,3:length(table)) ], id.vars='cluster')
    colnames(temp) <- c('cluster','phase','cells')
    temp$cluster <- factor(temp$cluster, levels=sample_data()$clusters$names)
    temp$phase <- factor(temp$phase, levels=c('G1','S','G2M'))
    temp$total <- 0

    for ( i in table$cluster ) {
      temp$total[ which(temp$cluster == i) ] <- table$total_cell_count[ which(table$cluster == i) ]
    }

    temp$pct <- temp$cells / temp$total

    plotly::plot_ly(temp,
      x         = ~cluster,
      y         = ~pct*100,
      type      = 'bar',
      color     = ~phase,
      colors    = c('#45aaf2','#f1c40f','#e74c3c'),
      text      = ~pct*100,
      hoverinfo = 'name+y') %>%
    plotly::layout(
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage (%)', hoverformat='.2f'),
      barmode   = 'stack',
      hovermode = 'compare') 
  })

  # alternative text for bar plot of clusters by cell cycle (Regev)
  output$clusters.by.cell.cycle.Regev.text <- renderText({ 'Data not available.' })

  observeEvent(input$clusters.by.cell.cycle.Regev.info, {
    showModal(
      modalDialog(
        title='Cell cycle analysis (Regev)', easyClose=TRUE, footer=NULL,
        p('Cell cycle distribution by cluster using the method embedded in the Seurat framework. For each cell, it calculates scores for both G2M and S phase based on lists of genes (see "Sample info" tab on the left) and assigns the cell cycle phase on the basis of these scores.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # UI element for bar plot of clusters by cell cycle (Cyclone)
  output$clusters.by.cell.cycle.Cyclone.UI <- renderUI(
    if ( !is.null(sample_data()$cells$cell.cycle.Cyclone) ) {
      plotly::plotlyOutput('clusters.by.cell.cycle.Cyclone.plot')
    } else {
      textOutput('clusters.by.cell.cycle.Cyclone.text')
    }
  )

  # bar plot of clusters by cell cycle (Cyclone)
  output$clusters.by.cell.cycle.Cyclone.plot <- plotly::renderPlotly({
    table           <- as.data.frame(table(sample_data()$cells$cluster))
    colnames(table) <- c('cluster', 'total_cell_count')
    table$cluster   <- as.character(table$cluster)
    for ( i in sample_data()$clusters$names ) {
      table[ which(table$cluster == i) , 'G1'  ] <- length(sample_data()$cells$cluster[ which(sample_data()$cells$cluster == i & sample_data()$cells$cell.cycle.Cyclone == 'G1')  ])
      table[ which(table$cluster == i) , 'S'   ] <- length(sample_data()$cells$cluster[ which(sample_data()$cells$cluster == i & sample_data()$cells$cell.cycle.Cyclone == 'S')   ])
      table[ which(table$cluster == i) , 'G2M' ] <- length(sample_data()$cells$cluster[ which(sample_data()$cells$cluster == i & sample_data()$cells$cell.cycle.Cyclone == 'G2M') ])
      table[ which(table$cluster == i) , '-'   ] <- length(sample_data()$cells$cluster[ which(sample_data()$cells$cluster == i & ( is.na(sample_data()$cells$cell.cycle.Cyclone) | sample_data()$cells$cell.cycle.Cyclone == '-') ) ])
    }
    temp <- reshape::melt(table[ , c(1,3:length(table)) ], id.vars='cluster')
    colnames(temp) <- c('cluster','phase','cells')
    temp$cluster <- factor(temp$cluster, levels=sample_data()$clusters$names)
    temp$phase <- factor(temp$phase, levels=c('G1','S','G2M','-'))
    temp$total <- 0

    for ( i in table$cluster ) {
      temp$total[ which(temp$cluster == i) ] <- table$total_cell_count[ which(table$cluster == i) ]
    }

    temp$pct <- temp$cells / temp$total

    plotly::plot_ly(temp, 
      x         = ~cluster,
      y         = ~pct*100,
      type      = 'bar',
      color     = ~phase,
      colors    = c('#45aaf2','#f1c40f','#e74c3c', '#7f8c8d'),
      text      = ~pct*100,
      hoverinfo = 'name+y') %>%
    plotly::layout(
      xaxis     = list(title=''),
      yaxis     = list(title='Percentage (%)', hoverformat='.2f'), 
      barmode   = 'stack',
      hovermode = 'compare')
  })

  # alternative text for bar plot of clusters by cell cycle (Cyclone)
  output$clusters.by.cell.cycle.Cyclone.text <- renderText({ 'Data not available.' })

  observeEvent(input$clusters.by.cell.cycle.Cyclone.info, {
    showModal(
      modalDialog(
        title='Cell cycle analysis (Cyclone)', easyClose=TRUE, footer=NULL,
        p('Cell cycle distribution by cluster using the machine learning-based Cyclone method published by Scialdone et al (2015). It assigns the cell cycle phase based on scores calculated using relative expression of lists of gene pairs. In contrast to the Seurat/Regev method, scores are calculated for G1 and G2M phase. Cells with a low score for both are assigned S phase.')
      )
    )
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Top expressed genes.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$samples$names
  ## - sample_data()$most.expressed.genes$by.sample
  ## - sample_data()$clusters$names
  ## - sample_data()$most.expressed.genes$by.cluster
  ##----------------------------------------------------------------------------##
  # by sample
  output$top.expressed.genes.by.sample.UI <- renderUI(
    if ( !is.null(sample_data()$most.expressed.genes$by.sample) ) {
      fluidRow(
        column(12,
          selectInput('top.expressed.genes.by.sample.input', label=NULL,
            choices = unname(sample_data()$samples$names)),
          DT::dataTableOutput('top.expressed.genes.by.sample.table.present')
        )
      )
    } else {
      textOutput('top.expressed.genes.by.sample.table.missing')
    }
  )

  output$top.expressed.genes.by.sample.table.present <- DT::renderDataTable(server=FALSE, {
    req(input$top.expressed.genes.by.sample.input)
    table <- sample_data()$most.expressed.genes$by.sample[ which(sample_data()$most.expressed.genes$by.sample$sample == input$top.expressed.genes.by.sample.input) , ]
    if (sum(table$pct > 1)) {
      table$pct <- table$pct / 100
    }
    table$pct <- round(table$pct, digits=4)
    colnames(table) <- c('Sample', 'Gene', '% of total expression')
    table$'% of total expression' <- formattable::percent(table$'% of total expression')
    table <- formattable::formattable(table, list(
      'Sample' = formattable::color_tile(colors[ which(sample_data()$samples$names == input$top.expressed.genes.by.sample.input) ],colors[ which(sample_data()$samples$names == input$top.expressed.genes.by.sample.input) ]),
      '% of total expression' = formattable::color_bar('pink')
    ))
    formattable::as.datatable(table,
      filter             = 'none',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      rownames           = FALSE,
      extensions         = c('Buttons'),
      class              = 'cell-border stripe',
      options            = list(
        dom        = 'Bfrtip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons    = list(
          'colvis', 
           list(
             extend  = 'collection',
             text    = 'Download',
             buttons = list(
              list(extend='csv',   filename='top_expressed_genes_per_sample', title='Top expressed genes per sample'),
              list(extend='excel', filename='top_expressed_genes_per_sample', title='Top expressed genes per sample'),
              list(extend='pdf',   filename='top_expressed_genes_per_sample', title='Top expressed genes per sample')
            )
          )
        )
      )
    ) %>% DT::formatStyle('% of total expression', textAlign='right')
  })

  output$top.expressed.genes.by.sample.table.missing <- renderText({ 'Data not available.' })

  observeEvent(input$top.expressed.genes.by.sample.info, {
    showModal(
      modalDialog(
        title='Top expressed genes per sample', easyClose=TRUE, footer=NULL,
        p('Table of top 100 most expressed genes in each sample. For example, if gene XY contributes with 5% to the total expression, that means 5% of all transcripts found in all cells of this sample come from that respective gene. These lists can help to identify/verify the dominant cell types.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # by cluster
  output$top.expressed.genes.by.cluster.UI <- renderUI(
    if ( !is.null(sample_data()$most.expressed.genes$by.cluster) ) {
      fluidRow(
        column(12,
          selectInput('top.expressed.genes.by.cluster.input', label=NULL, choices=sample_data()$clusters$names),
          DT::dataTableOutput('top.expressed.genes.by.cluster.table.present')
        )
      )
    } else {
      textOutput('top.expressed.genes.by.cluster.table.missing')
    }
  )

  output$top.expressed.genes.by.cluster.table.present <- DT::renderDataTable(server=FALSE, {
    req(input$top.expressed.genes.by.cluster.input)
    table <- sample_data()$most.expressed.genes$by.cluster[ which(sample_data()$most.expressed.genes$by.cluster$cluster == input$top.expressed.genes.by.cluster.input) , ]
    if (sum(table$pct > 1)) {
      table$pct <- table$pct / 100
    }
    table$pct <- round(table$pct, digits=4)
    colnames(table) <- c('Cluster', 'Gene', '% of total expression')
    table$'% of total expression' <- formattable::percent(table$'% of total expression')
    table <- formattable::formattable(table, list(
      'Cluster' = formattable::color_tile(colors[ which(sample_data()$clusters$names == input$top.expressed.genes.by.cluster.input) ],colors[ which(sample_data()$clusters$names == input$top.expressed.genes.by.cluster.input) ]),
      '% of total expression' = formattable::color_bar('pink')
    ))
    formattable::as.datatable(table,
      filter             = 'none',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      rownames           = FALSE,
      extensions         = c('Buttons'),
      class              = 'cell-border stripe',
      options            = list(
        dom        = 'Bfrtip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons    = list(
          'colvis', 
           list(
             extend  = 'collection',
             text    = 'Download',
             buttons = list(
              list(extend='csv',   filename='top_expressed_genes_per_cluster', title='Top expressed genes per cluster'),
              list(extend='excel', filename='top_expressed_genes_per_cluster', title='Top expressed genes per cluster'),
              list(extend='pdf',   filename='top_expressed_genes_per_cluster', title='Top expressed genes per cluster')
            )
          )
        )
      )
    ) %>% DT::formatStyle('% of total expression', textAlign='right')
  })

  output$top.expressed.genes.by.cluster.table.missing <- renderText({ 'Data not available.' })

  observeEvent(input$top.expressed.genes.by.cluster.info, {
    showModal(
      modalDialog(
        title='Top expressed genes per cluster', easyClose=TRUE, footer=NULL,
        p('Table of top 100 most expressed genes in each cluster. For example, if gene XY contributes with 5% to the total expression, that means 5% of all transcripts found in all cells of this cluster come from that respective gene. These lists can help to identify/verify the dominant cell types.')
      )
    )
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Marker genes.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$samples$names
  ## - sample_data()$samples$count
  ## - sample_data()$marker.genes$by.sample (optional)
  ## - sample_data()$clusters$names
  ## - sample_data()$clusters$count
  ## - sample_data()$marker.genes$by.cluster (optional)
  ##----------------------------------------------------------------------------##
  # by sample
  output$marker.genes.by.sample.UI <- renderUI(
    if ( sample_data()$samples$count > 1 & !is.null(sample_data()$marker.genes$by.sample) ) {
      fluidRow(
        column(12,
          selectInput('marker.genes.by.sample.input', label=NULL, choices=unname(sample_data()$samples$names)),
          DT::dataTableOutput('marker.genes.by.sample.table.present')
        )
      )
    } else {
      textOutput('marker.genes.by.sample.table.missing')
    }
  )

  output$marker.genes.by.sample.table.present <- DT::renderDataTable(server=FALSE, {
    req(input$marker.genes.by.sample.input)
    if ('on_cell_surface' %in% colnames(sample_data()$marker.genes$by.sample)) {
      table <- sample_data()$marker.genes$by.sample[ which(sample_data()$marker.genes$by.sample$sample == input$marker.genes.by.sample.input) , c(2,4,5,6,7,8) ]
      colnames(table) <- c('Gene', 'avg. logFC', '% cells in this sample', '% cells in other samples', 'adj. p-value', 'present on cell surface')
      table$'avg. logFC' <- round(table$'avg. logFC', digits=3)
      table$'% cells in this sample' <- formattable::percent(table$'% cells in this sample')
      table$'% cells in other samples' <- formattable::percent(table$'% cells in other samples')
      table <- table %>%
      mutate('adj. p-value'=formatC(table$'adj. p-value', format='e', digits=3)) %>%
      formattable::formattable(list(
        'avg. logFC' = formattable::color_tile('white', 'orange'),
        '% cells in this sample' = formattable::color_bar('pink'),
        '% cells in other samples' = formattable::color_bar('pink'),
        'present on cell surface' = formattable::formatter('span', style=x~style(color=ifelse(x, 'green', 'red')))
      ))
    } else if ( tolower(sample_data()$parameters$organism) %in% c('hg','mm') ) {
      if ( !exists('genes.surface') ) {
        genes.surface <- read.table(paste0('resources/genes_surface_', tolower(sample_data()$parameters$organism), '.txt'), sep='\t', header=FALSE, stringsAsFactors=FALSE)[,1]
      }
      table <- sample_data()$marker.genes$by.sample[ which(sample_data()$marker.genes$by.sample$sample == input$marker.genes.by.sample.input) , c(2,4,5,6,7) ]
      table$surface <- table$gene %in% genes.surface
      colnames(table) <- c('Gene', 'avg. logFC', '% cells in this sample', '% cells in other samples', 'adj. p-value', 'present on cell surface')
      table$'avg. logFC' <- round(table$'avg. logFC', digits=3)
      table$'% cells in this sample' <- formattable::percent(table$'% cells in this sample')
      table$'% cells in other samples' <- formattable::percent(table$'% cells in other samples')
      table <- table %>%
      mutate('adj. p-value'=formatC(table$'adj. p-value', format='e', digits=3)) %>%
      formattable::formattable(list(
        'avg. logFC' = formattable::color_tile('white', 'orange'),
        '% cells in this sample' = formattable::color_bar('pink'),
        '% cells in other samples' = formattable::color_bar('pink'),
        'present on cell surface' = formattable::formatter('span', style=x~style(color=ifelse(x, 'green', 'red')))
      ))
    } else {
      table <- sample_data()$marker.genes$by.sample[ which(sample_data()$marker.genes$by.sample$sample == input$marker.genes.by.sample.input) , c(2,4,5,6,7) ]
      colnames(table) <- c('Gene', 'avg. logFC', '% cells in this sample', '% cells in other samples', 'adj. p-value')
      table$'avg. logFC' <- round(table$'avg. logFC', digits=3)
      table$'% cells in this sample' <- formattable::percent(table$'% cells in this sample')
      table$'% cells in other samples' <- formattable::percent(table$'% cells in other samples')
      table <- table %>%
      mutate('adj. p-value'=formatC(table$'adj. p-value', format='e', digits=3)) %>%
      formattable::formattable(list(
        'avg. logFC' = formattable::color_tile('white', 'orange'),
        '% cells in this sample' = formattable::color_bar('pink'),
        '% cells in other samples' = formattable::color_bar('pink')
      ))
    }
    formattable::as.datatable(table,
      filter             = 'top',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      rownames           = FALSE,
      extensions         = c('Buttons'),
      class              = 'cell-border stripe',
      options            = list(
        dom        = 'Bfrtip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons    = list(
          'colvis', 
          list(
            extend  = 'collection',
            text    = 'Download',
            buttons = list(
              list(extend='csv',   filename='marker_genes_by_sample', title='Marker genes by sample'),
              list(extend='excel', filename='marker_genes_by_sample', title='Marker genes by sample'),
              list(extend='pdf',   filename='marker_genes_by_sample', title='Marker genes by sample')
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(columns=c('avg. logFC','% cells in this sample', '% cells in other samples', 'adj. p-value'), textAlign='right')
  })

  output$marker.genes.by.sample.table.missing <- renderText({ 'Only 1 sample in this data set or data not available.' })

  observeEvent(input$marker.genes.by.sample.info, {
    showModal(
      modalDialog(
        title='Marker genes per sample', easyClose=TRUE, footer=NULL,
        p('Shown here are the marker genes identified for each sample - resembling bulk RNA-seq. These genes should help to identify the cell type in this sample or find new markers to purify it. In this analysis, each sample is compared to all other samples combined. Only genes with a positive average log-fold change of at least 0.25 are reported - meaning only over-expressed genes are shown. Also, marker genes must be expressed in at least 70% of the cells of the respective sample. Statistical analysis is performed using a classical t-test as it has been shown to be very accurate in single cell RNA-seq. Finally, if data is available, the last column reports for each gene if it is associated with gene ontology term GO:0009986 which is an indicator that the respective gene is present on the cell surface (which could make it more interesting to purify a given population).')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # by cluster
  output$marker.genes.by.cluster.UI <- renderUI(
    if ( sample_data()$clusters$count > 1 & !is.null(sample_data()$marker.genes$by.cluster) ) {
      fluidRow(
        column(12,
          selectInput('marker.genes.by.cluster.input', label=NULL, choices=sample_data()$clusters$names),
          DT::dataTableOutput('marker.genes.by.cluster.table.present')
        )
      )
    } else {
      textOutput('marker.genes.by.cluster.table.missing')
    }
  )

  output$marker.genes.by.cluster.table.present <- DT::renderDataTable(server=FALSE, {
    req(input$marker.genes.by.cluster.input)
    if ('on_cell_surface' %in% colnames(sample_data()$marker.genes$by.cluster)) {
      table <- sample_data()$marker.genes$by.cluster[ which(sample_data()$marker.genes$by.cluster$cluster == input$marker.genes.by.cluster.input) , c(2,4,5,6,7,8) ]
      colnames(table) <- c('Gene', 'avg. logFC', '% cells in this cluster', '% cells in other clusters', 'adj. p-value', 'present on cell surface')
      table$'avg. logFC' <- round(table$'avg. logFC', digits=3)
      table$'% cells in this cluster' <- formattable::percent(table$'% cells in this cluster')
      table$'% cells in other clusters' <- formattable::percent(table$'% cells in other clusters')
      table <- table %>% mutate('adj. p-value'=formatC(table$'adj. p-value', format='e', digits=3)) %>%
      formattable::formattable(list(
        'avg. logFC' = formattable::color_tile('white', 'orange'),
        '% cells in this cluster' = formattable::color_bar('pink'),
        '% cells in other clusters' = formattable::color_bar('pink'),
        'present on cell surface' = formattable::formatter('span', style=x~style(color=ifelse(x, 'green', 'red')))
      ))
    } else if ( tolower(sample_data()$parameters$organism) %in% c('hg','mm') ) {
      if ( !exists('genes.surface') ) {
        genes.surface <- read.table(paste0('resources/genes_surface_', tolower(sample_data()$parameters$organism), '.txt'), sep='\t', header=FALSE, stringsAsFactors=FALSE)[,1]
      }
      table <- sample_data()$marker.genes$by.cluster[ which(sample_data()$marker.genes$by.cluster$cluster == input$marker.genes.by.cluster.input) , c(2,4,5,6,7) ]
      table$surface <- table$gene %in% genes.surface
      colnames(table) <- c('Gene', 'avg. logFC', '% cells in this cluster', '% cells in other clusters', 'adj. p-value', 'present on cell surface')
      table$'avg. logFC' <- round(table$'avg. logFC', digits=3)
      table$'% cells in this cluster' <- formattable::percent(table$'% cells in this cluster')
      table$'% cells in other clusters' <- formattable::percent(table$'% cells in other clusters')
      table <- table %>%
      mutate('adj. p-value'=formatC(table$'adj. p-value', format='e', digits=3)) %>%
      formattable::formattable(list(
        'avg. logFC' = formattable::color_tile('white', 'orange'),
        '% cells in this cluster' = formattable::color_bar('pink'),
        '% cells in other clusters' = formattable::color_bar('pink'),
        'present on cell surface' = formattable::formatter('span', style=x~style(color=ifelse(x, 'green', 'red')))
      ))
    } else {
      table <- sample_data()$marker.genes$by.cluster[ which(sample_data()$marker.genes$by.cluster$cluster == input$marker.genes.by.cluster.input) , c(2,4,5,6,7) ]
      colnames(table) <- c('Gene', 'avg. logFC', '% cells in this cluster', '% cells in other clusters', 'adj. p-value')
      table$'avg. logFC' <- round(table$'avg. logFC', digits=3)
      table$'% cells in this cluster' <- formattable::percent(table$'% cells in this cluster')
      table$'% cells in other clusters' <- formattable::percent(table$'% cells in other clusters')
      table <- table %>% mutate('adj. p-value'=formatC(table$'adj. p-value', format='e', digits=3)) %>%
      formattable::formattable(list(
        'avg. logFC' = formattable::color_tile('white', 'orange'),
        '% cells in this cluster' = formattable::color_bar('pink'),
        '% cells in other clusters' = formattable::color_bar('pink')
      ))
    }
    formattable::as.datatable(table,
      filter             = 'top',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      rownames           = FALSE,
      extensions         = c('Buttons'),
      class              = 'cell-border stripe',
      options            = list(
        dom        = 'Bfrtip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons    = list(
          'colvis', 
          list(
            extend  = 'collection',
            text    = 'Download',
            buttons = list(
              list(extend='csv',   filename='marker_genes_by_cluster', title='Marker genes by cluster'),
              list(extend='excel', filename='marker_genes_by_cluster', title='Marker genes by cluster'),
              list(extend='pdf',   filename='marker_genes_by_cluster', title='Marker genes by cluster')
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(columns=c('avg. logFC','% cells in this cluster', '% cells in other clusters', 'adj. p-value'), textAlign='right')
  })

  output$marker.genes.by.cluster.table.missing <- renderText({ 'Only 1 cluster in this data set or data not available.' })

  observeEvent(input$marker.genes.by.cluster.info, {
    showModal(
      modalDialog(
        title='Marker genes per cluster', easyClose=TRUE, footer=NULL,
        p('Shown here are the marker genes identified for each cluster. These genes should help to identify the cell type in this cluster or find new markers to purify it. In this analysis, each cluster is compared to all other clusters combined. Only genes with a positive average log-fold change of at least 0.25 are reported - meaning only over-expressed genes are shown. Also, marker genes must be expressed in at least 70% of the cells of the respective cluster. Statistical analysis is performed using a classical t-test as it has been shown to be very accurate in single cell RNA-seq. Finally, if data is available, the last column reports for each gene if it is associated with gene ontology term GO:0009986 which is an indicator that the respective gene is present on the cell surface (which could make it more interesting to purify a given population).')
      )
    )
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Enriched pathways
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$parameters$enrichr.dbs
  ## - sample_data()$samples$count
  ## - sample_data()$samples$names
  ## - sample_data()$clusters$count
  ## - sample_data()$clusters$names
  ## - sample_data()$marker.genes$by.sample.annotation (optional)
  ## - sample_data()$marker.genes$by.cluster.annotation (optional)
  ##----------------------------------------------------------------------------##
  # by sample
  output$enriched.pathways.by.sample.UI <- renderUI(
    if ( sample_data()$samples$count > 1 & !is.null(sample_data()$marker.genes$by.sample.annotation) ) {    
      tagList(
        fluidRow(
          column(4,
            selectInput('enriched.pathways.select.sample', label=NULL,
              choices = unname(sample_data()$samples$names))
          ),
          column(8,
            selectInput('enriched.pathways.select.db.for.sample', label=NULL,
              choices = sample_data()$parameters$enrichr.dbs)
          )
        ),
        fluidRow(
          column(12,
            DT::dataTableOutput('enriched.pathways.by.sample.table.present')
          )
        )
      )
    } else {
      textOutput('enriched.pathways.by.sample.table.missing')
    }
  )

  output$enriched.pathways.by.sample.table.present <- DT::renderDataTable(server=FALSE, {
    req(input$enriched.pathways.select.sample)
    req(input$enriched.pathways.select.db.for.sample)
    sample_data()$marker.genes$by.sample.annotation[[ input$enriched.pathways.select.sample ]][[ input$enriched.pathways.select.db.for.sample ]][ , c(1,2,3,4,8,9) ] %>%
    mutate(P.value=formatC(P.value, format='e', digits=3)) %>%
    mutate(Adjusted.P.value=formatC(Adjusted.P.value, format='e', digits=3)) %>%
    mutate(Combined.Score=formatC(Combined.Score, format='f', digits=2)) %>%
    formattable::formattable(list('Combined.Score' = formattable::color_bar('pink'))) %>%
    formattable::as.datatable(
      filter             = 'top',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      colnames           = c('Term', 'Overlap', 'p-value', 'adj. p-value', 'combined score', 'Genes'),
      rownames           = FALSE,
      extensions         = c('Buttons'),
      class              = 'cell-border stripe',
      options            = list(
        columnDefs = list(list(visible=FALSE, targets=c(2,5))),
        scrollX    = TRUE,
        dom        = 'Bfrtip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons    = list(
          'colvis', 
          list(
            extend  = 'collection',
            text    = 'Download',
            buttons = list(
              list(extend='csv',   filename='enriched_pathways_by_sample', title='Enriched pathways by sample'),
              list(extend='excel', filename='enriched_pathways_by_sample', title='Enriched pathways by sample'),
              list(extend='pdf',   filename='enriched_pathways_by_sample', title='Enriched pathways by sample')
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(columns=c('Combined.Score'), textAlign='right')
  })

  output$enriched.pathways.by.sample.table.missing <- renderText({ 'Only 1 sample in this data set or data not available.' })

  observeEvent(input$enriched.pathways.by.sample.info, {
    showModal(
      modalDialog(
        title='Enriched pathways by sample', easyClose=TRUE, footer=NULL,
        p('Using all marker genes identified for a respective sample, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the "Column visibility" button. For the details on the combined score is calculated, please refer to the Enrichr website and publication: http://amp.pharm.mssm.edu/Enrichr/.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # by cluster
  output$enriched.pathways.by.cluster.UI <- renderUI(
    if ( sample_data()$clusters$count > 1 & !is.null(sample_data()$marker.genes$by.cluster.annotation) ) {    
      tagList(
        fluidRow(
          column(4,
            selectInput('enriched.pathways.select.cluster', label=NULL,
              choices = sample_data()$clusters$names)
          ),
          column(8,
            selectInput('enriched.pathways.select.db.for.cluster', label=NULL,
              choices = sample_data()$parameters$enrichr.dbs)
          )
        ),
        fluidRow(
          column(12,
            DT::dataTableOutput('enriched.pathways.by.cluster.table.present')
          )
        )
      )
    } else {
      textOutput('enriched.pathways.by.cluster.table.missing')
    }
  )

  output$enriched.pathways.by.cluster.table.present <- DT::renderDataTable(server=FALSE, {
    req(input$enriched.pathways.select.cluster)
    req(input$enriched.pathways.select.db.for.cluster)
    sample_data()$marker.genes$by.cluster.annotation[[ input$enriched.pathways.select.cluster ]][[ input$enriched.pathways.select.db.for.cluster ]][ , c(1,2,3,4,8,9) ] %>%
    mutate(P.value=formatC(P.value, format='e', digits=3)) %>%
    mutate(Combined.Score=formatC(Combined.Score, format='f', digits=2)) %>%
    mutate(Adjusted.P.value=formatC(Adjusted.P.value, format='e', digits=3)) %>%
    formattable::formattable(list('Combined.Score' = formattable::color_bar('pink'))) %>%
    formattable::as.datatable(
      filter             = 'top',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      colnames           = c('Term', 'Overlap', 'p-value', 'adj. p-value', 'combined score', 'Genes'),
      rownames           = FALSE,
      extensions         = c('Buttons'),
      class              = 'cell-border stripe',
      options            = list(
        columnDefs = list(list(visible=FALSE, targets=c(2,5))),
        scrollX    = TRUE,
        dom        = 'Bfrtip', 
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 15,
        buttons    = list(
          'colvis',
          list(
            extend  = 'collection',
            text    = 'Download',
            buttons = list(
              list(extend='csv',   filename='enriched_pathways_by_cluster', title='Enriched pathways by cluster'),
              list(extend='excel', filename='enriched_pathways_by_cluster', title='Enriched pathways by cluster'),
              list(extend='pdf',   filename='enriched_pathways_by_cluster', title='Enriched pathways by cluster')
            )
          )
        )
      )
    ) %>% 
    DT::formatStyle(columns=c('Combined.Score'), textAlign='right')
  })

  output$enriched.pathways.by.cluster.table.missing <- renderText({ 'Only 1 cluster in this data set or data not available.' })

  observeEvent(input$enriched.pathways.by.cluster.info, {
    showModal(
      modalDialog(
        title='Enriched pathways by cluster', easyClose=TRUE, footer=NULL,
        p('Using all marker genes identified for a respective cluster, gene list enrichment analysis is performed using the Enrichr API, including gene ontology terms, KEGG and Wiki Pathways, BioCarta and many others. Terms are sorted based on the combined score. By default, the genes that overlap between the marker gene list and a term are not shown (for better visibility) but the column can be added using the "Column visibility" button. For the details on the combined score is calculated, please refer to the Enrichr website and publication: http://amp.pharm.mssm.edu/Enrichr/')
      )
    )
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Gene expression.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$projections
  ## - sample_data()$samples$names
  ## - sample_data()$clusters$names
  ## - sample_data()$expression
  ## - sample_data()$cells$sample
  ## - sample_data()$cells$cluster
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent.mt
  ## - sample_data()$cells$percent.ribo
  ##----------------------------------------------------------------------------##
 
  # reactive data
  genesToPlot <- reactive({
    genesToPlot <- list()
    if ( is.null(input$expression.genes.input) ) {
      genesToPlot$genes.to.display <- ''
    } else {
      genesToPlot$genes.to.display <- input$expression.genes.input
      genesToPlot$genes.to.display <- unlist(strsplit(genesToPlot$genes.to.display, ',| |;|\n'))
      genesToPlot$genes.to.display <- gsub(x=genesToPlot$genes.to.display, pattern=' ', replacement='', fixed=TRUE)
      genesToPlot$genes.to.display <- unique(genesToPlot$genes.to.display)
    }
    genesToPlot$genes.to.display.here    <- rownames(sample_data()$expression)[ match(tolower(genesToPlot$genes.to.display), tolower(rownames(sample_data()$expression))) ]
    genesToPlot$genes.to.display.present <- genesToPlot$genes.to.display.here[ which(!is.na(genesToPlot$genes.to.display.here)) ]
    genesToPlot$genes.to.display.missing <- genesToPlot$genes.to.display[ which(is.na(genesToPlot$genes.to.display.here)) ]
    return(genesToPlot)
  })

  # select genes to be displayed
  output$expression.genes.displayed <- renderText({
    genes.to.display.text <- paste0('<b>Showing expression for ', length(genesToPlot()$genes.to.display.present), ' gene(s):</b><br>',
                                    paste0(genesToPlot()$genes.to.display.present, collapse=', '), '<br><b>',
                                    length(genesToPlot()$genes.to.display.missing) , ' gene(s) are not in data set: </b><br>',
                                    paste0(genesToPlot()$genes.to.display.missing, collapse=', '))
    return(genes.to.display.text)
  })

  ##--------------------------------------------------------------------------##
 # UI
  output$expression.UI <- renderUI({
    tagList(
      selectInput('expression.projection.to.display', label='Projection:',
        choices  = names(sample_data()$projections)),
      textAreaInput('expression.genes.input', label='Gene(s):',
        value = '',
        placeholder = 'Insert genes here.'),
      shinyWidgets::pickerInput('expression.samples.to.display', label='Samples to display:',
        choices  = unname(sample_data()$samples$names),
        selected = unname(sample_data()$samples$names),
        options  = list('actions-box'=TRUE),
        multiple = TRUE),
      shinyWidgets::pickerInput('expression.clusters.to.display', label='Clusters to display:',
        choices  = sample_data()$clusters$names,
        selected = sample_data()$clusters$names,
        options  = list('actions-box'=TRUE),
        multiple = TRUE),
      selectInput('expression.plotting.order', label='Plotting order:',
        choices  = c('Random','Highest expression on top')),
      sliderInput('expression.projection.dot.size', label='Point size:',
        min=0, max=50, value=25, step=1),
      sliderInput('expression.projection.opacity', label='Point opacity:',
        min=0, max=1, value=1, step=0.05)
    )
  })

  output$expression.scales <- renderUI({
    projection.to.display <- if ( is.null(input$expression.projection.to.display) || is.na(input$expression.projection.to.display) ) names(sample_data()$projections)[1] else input$expression.projection.to.display
    range.x.min <- round(min(sample_data()$projections[[ projection.to.display ]][,1])*1.1)
    range.x.max <- round(max(sample_data()$projections[[ projection.to.display ]][,1])*1.1)
    range.y.min <- round(min(sample_data()$projections[[ projection.to.display ]][,2])*1.1)
    range.y.max <- round(max(sample_data()$projections[[ projection.to.display ]][,2])*1.1)
    tagList(
      sliderInput('expression.projection.scale.x.manual.range', label='X axis',
        min=range.x.min, max=range.x.max, value=c(range.x.min, range.x.max)),
      sliderInput('expression.projection.scale.y.manual.range', label='Y axis',
        min=range.y.min, max=range.y.max, value=c(range.y.min, range.y.max))
    )
  })

  ##--------------------------------------------------------------------------##

  # projection with scatterD3
  output$expression.projection <- scatterD3::renderScatterD3({
    req(input$expression.projection.to.display)
    req(input$expression.samples.to.display)
    req(input$expression.clusters.to.display)
    req(input$expression.plotting.order)
    projection.to.display <- input$expression.projection.to.display
    samples.to.display    <- input$expression.samples.to.display
    clusters.to.display   <- input$expression.clusters.to.display
    cells.to.display      <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', samples.to.display, '$', collapse='|')) == TRUE & 
                                   grepl(sample_data()$cells$cluster, pattern=paste0('^', clusters.to.display, '$', collapse='|')) == TRUE)
    to.plot               <- cbind(sample_data()$projections[[ projection.to.display ]][ cells.to.display , ],
                                   sample_data()$cells[ cells.to.display , ])
    if ( length(genesToPlot()$genes.to.display.present) == 0 ) {
      to.plot$level <- 0
    } else if ( length(genesToPlot()$genes.to.display.present) == 1 ) {
      to.plot$level <- sample_data()$expression[ which(rownames(sample_data()$expression) == genesToPlot()$genes.to.display.present) , cells.to.display ]
    } else {
      to.plot$level <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', genesToPlot()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
    }
    if ( input$expression.plotting.order == 'Random' ) {
      to.plot <- to.plot[ sample(1:length(to.plot$level), length(to.plot$level)) , ]
    } else if ( input$expression.plotting.order == 'Highest expression on top' ) {
      to.plot <- to.plot[ order(to.plot$level, decreasing=FALSE) , ]
    }
    scatterD3::scatterD3(
      x              = to.plot[ , 1 ],
      y              = to.plot[ , 2 ],
      xlab           = colnames(to.plot)[ 1 ],
      ylab           = colnames(to.plot)[ 2 ],
      xlim          = c(input$expression.projection.scale.x.manual.range[1],input$expression.projection.scale.x.manual.range[2]),
      ylim          = c(input$expression.projection.scale.y.manual.range[1],input$expression.projection.scale.y.manual.range[2]),
      point_size     = input$expression.projection.dot.size,
      col_var        = to.plot$level,
      col_lab        = 'Gene expression',
      col_continuous = TRUE,
      point_opacity  = input$expression.projection.opacity,
      transitions    = FALSE,
      legend_width   = 0,
      menu           = FALSE,
      tooltip_text   = paste0(
        '<b>Sample</b>: ', to.plot[ , 'sample'       ], '<br/>',
        '<b>Cluster</b>: ',   to.plot[ , 'cluster'      ], '<br/>',
        '<b>nUMI</b>: ',      to.plot[ , 'nUMI'         ], '<br/>',
        '<b>nGene</b>: ',     to.plot[ , 'nGene'        ], '<br/>',
        '<b>Expr. MT</b>: ',   format(to.plot[ , 'percent.mt'   ]*100, digits=1), '%<br/>',
        '<b>Expr. ribo</b>: ', format(to.plot[ , 'percent.ribo' ]*100, digits=1), '%<br/>'))
  })

  observeEvent(input$expression.projection.info, {
    showModal(
      modalDialog(
        title='Dimensional reduction', easyClose=TRUE, footer=NULL,
        p('Interactive projection of cells into 2-dimensional space based on their expression profile.', 
          tags$ul(tags$li('Both tSNE and UMAP are frequently used algorithms for dimensional reduction in single cell transcriptomics. While they generally allow to make similar conclusions, some differences exist between the two (please refer to Google).'),
                  tags$li('Cell color reflects the log-normalised expression of entered genes. If more than 1 gene is entered, the color reflects the average expression of all genes. Genes must be in separate lines or separated by a space, comma, or semicolon. Reported below the projection are the genes that are present and absent in this data set. Absent genes could either have been annotated with a different name or were not expressed in any of the cells. Matching of gene names is case-insensitive, that means Myc/MYC/myc are treated equally.'),
                  tags$li('Samples and clusters can be removed from the plot individually to highlight a contrast of interest.'),
                  tags$li('Cells can be plotted either randomly (which a more unbiased image) or in the order of expression (with highest expression plotted last), sometimes resulting in a more appealing figure.'),
                  tags$li('By default, the dot size is set to 15 without any transparency but both these attributes can be changed using the sliders on the left.'),
                  tags$li('The last 2 slider elements on the left can be used to resize the projection axes. This can be particularly useful when a projection contains a population of cell that is very far away from the rest and therefore creates a big empty space (which is not uncommon for UMAPs).')
          ),
          'The plot is interactive (drag and zoom) but depending on the computer of the user and the number of cells displayed it can become very slow.'
        )
      )
    )
  })

  observeEvent(input$expression.projection.export, {
    library('ggplot2')

    projection.to.display <- input$expression.projection.to.display
    samples.to.display    <- input$expression.samples.to.display
    clusters.to.display   <- input$expression.clusters.to.display
    cells.to.display      <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', samples.to.display, '$', collapse='|')) == TRUE & grepl(sample_data()$cells$cluster, pattern=paste0('^', clusters.to.display, '$', collapse='|')) == TRUE)
    to.plot               <- cbind(sample_data()$projections[[ projection.to.display ]][ cells.to.display , ],
                                   sample_data()$cells[ cells.to.display , ])

    xlim <- c(input$expression.projection.scale.x.manual.range[1], input$expression.projection.scale.x.manual.range[2])
    ylim <- c(input$expression.projection.scale.y.manual.range[1], input$expression.projection.scale.y.manual.range[2])

    if ( length(genesToPlot()$genes.to.display.present) == 0 ) {
      to.plot$level <- 0
      out.filename <- paste0(plot.export.path, 'scBrowser_', sample_data()$parameters$project.name, '_gene_expression_none')
    } else if ( length(genesToPlot()$genes.to.display.present) == 1 ) {
      to.plot$level <- sample_data()$expression[ which(rownames(sample_data()$expression) == genesToPlot()$genes.to.display.present) , cells.to.display ]
      out.filename <- paste0(plot.export.path, 'scBrowser_', sample_data()$parameters$project.name, '_gene_expression_', genesToPlot()$genes.to.display.present, '_', input$expression.projection.to.display)
    } else {
      to.plot$level <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', genesToPlot()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
      out.filename <- paste0(plot.export.path, 'scBrowser_', gsub(sample_data()$parameters$project.name, pattern=' ', replacement='_'), '_gene_expression_', genesToPlot()$genes.to.display.present[1], '_and_others_', input$expression.projection.to.display)
    }

    if ( input$expression.plotting.order == 'Random' ) {
      to.plot <- to.plot[ sample(1:length(to.plot$level), length(to.plot$level)) , ]
      out.filename <- paste0(out.filename, '_random_order.pdf')
    } else if ( input$expression.plotting.order == 'Highest expression on top' ) {
      to.plot <- to.plot[ order(to.plot$level, decreasing=FALSE) , ]
      out.filename <- paste0(out.filename, '_highest_expression_on_top.pdf')
    }

    p <- ggplot(to.plot, aes_q(x=as.name(colnames(to.plot)[1]), y=as.name(colnames(to.plot)[2]), colour=as.name('level'))) +
         geom_point() +
         viridis::scale_colour_viridis(name='Log-normalised\nexpression', guide=guide_colorbar(frame.colour='black', ticks.colour='black')) +
         lims(x=xlim, y=ylim) +
         theme_bw()

    pdf(NULL)
    ggsave(out.filename, p, height=8, width=11)

    if (file.exists(out.filename)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = 'Success!',
        text = paste0('Plot saved successfully as: ', out.filename),
        type = 'success'
      )
    } else {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = 'Error!',
        text = "Sorry, it seems something went wrong...",
        type = 'error'
      )
    }
  })

  ##--------------------------------------------------------------------------##

  # box plot by sample
  output$expression.by.sample <- plotly::renderPlotly({
    cells.to.display <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', sample_data()$samples$names, '$', collapse='|')) == TRUE )
    if ( length(genesToPlot()$genes.to.display.present) == 0 ) {
      expression.levels <- 0
    } else if ( length(genesToPlot()$genes.to.display.present) == 1 ) {
      expression.levels <- sample_data()$expression[ which(rownames(sample_data()$expression) == genesToPlot()$genes.to.display.present) , cells.to.display ]
    } else {
      expression.levels <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', genesToPlot()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
    }
    to_display <- data.frame('sample' = sample_data()$cells[ , 'sample' ], 'expression' = expression.levels)
    plotly::plot_ly(to_display,
      x          = ~sample,
      y          = ~expression,
      type       = 'box',
      color      = ~sample,
      colors     = colors[ 1:sample_data()$samples$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Expression level', hoverformat='.4f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$expression.by.sample.info, {
    showModal(
      modalDialog(
        title='Expression levels by sample', easyClose=TRUE, footer=NULL,
        p('Log-normalised expression of genes inserted above by sample. If more than 1 gene was given, this reflects the average across all cells of each sample.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot by cluster
  output$expression.by.cluster <- plotly::renderPlotly({
    cells.to.display <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', sample_data()$samples$names, '$', collapse='|')) == TRUE )
    if ( length(genesToPlot()$genes.to.display.present) == 0 ) {
      expression.levels <- 0
    } else if ( length(genesToPlot()$genes.to.display.present) == 1 ) {
      expression.levels <- sample_data()$expression[ which(rownames(sample_data()$expression) == genesToPlot()$genes.to.display.present) , cells.to.display ]
    } else {
      expression.levels <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', genesToPlot()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
    }
    to_display <- data.frame('cluster' = sample_data()$cells[ , 'cluster' ], 'expression' = expression.levels)
    plotly::plot_ly(to_display,
      x          = ~cluster,
      y          = ~expression,
      type       = 'box',
      color      = ~cluster,
      colors     = colors[ 1:sample_data()$clusters$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size=5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Expression level', hoverformat='.4f'),
      dragmode  =  'select',
      hovermode = 'compare')
  })

  observeEvent(input$expression.by.cluster.info, {
    showModal(
      modalDialog(
        title='Expression levels by cluster', easyClose=TRUE, footer=NULL,
        p('Log-normalised expression of genes inserted above by cluster. If more than 1 gene was given, this reflects the average across all cells of each cluster.')
      )
    )
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Gene set expression.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - sample_data()$parameters$organism
  ## - sample_data()$expression
  ## - sample_data()$projections
  ## - sample_data()$samples$names
  ## - sample_data()$clusters$names
  ## - sample_data()$cells$sample
  ## - sample_data()$cells$cluster
  ## - sample_data()$cells$nUMI
  ## - sample_data()$cells$nGene
  ## - sample_data()$cells$percent.mt
  ## - sample_data()$cells$percent.ribo
  ##----------------------------------------------------------------------------##
  # reactive data
  geneSets <- reactive({
    if ( sample_data()$parameters$organism == 'mm' ) {
      msigdbr::msigdbr(species='Mus musculus')
    } else if ( sample_data()$parameters$organism == 'hg' ) {
      msigdbr::msigdbr(species='Homo sapiens')
    } else {
      msigdbr::msigdbr(species='Mus musculus')
    }
  })

  # reactive data
  geneSetData <- reactive({
    geneSetData <- list()
    if ( is.null(input$geneSetExpression.select.geneSet) || is.na(input$geneSetExpression.select.geneSet) || input$geneSetExpression.select.geneSet == '-' ) {
      geneSetData$genes.to.display.present <- NULL
    } else {
      geneSetData$genes.to.display         <- geneSets()[ which(geneSets()$gs_name == input$geneSetExpression.select.geneSet) , 'gene_symbol' ]$gene_symbol
      geneSetData$genes.to.display         <- unique(geneSetData$genes.to.display)
      # print(head(tolower(rownames(sample_data()$expression))))
      geneSetData$genes.to.display.here    <- rownames(sample_data()$expression)[ match(tolower(geneSetData$genes.to.display), tolower(rownames(sample_data()$expression))) ]
      geneSetData$genes.to.display.present <- geneSetData$genes.to.display.here[ which(!is.na(geneSetData$genes.to.display.here)) ]
      # geneSetData$genes.to.display.present <- geneSetData$genes.to.display[ which(!is.na(geneSetData$genes.to.display.here)) ]
      geneSetData$genes.to.display.missing <- geneSetData$genes.to.display[ which(is.na(geneSetData$genes.to.display.here)) ]
    }
    return(geneSetData)
  })

  # show which genes are in the data set and which aren't
  output$geneSetExpression.genes.displayed <- renderText({
    genes.to.display.text <- paste0(
      '<br><b>Total unique genes in gene set:</b><br>',
      length(geneSetData()$genes.to.display),
      '<br><b>Showing expression for ', length(geneSetData()$genes.to.display.present), ' genes:</b><br>',
      paste0(geneSetData()$genes.to.display.present, collapse=', '),
      '<br><b>', length(geneSetData()$genes.to.display.missing) , ' gene(s) are not in data set: </b><br>',
      paste0(geneSetData()$genes.to.display.missing, collapse=', '))
    if ( sample_data()$parameters$organism != 'mm' & sample_data()$parameters$organism != 'hg' ) {
      return(paste0('<br><b><font color="red">Note:</b> Gene sets are available for human and mouse only. Organism for loaded samples is either not set or none of the two. Mouse gene sets are loaded and can be used.</font><br>',
                    genes.to.display.text))
    } else {
      return(genes.to.display.text)    
    }
  })

  ##--------------------------------------------------------------------------##

  # UI
  output$geneSetExpression.UI <- renderUI({
    tagList(
      selectInput('geneSetExpression.projection.to.display', label='Projection:',
        choices  = names(sample_data()$projections)),
      selectInput('geneSetExpression.select.geneSet', label='Gene set:',
        choices  = c('-', unique(geneSets()$gs_name)),
        selected = '-'),
      shinyWidgets::pickerInput('geneSetExpression.samples.to.display', label='Samples to display:',
        choices  = unname(sample_data()$samples$names),
        selected = unname(sample_data()$samples$names),
        options  = list('actions-box'=TRUE),
        multiple = TRUE),
      shinyWidgets::pickerInput('geneSetExpression.clusters.to.display', label='Clusters to display:',
        choices  = sample_data()$clusters$names,
        selected = sample_data()$clusters$names,
        options  = list('actions-box'=TRUE),
        multiple = TRUE),
      selectInput('geneSetExpression.plotting.order', label='Plotting order:',
        choices  = c('Random','Highest expression on top')),
      sliderInput('geneSetExpression.projection.dot.size', label='Point size:',
        min=0, max=50, value=25, step=1),
      sliderInput('geneSetExpression.projection.opacity', label='Point opacity:',
        min=0, max=1, value=1, step=0.05)
    )
  })

  output$geneSetExpression.scales <- renderUI({
    projection.to.display <- if ( is.null(input$geneSetExpression.projection.to.display) || is.na(input$geneSetExpression.projection.to.display) ) names(sample_data()$projections)[1] else input$geneSetExpression.projection.to.display
    range.x.min <- round(min(sample_data()$projections[[ projection.to.display ]][,1])*1.1)
    range.x.max <- round(max(sample_data()$projections[[ projection.to.display ]][,1])*1.1)
    range.y.min <- round(min(sample_data()$projections[[ projection.to.display ]][,2])*1.1)
    range.y.max <- round(max(sample_data()$projections[[ projection.to.display ]][,2])*1.1)
    tagList(
      sliderInput('geneSetExpression.projection.scale.x.manual.range', label='X axis',
        min=range.x.min, max=range.x.max, value=c(range.x.min, range.x.max)),
      sliderInput('geneSetExpression.projection.scale.y.manual.range', label='Y axis',
        min=range.y.min, max=range.y.max, value=c(range.y.min, range.y.max))
    )
  })

  ##--------------------------------------------------------------------------##

  # projection with scatterD3
  output$geneSetExpression.projection <- scatterD3::renderScatterD3({
    req(input$geneSetExpression.projection.to.display)
    req(input$geneSetExpression.samples.to.display)
    req(input$geneSetExpression.clusters.to.display)
    req(input$geneSetExpression.plotting.order)
    projection.to.display <- input$geneSetExpression.projection.to.display
    samples.to.display    <- input$geneSetExpression.samples.to.display
    clusters.to.display   <- input$geneSetExpression.clusters.to.display
    cells.to.display      <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', samples.to.display, '$', collapse='|')) == TRUE & grepl(sample_data()$cells$cluster, pattern=paste0('^', clusters.to.display, '$', collapse='|')) == TRUE)
    to.plot               <- cbind(sample_data()$projections[[ projection.to.display ]][ cells.to.display , ],
                                   sample_data()$cells[ cells.to.display , ])
    if ( length(geneSetData()$genes.to.display.present) == 0 ) {
      to.plot$level <- 0
    } else if ( length(geneSetData()$genes.to.display.present) == 1 ) {
      to.plot$level <- sample_data()$expression[ which(rownames(sample_data()$expression) == geneSetData()$genes.to.display.present) , cells.to.display ]
    } else {
      to.plot$level <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', geneSetData()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
    }
    if ( input$geneSetExpression.plotting.order == 'Random' ) {
      to.plot <- to.plot[ sample(1:length(to.plot$level), length(to.plot$level)) , ]
    } else if ( input$geneSetExpression.plotting.order == 'Highest expression on top' ) {
      to.plot <- to.plot[ order(to.plot$level, decreasing=FALSE) , ]
    }
    scatterD3::scatterD3(
      x              = to.plot[ , 1 ],
      y              = to.plot[ , 2 ],
      xlab           = colnames(to.plot)[ 1 ],
      ylab           = colnames(to.plot)[ 2 ],
      xlim           = c(input$geneSetExpression.projection.scale.x.manual.range[1],input$geneSetExpression.projection.scale.x.manual.range[2]),
      ylim           = c(input$geneSetExpression.projection.scale.y.manual.range[1],input$geneSetExpression.projection.scale.y.manual.range[2]),
      point_size     = input$geneSetExpression.projection.dot.size,
      col_var        = to.plot$level,
      col_lab        = 'Gene expression',
      col_continuous = TRUE,
      point_opacity  = input$geneSetExpression.projection.opacity,
      transitions    = FALSE,
      legend_width   = 0,
      menu           = FALSE,
      tooltip_text   = paste0(
        '<b>Sample</b>: ', to.plot[ , 'sample'       ], '<br/>',
        '<b>Cluster</b>: ',   to.plot[ , 'cluster'      ], '<br/>',
        '<b>nUMI</b>: ',      to.plot[ , 'nUMI'         ], '<br/>',
        '<b>nGene</b>: ',     to.plot[ , 'nGene'        ], '<br/>',
        '<b>Expr. MT</b>: ',   format(to.plot[ , 'percent.mt'   ]*100, digits=1), '%<br/>',
        '<b>Expr. ribo</b>: ', format(to.plot[ , 'percent.ribo' ]*100, digits=1), '%<br/>'))
  })

  observeEvent(input$geneSetExpression.projection.info, {
    showModal(
      modalDialog(
        title='Dimensional reduction', easyClose=TRUE, footer=NULL,
        p('Interactive projection of cells into 2-dimensional space based on their expression profile.', 
          tags$ul(tags$li('Both tSNE and UMAP are frequently used algorithms for dimensional reduction in single cell transcriptomics. While they generally allow to make similar conclusions, some differences exist between the two (please refer to Google).'),
                  tags$li('For human and murine data sets, all organism-specific gene sets from the MSigDB can be selected. If the experiment was performed in another organism, the murine gene sets will be available.'),
                  tags$li('Cell color reflects the average log-normalised expression of the genes in the selected gene set. Reported below the projection are the genes that are present and absent in this data set. Absent genes could either have been annotated with a different name or were not expressed in any of the cells. Matching of gene names is case-insensitive, that means Myc/MYC/myc are treated equally.'),
                  tags$li('Samples and clusters can be removed from the plot individually to highlight a contrast of interest.'),
                  tags$li('Cells can be plotted either randomly (which a more unbiased image) or in the order of expression (with highest expression plotted last), sometimes resulting in a more appealing figure.'),
                  tags$li('By default, the dot size is set to 15 without any transparency but both these attributes can be changed using the sliders on the left.'),
                  tags$li('The last 2 slider elements on the left can be used to resize the projection axes. This can be particularly useful when a projection contains a population of cell that is very far away from the rest and therefore creates a big empty space (which is not uncommon for UMAPs).')
          ),
          'The plot is interactive (drag and zoom) but depending on the computer of the user and the number of cells displayed it can become very slow.'
        )
      )
    )
  })

  observeEvent(input$geneSetExpression.projection.export, {
    library('ggplot2')

    projection.to.display <- input$geneSetExpression.projection.to.display
    samples.to.display    <- input$geneSetExpression.samples.to.display
    clusters.to.display   <- input$geneSetExpression.clusters.to.display
    cells.to.display      <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', samples.to.display, '$', collapse='|')) == TRUE & grepl(sample_data()$cells$cluster, pattern=paste0('^', clusters.to.display, '$', collapse='|')) == TRUE)
    to.plot               <- cbind(sample_data()$projections[[ projection.to.display ]][ cells.to.display , ],
                                   sample_data()$cells[ cells.to.display , ])

    xlim <- c(input$geneSetExpression.projection.scale.x.manual.range[1], input$geneSetExpression.projection.scale.x.manual.range[2])
    ylim <- c(input$geneSetExpression.projection.scale.y.manual.range[1], input$geneSetExpression.projection.scale.y.manual.range[2])

    if ( length(geneSetData()$genes.to.display.present) == 0 ) {
      to.plot$level <- 0
    } else if ( length(geneSetData()$genes.to.display.present) == 1 ) {
      to.plot$level <- sample_data()$expression[ which(rownames(sample_data()$expression) == geneSetData()$genes.to.display.present) , cells.to.display ]
    } else {
      to.plot$level <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', geneSetData()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
    }

    out.filename <- paste0(plot.export.path, 'scBrowser_', gsub(sample_data()$parameters$project.name, pattern=' ', replacement='_'), '_gene_set_expression_', input$geneSetExpression.select.geneSet, '_', input$geneSetExpression.projection.to.display)

    if ( input$geneSetExpression.plotting.order == 'Random' ) {
      to.plot <- to.plot[ sample(1:length(to.plot$level), length(to.plot$level)) , ]
      out.filename <- paste0(out.filename, '_random_order.pdf')
    } else if ( input$geneSetExpression.plotting.order == 'Highest expression on top' ) {
      to.plot <- to.plot[ order(to.plot$level, decreasing=FALSE) , ]
      out.filename <- paste0(out.filename, '_highest_expression_on_top.pdf')
    }

    p <- ggplot(to.plot, aes_q(x=as.name(colnames(to.plot)[1]), y=as.name(colnames(to.plot)[2]), colour=as.name('level'))) +
         geom_point() +
         viridis::scale_colour_viridis(name='Average\nlog-normalised\nexpression', guide=guide_colorbar(frame.colour='black', ticks.colour='black')) +
         lims(x=xlim, y=ylim) +
         theme_bw()

    pdf(NULL)
    ggsave(out.filename, p, height=8, width=11)

    if (file.exists(out.filename)) {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = 'Success!',
        text = paste0('Plot saved successfully as: ', out.filename),
        type = 'success'
      )
    } else {
      shinyWidgets::sendSweetAlert(
        session = session,
        title = 'Error!',
        text = "Sorry, it seems something went wrong...",
        type = 'error'
      )
    }
  })

  ##--------------------------------------------------------------------------##

  # box plot by sample
  output$geneSetExpression.by.sample <- plotly::renderPlotly({
    cells.to.display <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', sample_data()$samples$names, '$', collapse='|')) == TRUE )
    if ( length(geneSetData()$genes.to.display.present) == 0 ) {
      expression.levels <- 0
    } else if ( length(geneSetData()$genes.to.display.present) == 1 ) {
      expression.levels <- sample_data()$expression[ which(rownames(sample_data()$expression) == geneSetData()$genes.to.display.present) , cells.to.display ]
    } else {
      expression.levels <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', geneSetData()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
    }
    to_display <- data.frame('sample' = sample_data()$cells[ , 'sample' ], 'expression' = expression.levels)
    plotly::plot_ly(to_display,
      x          = ~sample,
      y          = ~expression,
      type       = 'box',
      color      = ~sample,
      colors     = colors[ 1:sample_data()$samples$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size = 5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title = ''),
      yaxis     = list(title = 'Expression level', hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$geneSetExpression.by.sample.info, {
    showModal(
      modalDialog(
        title='Expression levels by sample', easyClose=TRUE, footer=NULL,
        p('Average log-normalised expression of genes in selected gene set by sample.')
      )
    )
  })

  ##--------------------------------------------------------------------------##

  # box plot by cluster
  output$geneSetExpression.by.cluster <- plotly::renderPlotly({
    cells.to.display <- which(grepl(sample_data()$cells$sample, pattern=paste0('^', sample_data()$samples$names, '$', collapse='|')) == TRUE )
    if ( length(geneSetData()$genes.to.display.present) == 0 ) {
      expression.levels <- 0
    } else if ( length(geneSetData()$genes.to.display.present) == 1 ) {
      expression.levels <- sample_data()$expression[ which(rownames(sample_data()$expression) == geneSetData()$genes.to.display.present) , cells.to.display ]
    } else {
      expression.levels <- as.vector(colMeans(as.matrix(sample_data()$expression[ which(grepl(rownames(sample_data()$expression), pattern=paste0('^', geneSetData()$genes.to.display.present, '$', collapse='|'))) , cells.to.display ])))
    }
    to_display <- data.frame('cluster' = sample_data()$cells[ , 'cluster' ], 'expression' = expression.levels)
    plotly::plot_ly(to_display,
      x          = ~cluster,
      y          = ~expression,
      type       = 'box',
      color      = ~cluster,
      colors     = colors[ 1:sample_data()$clusters$count ],
      source     = 'subset',
      showlegend = FALSE,
      hoverinfo  = 'y',
      marker     = list(size = 5)) %>%
    plotly::layout(
      title     = '',
      xaxis     = list(title=''),
      yaxis     = list(title='Expression level', hoverformat='.2f'),
      dragmode  = 'select',
      hovermode = 'compare')
  })

  observeEvent(input$geneSetExpression.by.cluster.info, {
    showModal(
      modalDialog(
        title='Expression levels by cluster', easyClose=TRUE, footer=NULL,
        p('Average log-normalised expression of genes in selected gene set by cluster.')
      )
    )
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Gene id/symbol conversion.
  ##----------------------------------------------------------------------------##
  output$gene_info <- DT::renderDataTable({
    if ( input$geneIdConversion.organism == 'mouse' ) {
      conversion.table <- read.table('resources/mm10_gene_ID_name.txt',
        sep='\t', header=TRUE, stringsAsFactors=FALSE)
    } else if ( input$geneIdConversion.organism == 'human' ) {
      conversion.table <- read.table('resources/hg38_gene_ID_name.txt',
        sep='\t', header=TRUE, stringsAsFactors=FALSE)
    }
    DT::datatable(data = conversion.table,
      filter             = 'none',
      selection          = 'multiple',
      escape             = FALSE,
      autoHideNavigation = TRUE,
      rownames           = FALSE,
      options            = list(
        scrollX    = FALSE,
        dom        = 'Bfrtip',
        lengthMenu = c(15, 30, 50, 100),
        pageLength = 50))
  })


  ##----------------------------------------------------------------------------##
  ## Panel: Sample info.
  ##----------------------------------------------------------------------------##
  ## Expected data:
  ## - see below
  ##----------------------------------------------------------------------------##
  output$sample.info.general <- renderText({
    info <- paste0(
              '<ul>',
                '<li><b>Data version:</b> ', sample_data()$data.version,
                '<li><b>Data type:</b> ', sample_data()$data.type,
                '<li><b>Project name:</b> ', sample_data()$parameters$project.name,
                '<li><b>Organism:</b> ', sample_data()$parameters$organism,
                '<li><b>Reference version:</b> ', sample_data()$parameters$reference.version,
                '<li><b>Annotation:</b> ', sample_data()$parameters$annotation,
                '<li><b>Annotation type:</b> ', sample_data()$parameters$annotation.type,
                '<li><b>Min. cells:</b> ', sample_data()$parameters$min.cells,
                '<li><b>Variables to regress:</b> ', sample_data()$parameters$vars.to.regress,
                '<li><b>Number of PCs:</b> ', sample_data()$parameters$number.PCs,
                '<li><b>tSNE perplexity:</b> ', sample_data()$parameters$tSNE.perplexity,
                '<li><b>Cluster resolution:</b> ', if ( !is.null(sample_data()$parameters$cluster.resolution) ) sample_data()$parameters$cluster.resolution else 'Not available.',
                '<li><b>enrichR databases:</b> ', paste0(sample_data()$parameters$enrichr.dbs, collapse=', '),
                '<li><b>Mitochondrial genes:</b> ', paste0(sample_data()$gene.lists$genes.mt, collapse=', '),
                '<li><b>Ribosomal genes:</b> ', paste0(sample_data()$gene.lists$genes.ribo, collapse=', '),
                '<li><b>Genes used for apoptotic score:</b> ', paste0(sample_data()$gene.lists$genes.apoptosis, collapse=', '),
                '<li><b>S phase genes:</b> ', paste0(sample_data()$gene.lists$genes.S, collapse=', '),
                '<li><b>G2M phase genes:</b> ', paste0(sample_data()$gene.lists$genes.G2M, collapse=', '),
              '</ul>'
            )
    if ( !is.null(sample_data()$parameters$info.R) ) {
      info.R.raw <- sample_data()$parameters$info.R
      info.R <- c()
      for ( i in 1:length(info.R.raw) ) {
        info.R <- paste(info.R, '<br>', info.R.raw[i])
      }
      info <- paste0(info, '<br><b>R environment and packages use in analysis:</b><br><pre>', info.R, '</pre>')
    }
    return(info)
  })


  # output$sample.info.general <- renderText({
  #   info <- paste0('<b>Data version:</b> ', sample_data()$data.version, '<br>',
  #                  '<b>Data type:</b> ', sample_data()$data.type, '<br>',
  #                  '<b>Project name:</b> ', sample_data()$parameters$project.name, '<br>',
  #                  '<b>Organism:</b> ', sample_data()$parameters$organism, '<br>',
  #                  '<b>Reference version:</b> ', sample_data()$parameters$reference.version, '<br>',
  #                  '<b>Annotation:</b> ', sample_data()$parameters$annotation, '<br>',
  #                  '<b>Annotation type:</b> ', sample_data()$parameters$annotation.type, '<br>',
  #                  '<b>Min. cells:</b> ', sample_data()$parameters$min.cells, '<br>',
  #                  '<b>Variables to regress:</b> ', sample_data()$parameters$vars.to.regress, '<br>',
  #                  '<b>Number of PCs:</b> ', sample_data()$parameters$number.PCs, '<br>',
  #                  '<b>tSNE perplexity:</b> ', sample_data()$parameters$tSNE.perplexity, '<br>',
  #                  paste0('<b>enrichR databases:</b> ', paste0(sample_data()$parameters$enrichr.dbs, collapse=', ')), '<br>',
  #                  paste0('<b>Mitochondrial genes:</b> ', paste0(sample_data()$gene.lists$genes.mt, collapse=', ')), '<br>',
  #                  paste0('<b>Ribosomal genes:</b> ', paste0(sample_data()$gene.lists$genes.ribo, collapse=', ')), '<br>',
  #                  paste0('<b>Genes used for apoptotic score:</b> ', paste0(sample_data()$gene.lists$genes.apoptosis, collapse=', ')), '<br>',
  #                  paste0('<b>S phase genes:</b> ', paste0(sample_data()$gene.lists$genes.S, collapse=', ')), '<br>',
  #                  paste0('<b>G2M phase genes:</b> ', paste0(sample_data()$gene.lists$genes.G2M, collapse=', '))
  #           )
  #   if ( !is.null(sample_data()$parameters$info.R) ) {
  #     info.R.raw <- capture.output(sample_data()$parameters$info.R)
  #     info.R <- c()
  #     for ( i in 1:length(info.R.raw) ) {
  #       info.R <- paste(info.R, '<br>', info.R.raw[i])
  #     }
  #     info <- paste0(info, '<br><b>R environment and packages use in analysis:</b><br>', info.R)
  #   }
  #   return(info)
  # })


  output$sample.info.R <- renderPrint({
    if ( !is.null(sample_data()$parameters$info.R) ) {
      capture.output(sample_data()$parameters$info.R)
    } else {
      print('Not available')
    }
  })


  ##----------------------------------------------------------------------------##
  ## Panel: About.
  ##----------------------------------------------------------------------------##
  output$about <- renderText({
    '<b>Version:</b><br>
     1.0.4 (October 2018)<br>
     <br>
     <b>Author:</b><br>
     Roman Hillje<br>
     Department of Experimental Oncology<br>
     IEO, European Institute of Oncology IRCCS, Milan<br>
     <br>
     <b>Contact:</b><br>
     <a href="mailto:roman.hillje@ieo.it?subject=Single%20Cell%20Browser">roman.hillje@ieo.it</a><br>
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
  dashboardHeader(title='Single Cell Browser'),
  dashboardSidebar(
    tags$head(tags$style(HTML('.content-wrapper {overflow-x: scroll;}'))),
    sidebarMenu(
      sidebarMenuOutput('sidebar.menu')
    )
  ),
  dashboardBody(
    tags$script(HTML("$('body').addClass('fixed');")),
    tabItems(
      tabItem(tabName='loadData',
        fluidRow(
          column(12,
            titlePanel('Load data'),
            fileInput(inputId     = 'RDS_file',
                      label       = 'Choose RDS file...',
                      multiple    = FALSE,
                      accept      = c('.rds'),
                      width       = NULL,
                      buttonLabel = 'Browse...',
                      placeholder = 'No file selected')
          )
        )
      ),
      tabItem(tabName='overview',
        tagList(
          fluidRow(
            column(width=3, offset=0, style='padding: 0px;',
              box(title='Input parameters', status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                  tagList(
                    uiOutput('overview.UI'),
                    uiOutput('overview.scales')
                  )
              )
            ),
            column(width=9, offset=0, style='padding: 0px;',
              box(title=tagList(p('Dimensional reduction', style='padding-right: 5px; display: inline'),
                                actionButton(inputId='overview.projection.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.', style='margin-right: 5px'),
                                actionButton(inputId='overview.projection.export', label='export to PDF', icon=NULL, class='btn-xs', title='Export dimensional reduction to PDF file.')),
                                #shinyFiles::shinyDirButton(id='overview.projection.export', label='export to PDF', icon=NULL, class='btn-xs', title='Choose directory to export plot to.')),
                  status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                  scatterD3::scatterD3Output('overview.projection', height='720px')
              )
            )
          )
        )
      ),
      tabItem(tabName='samples',
        box(title=tagList(p('Overview of samples', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.table.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            DT::dataTableOutput('samples.table')
        ),
        box(title=tagList(p('Samples by clusters', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.by.cluster.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('samples.by.cluster.UI')
        ),
        box(title=tagList(p('Number of transcripts', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.box.nUMI.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('samples.box.nUMI')
        ),
        box(title=tagList(p('Number of expressed genes', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.box.nGene.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('samples.box.nGene')
        ),
        box(title=tagList(p('Mitochondrial gene expression', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.box.percent.mt.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('samples.box.percent.mt')
        ),
        box(title=tagList(p('Ribosomal gene expression', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.box.percent.ribo.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('samples.box.percent.ribo')
        ),
        box(title=tagList(p('Cell cycle analysis (Regev)', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.by.cell.cycle.Regev.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('samples.by.cell.cycle.Regev.UI')
        ),
        box(title=tagList(p('Cell cycle analysis (Cyclone)', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='samples.by.cell.cycle.Cyclone.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('samples.by.cell.cycle.Cyclone.UI')
        )
      ),
      tabItem(tabName='clusters',
        box(title=tagList(p('Overview of clusters', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.table.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            DT::dataTableOutput('clusters.table')
        ),
        box(title=tagList(p('Cluster tree', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.tree.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('clusters.tree.UI')
        ),
        box(title=tagList(p('Clusters by samples', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.by.sample.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('clusters.by.sample.UI')
        ),
        box(title=tagList(p('Number of transcripts', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.box.nUMI.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('clusters.box.nUMI')
        ),
        box(title=tagList(p('Number of expressed genes', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.box.nGene.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('clusters.box.nGene')
        ),
        box(title=tagList(p('Mitochondrial gene expression', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.box.percent.mt.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('clusters.box.percent.mt')
        ),
        box(title=tagList(p('Ribosomal gene expression', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.box.percent.ribo.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            plotly::plotlyOutput('clusters.box.percent.ribo')
        ),
        box(title=tagList(p('Cell cycle analysis (Regev)', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.by.cell.cycle.Regev.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('clusters.by.cell.cycle.Regev.UI')
        ),
        box(title=tagList(p('Cell cycle analysis (Cyclone)', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='clusters.by.cell.cycle.Cyclone.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('clusters.by.cell.cycle.Cyclone.UI')
        )
      ),
      tabItem(tabName='topExpressedGenes',
        box(title=tagList(p('Top expressed genes per sample', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='top.expressed.genes.by.sample.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('top.expressed.genes.by.sample.UI')
        ),
        box(title=tagList(p('Top expressed genes per cluster', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='top.expressed.genes.by.cluster.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('top.expressed.genes.by.cluster.UI')
        )
      ),
      tabItem(tabName='markerGenes',
        box(title=tagList(p('Marker genes per sample', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='marker.genes.by.sample.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('marker.genes.by.sample.UI')
        ),
        box(title=tagList(p('Marker genes per cluster', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='marker.genes.by.cluster.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('marker.genes.by.cluster.UI')
        )
      ),
      tabItem(tabName='enrichedPathways',
        box(title=tagList(p('Enriched pathways by sample', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='enriched.pathways.by.sample.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('enriched.pathways.by.sample.UI')
        ),
        box(title=tagList(p('Enriched pathways by cluster', style='padding-right: 5px; display: inline'),
                          actionButton(inputId='enriched.pathways.by.cluster.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
            status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
            uiOutput('enriched.pathways.by.cluster.UI')
        )
      ),
      tabItem(tabName='geneExpression',
        tagList(
          fluidRow(
            column(width=3, offset=0, style='padding: 0px;',
              box(title='Input parameters', status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                  tagList(
                    uiOutput('expression.UI'),
                    uiOutput('expression.scales')
                  )
              )
            ),
            column(width=9, offset=0, style='padding: 0px;',
              box(title=tagList(p('Dimensional reduction', style='padding-right: 5px; display: inline'),
                                actionButton(inputId='expression.projection.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.', style='margin-right: 5px'),
                                actionButton(inputId='expression.projection.export', label='export to PDF', icon=NULL, class='btn-xs', title='Export dimensional reduction to PDF file.')),
                  status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                  tagList(
                    scatterD3::scatterD3Output('expression.projection', height='720px'),
                    tags$br(),
                    htmlOutput('expression.genes.displayed')
                  )
              )
            )
          ),
          fluidRow(
            box(title=tagList(p('Expression levels by sample', style='padding-right: 5px; display: inline'),
                              actionButton(inputId='expression.by.sample.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
                status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                plotly::plotlyOutput('expression.by.sample')
            )
          ),
          fluidRow(
            box(title=tagList(p('Expression levels by cluster', style='padding-right: 5px; display: inline'),
                              actionButton(inputId='expression.by.cluster.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
                status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                plotly::plotlyOutput('expression.by.cluster')
            )
          )
        )
      ),
      tabItem(tabName='geneSetExpression',
        tagList(
          fluidRow(
            column(width=3, offset=0, style='padding:0px;',
              box(title='Input parameters', status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                  tagList(
                    uiOutput('geneSetExpression.UI'),
                    uiOutput('geneSetExpression.scales')
                  )
              )
            ),
            column(width=9, offset=0, style='padding:0px;',
              box(title=tagList(p('Dimensional reduction', style='padding-right: 5px; display: inline'),
                                actionButton(inputId='geneSetExpression.projection.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.', style='margin-right: 5px'),
                                actionButton(inputId='geneSetExpression.projection.export', label='export to PDF', icon=NULL, class='btn-xs', title='Export dimensional reduction to PDF file.')),
                  status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                  tagList(
                    scatterD3::scatterD3Output('geneSetExpression.projection', height='720px'),
                    tags$br(),
                    htmlOutput('geneSetExpression.genes.displayed')
                  )
              )
            )
          ),
          fluidRow(
            box(title=tagList(p('Expression levels by sample', style='padding-right: 5px; display: inline'),
                              actionButton(inputId='geneSetExpression.by.sample.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
                status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                plotly::plotlyOutput('geneSetExpression.by.sample')
            )
          ),
          fluidRow(
            box(title=tagList(p('Expression levels by cluster', style='padding-right: 5px; display: inline'),
                              actionButton(inputId='geneSetExpression.by.cluster.info', label='info', icon=NULL, class='btn-xs', title='Show additional information for this panel.')),
                status='primary', solidHeader=TRUE, width=12, collapsible=TRUE,
                plotly::plotlyOutput('geneSetExpression.by.cluster')
            )
          )
        )
      ),
      tabItem(tabName='geneIdConversion',
        box(title='Convert gene ID <-> gene symbol', status='primary', solidHeader=TRUE, width=12, collapsible=FALSE,
          tagList(
            selectInput('geneIdConversion.organism', 'Organism:', choices=c('mouse', 'human')),
            DT::dataTableOutput('gene_info')
          )
        )
      ),
      tabItem(tabName='info',
        fluidPage(
          fluidRow(
            column(12,
              titlePanel('Information about samples and analysis'),
              htmlOutput('sample.info.general')
              #htmlOutput('sample.info.R')
              # htmlOutput('info.data.version'),      # data.version
              # htmlOutput('info.data.type'),         # data.type
              # htmlOutput('info.project.name'),      # project.name
              # htmlOutput('info.organism'),          # organism
              # htmlOutput('info.reference.version'), # reference.version
              # htmlOutput('info.annotation'),        # annotation
              # htmlOutput('info.annotation.type'),   # annotation.type
              # htmlOutput('info.min.genes'),         # min.genes
              # htmlOutput('info.vars.to.regress'),   # vars.to.regress
              # htmlOutput('info.number.PCs'),        # number.PCs
              # htmlOutput('info.tsne.perplexity'),   # tsne.perplexity
              # htmlOutput('info.enrichr.dbs'),       # enrichr.dbs
              # htmlOutput('info.genes.mt'),          # genes.mt
              # htmlOutput('info.genes.ribo'),        # genes.ribo
              # htmlOutput('info.genes.apoptosis'),   # genes.apoptosis
              # htmlOutput('info.genes.S'),           # genes.S
              # htmlOutput('info.genes.G2M')          # genes.G2M
            )
          )
        )
      ),
      tabItem(tabName='about',
        fluidPage(
          fluidRow(
            column(12,
              titlePanel('About this application'),
              htmlOutput('about')
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


