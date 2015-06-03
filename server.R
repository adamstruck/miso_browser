##--------------------------------------------------------------
## Setup enviroment
##--------------------------------------------------------------
parseCoordinates <- function(text) {
    coords <- unlist(lapply(strsplit(text, "[:-]| +"), function(x) { x[x != ""] }))
    if (length(coords) == 3 & as.numeric(coords[2]) < as.numeric(coords[3])) {
        coordinates <- c(as.character(coords[1]), as.numeric(coords[2]), as.numeric(coords[3]))
        return(coordinates)
    } else {
        return(NULL)
    }
}

##--------------------------------------------------------------
## Define the server behavior
##--------------------------------------------------------------
function(input, output, session) {
    ##
    ## Reactive inputs
    ##
    ## rinputs <- reactiveValues(coordinates = NULL, gene = NULL, isoform = NULL, miso_event = NULL)

    ## observe({
    ##     rinputs$coordinates <- parseCoordinates(input$coordinates)
    ## })
    
    ## selectGene <- eventReactive(
    ##     rinputs$coordinates,
    ## {
    ##     coords <- rinputs$coordinates
    ##     coords_granges <- GRanges(seqnames = coords[1], ranges = IRanges(as.numeric(coords[2]), as.numeric(coords[3])))
    ##     hits <- findOverlaps(coords_granges, MISOevent_lookup_granges)
    ##     event_lookup_filtered <- MISOevent_lookup_granges[subjectHits(hits),] %>% as.data.frame()
    ##     genes <- event_lookup_filtered$symbol %>% unique()
    ##     genes
    ## }
    ## )
    
    ## selectIsoform <- eventReactive(
    ##     input$gene_symbol,
    ## {
    ##     gene_filter <- interp(~param == input$gene_symbol, param = as.name("symbol"))
    ##     iso_lookup <- filter_(MISOevent_lookup, .dots = gene_filter) %>% select(isoform, desc_iso)
    ##     iso_lookup
    ## }
    ## )
    
    ## miso_event <- eventReactive(
    ##     input$isoform,
    ## {
    ##     iso_lookup <- selectIsoform()
    ##     iso_lookup <- dplyr::filter(iso_lookup, desc_iso == input$isoform)
    ##     ## desc-iso not unique... need to fix
    ##     iso_lookup$isoform[1]
    ## }
    ## )
    
    ##
    ## Reactive miso data selection 
    ##

    MISOdataInput <- reactive({
        MISOdata <- miso_sqlite      
        miso_event <- input$misoEvent
        ## isoforms <- selectIsoform()$isoform
        isoforms <- NULL
        if (!is.null(miso_event) || miso_event != "") {
            event_filter <- interp(~param == miso_event, param = as.name("isoform"))
            MISOdata <- filter_(MISOdata, .dots = event_filter)
        } else if (!is.null(isoforms) || isoforms[1] != "") {
            if (length(isoforms) == 1) {
                event_filter <- interp(~param == isoforms, param = as.name("isoform"))
            } else {
                event_filter <- interp(~param %in% isoforms, param = as.name("isoform"))
            }
            MISOdata <- filter_(MISOdata, .dots = event_filter)
        } else {
            MISOdata <- head(MISOdata)[0,]
        }
        selected_cols <- miso_cols()
        MISOdata <- select_(MISOdata, .dots = selected_cols)
        MISOdata_df <- collect(MISOdata)
        MISOdata_df
    })

    ##
    ## Render gene selection ui
    ##
    ## output$select_gene <- renderUI({
    ##     if (is.null(rinputs$coordinates)) {
    ##         selectizeInput(inputId = "gene_symbol",
    ##                        label = "Gene:",
    ##                        choices = gene_list,
    ##                        multiple = FALSE,
    ##                        options = list(
    ##                            placeholder = "",
    ##                            onInitialize = I('function() { this.setValue(""); }'),
    ##                            maxOptions = 5)
    ##                        )
    ##     } else {
    ##         selectizeInput(inputId = "gene_symbol",
    ##                        label = "Gene:",
    ##                        choices = selectGene(),
    ##                        multiple = FALSE,
    ##                        options = list(
    ##                            placeholder = "",
    ##                            onInitialize = I('function() { this.setValue(""); }'),
    ##                            maxOptions = 5)
    ##                        )
    ##     } 
    ## })
    
    ## ##
    ## ## Render isoform selection ui
    ## ##
    ## output$select_isoform <- renderUI({
    ##     selectizeInput(inputId = "isoform",
    ##                    label = "Select Isoform:",
    ##                    choices = selectIsoform()$desc_iso,
    ##                    multiple = FALSE,
    ##                    options = list(
    ##                        placeholder = "",
    ##                        onInitialize = I('function() { this.setValue(""); }'))                       
    ##                    )
    ## })

    ##
    ## Plot psi values
    ##
    output$jitter_plot <- renderPlot({
        miso_event <- input$misoEvent
        plot_main <- ggplot(MISOdataInput(), aes_string(x = "tissue", y = "miso_posterior_mean"))
        if (miso_event != "") {
            plot_main <- plot_main +
                geom_point(aes_string(color = "tissue", shape = "diagnosis"),
                           alpha = 1, size = 3, position = position_jitter(width=0.1, height=0)) +
                scale_colour_discrete(guide = FALSE) +
                labs(x = "Tissue", y = "PSI") +
                theme_bw(15) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            } else {
                plot_main <- plot_main +
                    geom_blank() +
                    labs(x = "Tissue", y = "PSI") +
                    theme_bw(15) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }
        print(plot_main)
    })

    output$box_plot <- renderPlot({
        miso_event <- input$misoEvent
        plot_main <- ggplot(MISOdataInput(), aes_string(x = "diagnosis", y = "miso_posterior_mean"))
        if (miso_event != "") {
            plot_main <- plot_main +
                geom_boxplot(aes_string(fill = "diagnosis")) +
                labs(x = "Diagnosis", y = "PSI") +
                scale_fill_discrete(guide = FALSE) +
                facet_wrap(~tissue) +
                theme_bw(15)
        } else {
            plot_main <- plot_main +
                geom_blank() +
                labs(x = "Diagnosis", y = "PSI") +
                theme_bw(15)
        }
        print(plot_main)
    })
    
    ##
    ## Reactive column selection
    ##
    miso_cols <- reactive({
        if (is.null(input$miso_show_cols)) {
            MISOdata_columns
        } else {
            input$miso_show_cols
        }
    })

    metadata_cols <- reactive({
        if (is.null(input$metadata_show_cols)) {
            sample_metadata_columns
        } else {
            input$metadata_show_cols
        }
    })

    ##
    ## Reactive sample metadata selection
    ##
    sample_metadata_input <- reactive({
        selected_cols <- metadata_cols()
        select_(sample_metadata, .dots = selected_cols)
    })

    ##
    ## Conditional panels
    ##
    output$show_vars <- renderUI({
        if (input$raw_data_box == "MISO Data") {
            wellPanel(
                style = "overflow: auto",
                checkboxGroupInput("miso_show_cols",
                                   "Columns to display:",
                                   MISOdata_columns,
                                   selected = MISOdata_columns)
            )
        } else if (input$raw_data_box == "Sample Metadata") {
            wellPanel(
                style = "overflow: auto",
                checkboxGroupInput("metadata_show_cols",
                                   "Columns to display:",
                                   sample_metadata_columns,
                                   selected = sample_metadata_columns)
            )
        }
    })

    output$download <- renderUI({
        if (input$raw_data_box == "MISO Data") {
            wellPanel(
                style = "text-align: center; overflow: auto",
                downloadButton(outputId = "download_miso_data",
                               label = "Download Selected"),
                br(),
                br()
            )
        } else if (input$raw_data_box == "Sample Metadata") {
            wellPanel(
                style = "text-align: center; overflow: auto",
                downloadButton(outputId = "download_sample_metadata",
                               label = "Download Selected"),
                br(),
                br()
            )
        }
    })

    ##
    ## Output miso data table
    ##
    output$miso_data <- DT::renderDataTable({
        ## input$submit
        ## isolate({
            data <- MISOdataInput()
            DT::datatable(data,
                          rownames = FALSE,
                          options = list(searching = FALSE,
                                         lengthMenu = c(25, 50, 100, 200),
                                         pageLength = 25))
        ## })
    })

    ##
    ## Output sample metadata table
    ##
    output$sample_metadata <- DT::renderDataTable({
        data <- sample_metadata_input()
        DT::datatable(data,
                      rownames = FALSE,
                      options = list(searching = FALSE,
                                     lengthMenu = list(c(25, 50, -1), c("25", "50", "All")),
                                     pageLength = 25))
    })

    ##
    ## Download handler
    ##
    output$download_miso_data <- downloadHandler(
        filename = "miso_data_selected.csv",
        content = function(file) {
        write.csv(MISOdataInput(),
                  file,
                  col.names=TRUE,
                  row.names=FALSE,
                  quote=FALSE)
    }
    )

    ##
    ## Download handler
    ##
    output$download_sample_metadata <- downloadHandler(
        filename = "sample_metadata_selected.csv",
        content = function(file) {
        write.csv(sample_metadata_input(),
                  file,
                  col.names=TRUE,
                  row.names=FALSE,
                  quote=FALSE)
    }
    )
}
