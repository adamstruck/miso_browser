##-----------------------------------
## Utility Functions
##-----------------------------------
parseCoordinates <- function(text) {
    coords <- unlist(lapply(strsplit(text, "[:-]| +"), function(x) { x[x != ""] }))
    if (length(coords) == 3 & as.numeric(coords[2]) < as.numeric(coords[3])) {
        coordinates <- c(as.character(coords[1]), as.numeric(coords[2]), as.numeric(coords[3]))
        return(coordinates)
    } else {
        return(NULL)
    }
}

##-----------------------------------
## Define the server behavior
##-----------------------------------
function(input, output, session) {
    ##-----------------------------------
    ## Reactive inputs
    ##-----------------------------------
    observeEvent(input$reset_inputs, {
        updateSelectInput(session, "misoSelect", choices = "Choose a MISO splicing event")
        updateTextInput(session, "misoEvent", value = "")
    })
    
    ## rinputs <- reactiveValues(coordinates = NULL)    
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
    
    ## Reactive miso data selection 
    MISOdataInput <- reactive({
        MISOdata <- miso_sqlite      
        miso_event <- input$misoEvent
        ## isoforms <- selectIsoform()$isoform
        isoforms <- NULL
        if (!is.null(miso_event) && miso_event != "") {
            event_filter <- interp(~param == miso_event, param = as.name("isoform"))
            MISOdata <- filter_(MISOdata, .dots = event_filter)
        } else if (!is.null(isoforms) && isoforms != "") {
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
        MISOdata <- MISOdata %>% filter(tissue %in% c("Tibialis", "Quad", "Heart"), diagnosis %in% c("DM1", "Control"))
        MISOdata_df <- collect(MISOdata)
        MISOdata_df
    })

    ## Reactive sample metadata selection
    sample_metadata_input <- reactive({
        MISOdata <- MISOdataInput()
        samples <- unique(MISOdata$sampleID)
        if (dim(MISOdata)[1] > 0 && length(samples) > 1) {
            selected_cols <- metadata_cols()        
            sample_metadata %>%
                filter(sampleID %in% samples) %>% 
                select_(.dots = selected_cols)
        } else {
            sample_metadata[0,]
        }
    })

    ## Render gene selection ui
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
    
    ## Render isoform selection ui
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

    
    ##-----------------------------------
    ## Plotting
    ##-----------------------------------
    output$psi_diagnosis_plot <- renderPlot({
        plot_data <- MISOdataInput()
        if (dim(plot_data)[1] > 0) {            
            plot_main <- ggplot(plot_data)
            plot_main <- plot_main +
                geom_point(aes(x = 2^as.numeric(factor(diagnosis))-0.5,
                               y = miso_posterior_mean),
                           position = position_jitter(width = 0.1),
                           size = 3,
                           alpha=0.75) +
                geom_boxplot(aes(x = 2^as.numeric(factor(diagnosis)),
                                 y = miso_posterior_mean,
                                 fill = diagnosis),
                             width = .75,
                             size = 0.75,
                             alpha = 0.75,
                             outlier.shape = NA) +
                labs(y = "PSI", x = "Diagnosis") +
                theme_bw(25) +
                scale_x_continuous(breaks = c(1.75, 3.75), labels = c("Control", "DM1"), limits=c(1.25, 4.5)) +
                scale_fill_discrete(guide = FALSE) +
                facet_wrap(~tissue, nrow=1, scales = "fixed")            
        } else {
            plot_main <- ggplot(plot_data, aes_string(x = "diagnosis", y = "miso_posterior_mean"))
            plot_main <- plot_main +
                geom_blank() +
                labs(title = "No Data Available" , x = "Diagnosis", y = "PSI") +
                theme_bw(25) +
                theme(plot.title = element_text(vjust=-35, size = 80, colour = "red"))
        }
        print(plot_main)
    })

    
    output$psi_sample_plot <- renderPlot({
        plot_data <- MISOdataInput()
        if (dim(plot_data)[1] > 0) {            
            plot_main <- plot_data %>% arrange(miso_posterior_mean) %>%
                ggplot(aes(x = factor(sampleID, levels = sampleID),
                           y = miso_posterior_mean))
            plot_main <- plot_main +
                geom_point(aes(colour = diagnosis),
                               size = 3,
                           alpha=0.75) +
                geom_errorbar(aes(ymax = ci_high, ymin = ci_low),
                              width = 0.25,
                              alpha = 0.75) + 
                labs(x = "Sample", y = "PSI") +
                facet_grid(~tissue, scales = "free_x", space = "free_x") +
                theme_bw(25) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                      legend.position = "top")
        } else {
            plot_main <- ggplot(plot_data, aes_string(x = "sampleID", y = "miso_posterior_mean"))
            plot_main <- plot_main +
                geom_blank() +                
                labs(title = "No Data Available", x = "Sample", y = "PSI") +
                theme_bw(25) +
                theme(plot.title = element_text(vjust=-35, size = 80, colour = "red"))
        }
        print(plot_main)
    })


    ##-----------------------------------
    ## Showing / Downloading Data
    ##-----------------------------------    

    ## Reactive column selection
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
    
    ## Conditional variable selection panels
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

    
    ## Output miso data table
    output$miso_data_ui <- renderTable({
        MISOdataInput()
    })

    ## output$miso_data_ui <- renderUI(DT::dataTableOutput(outputId = "miso_data"))
    
    ## output$miso_data <- DT::renderDataTable({
    ##     DT::datatable(MISOdataInput(),
    ##                   rownames = FALSE,
    ##                   options = list(searching = FALSE,
    ##                                  lengthMenu = c(25, 50, 100, 200),
    ##                                  pageLength = 25))
    ## })

    ##
    ## Output sample metadata table
    ##
    output$sample_metadata_ui <- renderTable({
        sample_metadata_input()
    })
    
    ## output$sample_metadata_ui <- renderUI({DT::dataTableOutput(outputId = "sample_metadata")})
    
    ## output$sample_metadata <- DT::renderDataTable({
    ##     DT::datatable(sample_metadata_input(),
    ##                   rownames = FALSE,
    ##                   options = list(searching = FALSE,
    ##                                  lengthMenu = list(c(25, 50, -1), c("25", "50", "All")),
    ##                                  pageLength = 25))
    ## })
    
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
