library("shiny")
library("BiocStyle")
library("bs4Dash")
library("shinyjs")
library("shinycssloaders")
library("plotly")
library("SummarizedExperiment")
library("visNetwork")
library("ggraph")
library("tidygraph")
library("ggkegg")
source("/Users/edoardofilippi/Development/Projects/BasicApp/enrichment_map.r")
source("/Users/edoardofilippi/Development/Projects/BasicApp/gs_heatmap.R")
MovidaApp <- function(movida_list = list(),
                      project_id = "") {
  # https://projects.lukehaas.me/css-loaders/
  # or even think of https://cran.r-project.org/web/packages/shinycustomloader/README.html
  #   oopt <- options(spinner.type = 6, spinner.color = .biocgreen)
  #   # play nice with other previously chosen options
  #   on.exit(options(oopt))

  usage_mode <- "shiny_mode"
  # UI definition -----------------------------------------------------------
  movida_ui <- bs4Dash::bs4DashPage(
    title = "Movida",
    header = bs4Dash::bs4DashNavbar(
      tagList(
        tags$code(tags$h3(style = "color: blue; font-family: Arial, sans-serif;", "MOViDA")),
        div(
          style = "margin-left: 20px; width: 40%;", # Set width in the div wrapper (e.g., 50% width)
          uiOutput("contrast_selector")
        )
      )
    ),
    sidebar = bs4DashSidebar(
      skin = "light",
      status = "primary",
      title = "Movida",
      bs4SidebarMenu(
        id = "gt_tabs",
        bs4SidebarMenuItem("Data Loading", tabName = "tab_data", icon = icon("upload")),
        # bs4SidebarMenuItem("QC", tabName = "tab_qc", icon = icon("check-circle")),
        bs4SidebarMenuItem("Overview", tabName = "tab_overview", icon = icon("eye")),
        bs4SidebarMenuItem("Metabolomics", tabName = "tab_metabolomics", icon = icon("flask")),
        bs4SidebarMenuItem("Plots", tabName = "tab_plots", icon = icon("chart-bar")),
        bs4SidebarMenuItem("KeggMaps", tabName = "tab_kegg", icon = icon("map")),
        bs4SidebarMenuItem("Integration", tabName = "tab_late_integration", icon = icon("puzzle-piece")),
        bs4SidebarMenuItem("Multi-Omics", tabName = "tab_multiomics", icon = icon("layer-group"))
      )
    ),
    body = bs4Dash::bs4DashBody(
      rintrojs::introjsUI(),
      bs4TabItems(

        # Data Loading panel -----------------------------------------------------------
        bs4TabItem(
          tabName = "tab_data",
          uiOutput("ui_panel_welcome")
        ),

        # Overview panel --------------------------------------------------------
        bs4TabItem(
          tabName = "tab_overview",
          uiOutput("ui_panel_overview")
        ),

        # Metabolomics panel --------------------------------------------------------
        bs4TabItem(
          tabName = "tab_metabolomics",
          uiOutput("ui_panel_metabolomics")
        ),

        # Plots panel ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_plots",
          uiOutput("ui_panel_plots")
        ),

        # KeggMaps panel ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_kegg",
          uiOutput("ui_panel_kegg")
        ),

        # Late Integration panel ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_late_integration",
          uiOutput("ui_late_integration")
        ),

        # Multi-Omics panel ------------------------------------------------------
        bs4TabItem(
          tabName = "tab_multiomics",
          uiOutput("ui_panel_multiomics")
        )
      )
    ),
    controlbar = bs4Dash::bs4DashControlbar(
      collapsed = TRUE,
      uiOutput("ui_controlbar")
    ),
    footer = bs4DashFooter(
      left = NULL,
      right = NULL
    )
  )


  # Server definition -----------------------------------------------------------
  movida_server <- function(input, output, session) {
    #### Basic logic ------------------------------------------------------
    reactive_values <- shiny::reactiveValues(
      annotation_df_prot = movida_list$annotation_df_prot,
      annotation_df_trans = movida_list$annotation_df_trans,
      annotation_df_metabo = movida_list$annotation_df_metabo,
      results_prot = movida_list$results_prot,
      results_trans = movida_list$results_trans,
      results_metabo = movida_list$results_metabo,
      se_prot = movida_list$se_prot,
      se_trans = movida_list$se_trans,
      se_metabo = movida_list$se_metabo,
      metadata = movida_list$metadata
    )

    # Compute the intersection reactively
    common_contrasts <- reactive({
      req(reactive_values$results_prot, reactive_values$results_trans)
      intersect(
        names(reactive_values$results_prot),
        names(reactive_values$results_trans)
      )
    })

    # Render the select input dynamically
    output$contrast_selector <- renderUI({
      req(common_contrasts())

      selectInput(
        inputId = "selected_contrast",
        label = "Select Contrast",
        choices = common_contrasts(),
        selected = common_contrasts()[1],
        width = "100%"
      )
    })
    #### Data loading panel ------------------------------------------------------
    output$ui_panel_welcome <- renderUI({
      tagList(
        fluidRow(
          bs4Dash::column(
            width = 10,
            offset = 1,
            tags$div(
              tags$p("In this section there will be:"),
              tags$p("A welcome page with a short description of the functionalities"),
              tags$p("Possibility to upload data if it was not done programmatically")
            )
          )
        )
      )
    })

    #### Overview panel ------------------------------------------------------
    output$ui_panel_overview <- renderUI({
      tagList(
        fluidRow(
          bs4Dash::column(
            width = 12,
            offset = 0,
            bs4Dash::bs4Card(
              id = "tab_overview_card_trans",
              title = "rnaSeq",
              elevation = 2,
              width = 12,
              closable = FALSE,
              bs4Dash::tabsetPanel(
                id = "tsp1",
                type = "pills",
                selected = "DE",
                side = "right",
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px", # You can adjust max-height as needed
                  title = "DE",
                  div(
                    style = "margin-top: 50px; transform: scale(0.8); transform-origin: top left",
                    DT::dataTableOutput("gs_de_table_trans")
                  )
                ),
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px", # You can adjust max-height as needed
                  title = "Enrichment",
                  div(
                    style = "margin-top: 50px; transform: scale(0.8); transform-origin: top left",
                    DT::dataTableOutput("gs_enrich_table_trans")
                  )
                ),
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px", # You can adjust max-height as needed
                  title = "PCA?"
                )
              )
            )
          )
        ),
        fluidRow(
          bs4Dash::column(
            width = 12,
            offset = 0,
            bs4Dash::bs4Card(
              id = "tab_overview_card_prot",
              title = "Proteomics",
              elevation = 2,
              width = 12,
              closable = FALSE,
              bs4Dash::tabsetPanel(
                id = "tsp1",
                type = "pills",
                selected = "DE",
                side = "right",
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px", # You can adjust max-height as needed
                  title = "DE",
                  div(
                    style = "margin-top: 50px; transform: scale(0.8); transform-origin: top left",
                    DT::dataTableOutput("gs_de_table_prot")
                  )
                ),
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px", # You can adjust max-height as needed
                  title = "Enrichment",
                  div(
                    style = "margin-top: 50px; transform: scale(0.8); transform-origin: top left",
                    DT::dataTableOutput("gs_enrich_table_prot")
                  )
                ),
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px", # You can adjust max-height as needed
                  title = "PCA?"
                )
              )
            )
          )
        )
      )
    })

    # Render filtered table based on the selected contrast
    output$gs_de_table_trans <- DT::renderDT({
      req(input$selected_contrast) # Ensure the input is available

      # Filter the table based on the selected contrast
      filtered_table <- reactive_values$results_trans[[input$selected_contrast]]$etbl_res_DE

      DT::datatable(filtered_table,
        escape = FALSE,
        rownames = TRUE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "excel",
              text = "Export to Excel",
              title = "DE_Proteomics_Export"
            )
          ),
          pageLength = 10
        ),
        class = "display nowrap"
      )
    })
    # Render filtered table based on the selected contrast
    output$gs_enrich_table_trans <- DT::renderDT({
      req(input$selected_contrast) # Ensure the input is available

      # Filter the table based on the selected contrast
      filtered_table <- reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl
      # Return the filtered table as a datatable
      DT::datatable(filtered_table,
        escape = FALSE,
        rownames = TRUE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "excel",
              text = "Export to Excel",
              title = "DE_Proteomics_Export"
            )
          ),
          pageLength = 10
        ),
        class = "display nowrap"
      )
    })
    # Render filtered table based on the selected contrast
    output$gs_de_table_prot <- DT::renderDT({
      req(input$selected_contrast) # Ensure the input is available

      # Filter the table based on the selected contrast
      filtered_table <- reactive_values$results_trans[[input$selected_contrast]]$etbl_res_DE

      # Return the filtered table as a datatable
      DT::datatable(filtered_table,
        escape = FALSE,
        rownames = TRUE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "excel",
              text = "Export to Excel",
              title = "DE_Proteomics_Export"
            )
          ),
          pageLength = 10
        ),
        class = "display nowrap"
      )
    })

    # Render filtered table based on the selected contrast
    output$gs_enrich_table_prot <- DT::renderDT({
      req(input$selected_contrast) # Ensure the input is available

      # Filter the table based on the selected contrast
      filtered_table <- reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl

      DT::datatable(filtered_table,
        escape = FALSE,
        rownames = TRUE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(
              extend = "excel",
              text = "Export to Excel",
              title = "DE_Proteomics_Export"
            )
          ),
          pageLength = 10
        ),
        class = "display nowrap"
      )
    })

    #### Metabolomics panel ------------------------------------------------------
    output$ui_panel_metabolomics <- renderUI({
      tagList(
        fluidRow(
          bs4Dash::column(
            width = 10,
            offset = 1,
            tags$div(
              tags$p("In this section there will be:"),
              tags$p("Since Metabolomics are quite different from the other two there can be some metabolomic specific analysiss, for example:"),
              tags$ul(
                tags$li("Metabolite class activity, Metabolite hub analysis. (PaintOmics4, Liu et al., 2022)"),
                tags$li("PCA, PLSDA, Volcano, Clustering, Enrichment... (Metabolon like analysis)"),
              )
            )
          )
        )
      )
    })

    #### Plots panel ------------------------------------------------------
    output$ui_panel_plots <- renderUI({
      req(input$selected_contrast) # Ensure the input is available

      tagList(
        fluidRow(
          bs4Dash::column(
            width = 12,
            offset = 0,
            bs4Dash::bs4Card(
              id = "tab_overview_card_trans",
              title = "rnaSeq",
              elevation = 2,
              width = 12,
              closable = FALSE,
              bs4Dash::tabsetPanel(
                id = "tsp1",
                type = "pills",
                selected = "EnrichmentMap",
                side = "right",
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px;",
                  title = "EnrichmentMap",
                  div(
                    style = "margin-top: 20px; padding: 10px; transform-origin: top left;",
                    visNetwork::visNetworkOutput("gs_de_network_trans", height = "500px")
                  )
                ),
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 800px",
                  title = "Pathway Heatmap",
                  div(
                    style = "margin-top: 40px; transform: scale(0.8); transform-origin: top left",
                    tagList(
                      div(
                        style = "margin-top: 40px;",
                        tagList(
                          # Selector UI (left)
                          div(
                            style = "display: inline-block; width: 48%; vertical-align: top;",
                            selectInput(
                              inputId = "pathway_select_trans",
                              label = "Select Pathway",
                              choices = reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl$GO.ID,
                              selected = reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl$GO.ID[1],
                              width = "100%"
                            )
                          ),
                          # Display mode radio buttons (right)
                          div(
                            style = "display: inline-block; width: 48%; vertical-align: top; margin-left:20px",
                            radioButtons(
                              inputId = "contrast_display_mode_trans",
                              label = "Display Mode",
                              choices = list(
                                "Show All" = "all",
                                "Show Only Selected Contrast" = "selected"
                              ),
                              selected = "all",
                              inline = TRUE
                            )
                          )
                        )
                      ),
                      div(
                        shinycssloaders::withSpinner(
                          div(
                            style = "display: flex; justify-content: center; margin-left: 200px; transform: scale(1.2); margin-top: 50px;",
                            plotOutput("gs_heatmap_trans")
                          ),
                          type = 4
                        )
                      )
                    )
                  )
                ),
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px",
                  title = "Heatmap",
                    tags$div(
                    style = "margin-top: 20px;",
                    tags$p("In this section there will be a heatmap where the genes are taken from either other tabs or given as input freely"),
                    tags$p("Another idea is to add a heatmap with multiple pathways"),
                    tags$p("We could also add the horizon plot or anything that can be useful")
                    )
                  # div(
                  #   style = "margin-top: 50px; transform: scale(0.8); transform-origin: top left",
                  #   DT::dataTableOutput("gs_heatmap_trans_genes")
                  # )
                )
              )
            )
          )
        ),
        fluidRow(
          bs4Dash::column(
            width = 12,
            offset = 0,
            bs4Dash::bs4Card(
              id = "tab_overview_card_trans",
              title = "Proteomics",
              elevation = 2,
              width = 12,
              closable = FALSE,
              bs4Dash::tabsetPanel(
                id = "tsp1",
                type = "pills",
                selected = "EnrichmentMap",
                side = "right",

                # EnrichmentMap Tab
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 800px;",
                  title = "EnrichmentMap",
                  div(
                    style = "margin-top: 20px; padding: 10px; transform-origin: top left; transform: scale(1);",
                    visNetwork::visNetworkOutput("gs_de_network_prot", height = "500px")
                  )
                ),

                # Pathway Heatmap Tab
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 800px",
                  title = "Pathway Heatmap",
                  div(
                    style = "margin-top: 40px;",
                    tagList(
                      # Selector UI (left)
                      div(
                        style = "display: inline-block; width: 48%; vertical-align: top;",
                        selectInput(
                          inputId = "pathway_select_prot",
                          label = "Select Pathway",
                          choices = reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl$GO.ID,
                          selected = reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl$GO.ID[1],
                          width = "100%"
                        )
                      ),
                      # Display mode radio buttons (right)
                      div(
                        style = "display: inline-block; width: 48%; vertical-align: top; margin-left:20px",
                        radioButtons(
                          inputId = "contrast_display_mode_prot",
                          label = "Display Mode",
                          choices = list(
                            "Show All" = "all",
                            "Show Only Selected Contrast" = "selected"
                          ),
                          selected = "all",
                          inline = TRUE
                        )
                      )
                    )
                  ),
                  div(
                    shinycssloaders::withSpinner(
                      div(
                        style = "display: flex; justify-content: center; margin-left: 200px; transform: scale(1.2); margin-top: 50px;",
                        plotOutput("gs_heatmap_prot")
                      ),
                      type = 4
                    )
                  )
                ),

                # Gene Table Tab
                shiny::tabPanel(
                  style = "overflow-x: scroll; overflow-y: auto; max-height: 500px",
                  title = "Heatmap",
                  div(
                    style = "margin-top: 50px; transform: scale(0.8); transform-origin: top left",
                    DT::dataTableOutput("gs_heatmap_prot_genes")
                  )
                )
              )
            )
          )
        )
      )
    })

    output$gs_de_network_trans <- renderVisNetwork({
      req(input$selected_contrast) # Ensure the input is available

      # Filter the table based on the selected contrast
      enrichresults_table <- GeneTonic::shake_topGOtableResult(reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl)

      # Assuming enrichment_map() returns an igraph object
      g <- enrichment_map(enrichresults_table) # Call the function that returns the igraph object

      # Convert the igraph object to visNetwork format
      network_data <- toVisNetworkData(g)

      # Generate the network visualization
      visNetwork::visNetwork(nodes = network_data$nodes, edges = network_data$edges) %>%
        visNetwork::visIgraphLayout() # Layout for igraph visualization
    })
    output$gs_de_network_prot <- renderVisNetwork({
      req(input$selected_contrast) # Ensure the input is available

      # Filter the table based on the selected contrast
      enrichresults_table <- GeneTonic::shake_topGOtableResult(reactive_values$results_prot[[input$selected_contrast]]$topGO_tbl)

      # Assuming enrichment_map() returns an igraph object
      g <- enrichment_map(enrichresults_table) # Call the function that returns the igraph object

      # Convert the igraph object to visNetwork format
      network_data <- toVisNetworkData(g)

      # Generate the network visualization
      visNetwork::visNetwork(nodes = network_data$nodes, edges = network_data$edges) %>%
        visNetwork::visIgraphLayout() # Layout for igraph visualization
    })

    # Observe the change in the selected pathway and update the visualization
    observeEvent(list(input$pathway_select_trans, input$contrast_display_mode_trans), {
      # Get the selected pathway
      selected_pathway <- input$pathway_select_trans

      req(input$selected_contrast) # Ensure the input is available
      req(input$contrast_display_mode_trans) # Ensure the input is available
      # Filter the table based on the selected contrast
      shaked_enrich_table <- GeneTonic::shake_topGOtableResult(reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl)

      if (selected_pathway %in% shaked_enrich_table$gs_id) {
        gene_list <- shaked_enrich_table[shaked_enrich_table$gs_id == selected_pathway, "gs_genes"]
        g_split <- strsplit(gene_list, split = ",", fixed = TRUE)[[1]]

        se <- reactive_values$se_trans

        group_label_map <- c(
          "ATII" = "AT2",
          "Finerenon" = "finerenon",
          "SGLT2I" = "sglt2i",
          "ATII_Finerenon" = "AT2finerenon",
          "ATII_SGLT2I" = "AT2sglt2i",
          "ctrl" = "ctrl"
        )

        if (input$contrast_display_mode_trans == "selected") {
          s <- input$selected_contrast
          parts <- strsplit(s, "_vs_")[[1]]
          group1_label <- group_label_map[[parts[1]]]
          group2_label <- group_label_map[[parts[2]]]

          se <- se[, colData(se)$group %in% c(group1_label, group2_label)]
        }


        thisset_members_ids <- reactive_values$annotation_df_trans$ensembl_gene_id[match(g_split, reactive_values$annotation_df_trans$gene_symbol)]
        order <- rownames(colData(se)[order(colData(se)$group, decreasing = FALSE), ])
        order_numbers <- match(order, rownames(colData(se)))
        # message(paste0( "the class is:" ,class(se)))
        # assay(se)
        vst_data <- DESeq2::vst(assay(se))

        if (length(intersect(thisset_members_ids, rownames(vst_data)))) {
          message("creating plot")
        }

        output$gs_heatmap_trans <- renderPlot(
          {
            gs_heatmap(
              mydata = vst_data, # [order,],
              se = se, # [order,],
              res_enrich = shaked_enrich_table,
              #  genelist = reactive_values$annotation_df_trans$gene_id[match(gene_signature, reactive_values$annotation_df_trans$gene_name)],
              geneset_id = selected_pathway,
              annotation_obj = reactive_values$annotation_df_trans,
              # gtl = gtl,
              scale_row = T,
              show_row_dend = F,
              show_column_names = T,
              column_order = order_numbers,
              # row_km = 2,
              anno_col_info = "group"
            )
          },
          height = 400,
          width = 800
        ) # Adjust height and width to make the image bigger
      } else {
        message("GO term not found in the enrichment results")
      }
    })

    # Observe the change in the selected pathway and update the visualization
    observeEvent(input$pathway_select_prot, {
      # Get the selected pathway
      selected_pathway <- input$pathway_select_prot

      req(input$selected_contrast) # Ensure the input is available

      # Filter the table based on the selected contrast
      shaked_enrich_table <- GeneTonic::shake_topGOtableResult(reactive_values$results_prot[[input$selected_contrast]]$topGO_tbl)

      if (selected_pathway %in% shaked_enrich_table$gs_id) {
        gene_list <- shaked_enrich_table[shaked_enrich_table$gs_id == selected_pathway, "gs_genes"]
        g_split <- strsplit(gene_list, split = ",", fixed = TRUE)[[1]]

        thisset_members_ids <- reactive_values$annotation_df_prot$ensembl_gene_id[match(g_split, reactive_values$annotation_df_prot$gene_symbol)]
        order <- rownames(colData(reactive_values$se_prot)[order(colData(reactive_values$se_prot)$group, decreasing = FALSE), ])
        order_numbers <- match(order, rownames(colData(reactive_values$se_prot)))
        # message(paste0( "the class is:" ,class(reactive_values$se_prot)))
        # assay(reactive_values$se_prot)
        vst_data <- edgeR::cpm(assay(reactive_values$se_prot))
        if (length(intersect(thisset_members_ids, rownames(vst_data)))) {
          # message("creating plot")
          output$gs_heatmap_prot <- renderPlot(
            {
              gs_heatmap(
                mydata = vst_data, # [order,],
                se = reactive_values$se_prot, # [order,],
                res_enrich = shaked_enrich_table,
                #  genelist = reactive_values$annotation_df_prot$gene_id[match(gene_signature, reactive_values$annotation_df_prot$gene_name)],
                geneset_id = selected_pathway,
                annotation_obj = reactive_values$annotation_df_prot,
                # gtl = gtl,
                scale_row = T,
                show_row_dend = F,
                show_column_names = T,
                column_order = order_numbers,
                # row_km = 2,
                anno_col_info = "group"
              )
            },
            height = 400,
            width = 800
          )
        } # Adjust height and width to make the image bigger
      } else {
        message("GO term not found in the enrichment results")
      }
    })


    #### KEGG map panel ------------------------------------------------------
    output$ui_panel_kegg <- renderUI({
      tagList(
        fluidRow(
          bs4Dash::column(
            width = 3,
            div(
              style = "padding: 10px;",
              selectInput(
                inputId = "kegg_pathway_select",
                label = "Choose KEGG Pathway",
                choices = c("hsa01100", "hsa01200", "hsa01210", "hsa01212", "hsa01230", "hsa01232", "hsa01250", "hsa01240", "hsa01320", "hsa00010", "hsa00020", "hsa00030", "hsa00040", "hsa00051", "hsa00052", "hsa00053", "hsa00500", "hsa00620", "hsa00630", "hsa00640", "hsa00650", "hsa00562", "hsa00190", "hsa00910", "hsa00920", "hsa00061", "hsa00062", "hsa00071", "hsa00100", "hsa00120", "hsa00140", "hsa00561", "hsa00564", "hsa00565", "hsa00600", "hsa00590", "hsa00591", "hsa00592", "hsa01040", "hsa00230"),
                selected = "hsa04110",
                width = "100%"
              ), # This will render the selectInput dynamically
              actionButton("update_kegg", "Update Plot", class = "btn btn-primary")
            )
          ),
          bs4Dash::column(
            width = 9,
            ggiraph::girafeOutput("kegg_map_plot", width = "100%", height = "500px")
          )
        )
      )
    })

    # This is the plot rendering logic. When the button is clicked, this will update.
    observeEvent(input$update_kegg, {
      # Make sure the selected pathway is available before proceeding
      req(input$kegg_pathway_select)

      pathway_id <- input$kegg_pathway_select # Get the selected pathway ID

      # Generate the plot based on the selected pathway
      p <- pathway(pathway_id) |> ## Obtain and parse the KEGG pathway dynamically
        activate(nodes) |>
        mutate(
          convert_hsa = convert_id("hsa"),
          convert_map = convert_id("pathway")
        ) |>
        ggraph(x = x, y = y) +
        geom_edge_parallel(
          arrow = arrow(length = unit(1, "mm")),
          aes(linetype = subtype_name),
          end_cap = circle(7.5, "mm")
        ) +
        geom_node_rect(
          aes(filter = type == "gene", fill = I(bgcolor)),
          color = "black"
        ) +
        geom_node_text(
          aes(label = convert_hsa),
          size = 2, family = "serif"
        ) +
        geom_node_text(
          aes(label = convert_hsa, filter = !is.na(convert_hsa) & convert_hsa == "TP53"),
          size = 2, color = "red", family = "serif"
        ) +
        theme_void()

      # Render the plot with ggiraph for interaction
      output$kegg_map_plot <- ggiraph::renderGirafe({
        ggiraph::girafe(
          ggobj = p, # Pass the `p` plot object here
          options = list(
            ggiraph::opts_hover(css = "fill:red;"),
            ggiraph::opts_sizing(rescale = TRUE),
            ggiraph::opts_selection(type = "multiple")
          )
        )
      })
    })

    #### Multiomics panel ------------------------------------------------------
    output$ui_late_integration <- renderUI({
      tagList(
        fluidRow(
          bs4Dash::column(
            width = 10,
            offset = 1,
            tags$div(
              tags$p("In this section there will be:"),
              tags$p("The idea here is to do p-value based integration at pathway or gene level, the goal is to discover common patways and related entities from the different datasets. "),
              tags$p("This can also increase the confidence in the findings when all the layers point in the same direction. "),
              tags$p("There is also a possibility to work with graph transformation methods, using different strategies to contruct a gene / compound network from existing knoweledge about patwhay topology ")
            )
          )
        )
      )
    })

    #### Integration panel ------------------------------------------------------
    output$ui_panel_multiomics <- renderUI({
      tagList(
        fluidRow(
          bs4Dash::column(
            width = 10,
            offset = 1,
            tags$div(
              tags$p("In this section there will be:"),
              tags$p("One or possibly more integration methods, that will result in a clustering and / or dimensionality reduction, the results will be then used to used to select the entities that are realted to each other in each dataset. "),
              tags$p("This page will be able to interact with the others to then use the insights acquired to generate plots / export tables. "),
              tags$p("Methods that could be implemented:"),
              tags$ul(
                tags$li("MOFA"), # These are latent space based methods
                tags$li("moClusters"),
                tags$li("DIABLO"),
                tags$li("Network based: Similarity network fusion"), # We can also try to implement a network based method for example this, which builds different similarity networks in the layers and then aggregates them
                tags$li("Machine learning based: Multiple kernel learning?") # We could aybe also explore machine learnign based based methods
              )
            )
          )
        )
      )
    })
  }



  # App definition -----------------------------------------------------------
  shinyApp(ui = movida_ui, server = movida_server)
}


# idea-> come input do tutto gia calcolato poi ho le pagine che ho definito

# quindi una per le analisi principale, dove metto una data table in cui poi posso selezionare il contrast e vedere le due tabelle
# un altra pagina dove vedo l'enrichment plot e alteri lot anche (da selezionare, magari poi posso togliere le cards)
# metabolomics
# un altra pagina dove ci sono le kegg maps, le scelgo dai pathway che sono enriched e quelli che non sono enriched
# poi posso aggiungere le pagine per l'integration, quindi pvalue based integration on the results both of DE and enrichment. these can also be network based where a network of related elements of different omics layers are shown
# clustering / early integration, these methods work by some by creating creating a some kind of feature space where the three data types are comparable and downstrea analysis can be done, this usually means clustering, either they are even more black blox, so machine learning based, multikernel

# output$gs_de_network_trans <- renderVisNetwork({

#     req(input$selected_contrast)  # Ensure the input is available

#     # Filter the table based on the selected contrast
#     enrichresults_table <- GeneTonic::shake_topGOtableResult(reactive_values$results_trans[[input$selected_contrast]]$topGO_tbl)

#     # Assuming enrichment_map() returns an igraph object
#     g <- enrichment_map(enrichresults_table)  # Call the function that returns the igraph object

#     # Convert the igraph object to visNetwork format
#     network_data <- toVisNetworkData(g)

#     # Generate the network visualization
#     visNetwork::visNetwork(nodes = network_data$nodes, edges = network_data$edges) %>%
#         visNetwork::visIgraphLayout()  # Layout for igraph visualization
# })

# output$kegg_map_plot <- ggiraph::renderGirafe({
#     library(dplyr)
#     library(tidygraph)
#     library(ggraph)
#     library(ggkegg)   # Or the appropriate pathway tool you're using

#     req(input$kegg_pathway_select)  # Ensure the user has selected a pathway


#     pathway_id <- input$kegg_pathway_select

#     p <- pathway(pathway_id) |>  ## Obtain and parse the KEGG pathway
#         activate(nodes) |>
#         mutate(
#             convert_hsa = convert_id("hsa"),
#             convert_map = convert_id("pathway")
#         ) |>
#         ggraph(x = x, y = y) +
#         geom_edge_parallel(
#             arrow = arrow(length = unit(1, "mm")),
#             aes(linetype = subtype_name),
#             end_cap = circle(7.5, "mm")
#         ) +
#         geom_node_rect(
#         aes(filter = type == "gene", fill = I(bgcolor)),
#             color = "black"
#         ) +
#         geom_node_text(
#         aes(label = convert_hsa),
#             size = 2, family = "serif"
#         ) +
#         geom_node_text(
#         aes(label = convert_hsa, filter = !is.na(convert_hsa) & convert_hsa == "TP53"),
#             size = 2, color = "red", family = "serif"
#         ) +
#         theme_void()

#         # Render as interactive ggiraph
#         ggiraph::girafe(ggobj = p,
#                 options = list(
#                 ggiraph::opts_hover(css = "fill:red;"),
#                 ggiraph::opts_sizing(rescale = TRUE),
#                 ggiraph::opts_selection(type = "multiple")
#                 ),
#             height_svg = 7,
#             width_svg = 10
#             )
#     })

# #### Plots panel ------------------------------------------------------
# output$ui_panel_plots <- renderUI({
#     tagList(
#         fluidRow(
#             bs4Dash::column()
#         )
#     )
# })

# visNetworkOutput("kegg_map_plot_vis", height = "600px")
# bs4Dash::bs4Card(
#     id = "tab_overview_card_trans",
#     title = "Kegg map",
#     elevation = 2,
#     width = 12,
#     closable = FALSE,

# )

# Reactive expression for KEGG pathways
# kegg_choices <- reactive({

#     # Set config to skip SSL verification
#     config <- config(ssl_verifypeer = 0L)
#     res <- httr::GET("http://rest.kegg.jp/list/pathway/hsa")
#     if (res$status_code != 200) return(NULL)

#     text <- content(res, as = "text")
#     lines <- strsplit(text, "\n")[[1]]
#     lines <- lines[nzchar(lines)]  # remove empty lines

#     pathways <- strsplit(lines, "\t")
#     ids <- sapply(pathways, function(x) gsub("path:", "", x[1]))
#     names <- sapply(pathways, function(x) x[2])

#     setNames(ids, names)
# })

# output$kegg_selector <- renderUI({
#     # req(kegg_choices())  # your reactive pathway list

#     selectInput(
#         inputId = "kegg_pathway_select",
#         label = "Choose KEGG Pathway",
#         choices = c("hsa01100","hsa01200","hsa01210","hsa01212","hsa01230","hsa01232","hsa01250","hsa01240","hsa01320","hsa00010","hsa00020","hsa00030","hsa00040","hsa00051","hsa00052","hsa00053","hsa00500","hsa00620","hsa00630","hsa00640","hsa00650","hsa00562","hsa00190","hsa00910","hsa00920","hsa00061","hsa00062","hsa00071","hsa00100","hsa00120","hsa00140","hsa00561","hsa00564","hsa00565","hsa00600","hsa00590","hsa00591","hsa00592","hsa01040","hsa00230"),
#         selected = "hsa04110",
#         width = "100%"
#     )
#     })
