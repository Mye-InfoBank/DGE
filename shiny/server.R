server <- function(input, output, session) {
  # Use reactive value for selected cell type instead of hidden input
  selected_celltype <- reactiveVal("")
  
  # Load pseudobulk metadata for selected celltype
  pseudobulk_metadata <- reactive({
    req(input$result_dir, selected_celltype(), selected_result_path())
    metadata_file <- file.path(selected_result_path(), selected_celltype(), "DESeq2_pseudobulk_metadata.txt")
    if (file.exists(metadata_file)) {
      data <- read.table(metadata_file, header = TRUE, sep = "\t", check.names = FALSE)
      return(data)
    }
    return(NULL)
  })
  
  celltype_params <- reactive({
    req(input$result_dir, selected_celltype(), selected_result_path())
    params_file <- file.path(selected_result_path(), selected_celltype(), "DESeq2_run_parameters.tsv")
    if (file.exists(params_file)) {
      data <- read.table(params_file, header = TRUE, sep = "\t", check.names = FALSE)
      return(data)
    }
    return(NULL)
  })

  # Get summary statistics for pseudobulk samples
  pseudobulk_stats_reactive <- reactive({
    req(pseudobulk_metadata(), input$contrast, celltype_params())
    meta <- pseudobulk_metadata()
    params_df <- celltype_params()
    
    # Split contrast name to get groups (assuming format "A_vs_B" or "A_vs_rest")
    groups <- strsplit(input$contrast, "_vs_")[[1]]
    group1 <- groups[1]
    group2 <- groups[2]
    
    # Initialize variables
    group1_samples <- NULL
    group2_samples <- NULL
    
    # For one-vs-all mode 
    if (params_df$one_vs_all) {
      condition_col <- params_df$category_column
      meta$tmp2 <- ifelse(meta[[condition_col]] == group1, group1, 'rest')
      group1_samples <- meta[meta$tmp2 == group1, ]
      group2_samples <- meta[meta$tmp2 == "rest", ]
      
    } else {
      # For pairwise mode - use the actual condition column from parameters
      condition_col <- params_df$category_column
      
      if (!is.null(condition_col) && condition_col %in% names(meta)) {
        group1_samples <- meta[meta[[condition_col]] == group1, ]
        group2_samples <- meta[meta[[condition_col]] == group2, ]
      }
    }
    
    # If we couldn't identify groups, return NULL
    if (is.null(group1_samples) || is.null(group2_samples) || 
        nrow(group1_samples) == 0 || nrow(group2_samples) == 0) {
      return(NULL)
    }
    
    list(
      group1_name = group1,
      group2_name = group2,
      group1_n_samples = nrow(group1_samples),
      group2_n_samples = nrow(group2_samples),
      group1_n_cells_per_sample = group1_samples$ncells,
      group2_n_cells_per_sample = group2_samples$ncells
    )
  })
  # Add a counter for cell type clicks to handle repeated clicks
  celltype_click_counter <- reactiveVal(0)
  
  # Add loading state for overview data
  overview_loading <- reactiveVal(FALSE)
  
  # Helper to format seconds as mm:ss
  format_eta <- function(seconds_total) {
    seconds_total <- as.integer(round(seconds_total))
    minutes <- seconds_total %/% 60
    seconds <- seconds_total %% 60
    if (minutes > 0) sprintf("%dm %02ds", minutes, seconds) else sprintf("%ds", seconds)
  }
  
  # Helper to sanitize pathway names to file names
  sanitize_pathway <- function(x) {
    toupper(gsub("[^A-Za-z0-9_]+", "_", x))
  }
  
  # Refresh available directories on demand
  observeEvent(input$refresh, {
    updateSelectInput(session, "result_dir", choices = list_result_dirs(), selected = NULL)
  })
  
  # Also refresh at startup, but don't select anything by default
  observe({
    updateSelectInput(session, "result_dir", choices = list_result_dirs(), selected = NULL)
  })
  
  # Reset selected cell type when result directory changes
  observeEvent(input$result_dir, {
    selected_celltype("")
    celltype_click_counter(0)
    overview_loading(FALSE)  # Reset loading state
  })
  
  # Back to overview button
  observeEvent(input$back_to_overview, {
    selected_celltype("")
    celltype_click_counter(0)
  })

  # Overview data - with loading state and progress bar
  overview_data <- reactive({
    if (is.null(input$result_dir) || !nzchar(input$result_dir)) return(list())
    
    overview_loading(TRUE)
    on.exit({ overview_loading(FALSE) }, add = TRUE)
    
    withProgress(message = "Loading results", value = 0, {
      # local stateful progress function
      last_val <- 0
      pfun <- function(value, detail = NULL) {
        delta <- max(0, min(1, value)) - last_val
        if (delta > 0) incProgress(amount = delta, detail = detail)
        last_val <<- max(0, min(1, value))
      }
      pfun(0)
      res <- scan_result_overview(input$result_dir, progress = pfun)
      pfun(1, "Done")
      res
    })
  })

 
  # Parse DESeq2 parameter settings from the start of the latest log
  deseq2_params <- reactive({
    p <- latest_log_paths()$deseq2
    if (is.na(p)) return(data.frame())
    lines <- tryCatch(readLines(p, warn = FALSE, n = 200), error = function(e) character())
    if (!length(lines)) return(data.frame())
    # Keep lines that look like key: value or key = value and are not comments or obvious runtime logs
    cand <- lines[grepl("^\\s*[^#\\[][^:=]{0,60}(\\:|=)\\s*.+$", lines)]
    # Stop at first empty line after we've started collecting
    if (length(cand) == 0) return(data.frame())
    # Normalize and split into key/value
    kv <- do.call(rbind, lapply(cand, function(x) {
      y <- sub("^\\s+|\\s+$", "", x)
      # split on first ':' or '='
      m <- regexpr(":|=", y)
      if (m[1] <= 0) return(c(NA, NA))
      key <- trimws(substr(y, 1, m[1]-1))
      val <- trimws(substr(y, m[1]+1, nchar(y)))
      c(key, val)
    }))
    if (is.null(kv) || nrow(kv) == 0) return(data.frame())
    df <- as.data.frame(kv, stringsAsFactors = FALSE)
    names(df) <- c("Parameter", "Value")
    # Deduplicate keeping first occurrence
    df <- df[!duplicated(df$Parameter) & nzchar(df$Parameter) & nzchar(df$Value), , drop = FALSE]
    head(df, 30)
  })
  
  # Overview panel
  output$overview_panel <- renderUI({
    if (is.null(input$result_dir) || !nzchar(input$result_dir) || nzchar(selected_celltype())) {
      return(NULL)
    }
    
    # Show loading state
    if (overview_loading()) {
      return(
        div(class = "card-like",
          h3("Loading Results Overview..."),
          div(style = "text-align: center; padding: 40px;",
            tags$p("Scanning directory and analyzing data..."),
            tags$div(class = "spinner-border text-primary", role = "status",
              tags$span(class = "visually-hidden", "Loading...")
            )
          )
        )
      )
    }
    
    tagList(
      div(class = "card-like",
        h3("Results Overview"),
        uiOutput("overview_path_info"),
        fluidRow(
          column(4, value_box("Cell Types", textOutput("vb_celltypes"))),
          column(4, value_box("Total Contrasts", textOutput("vb_total_contrasts"))),
          column(4, value_box("Total Sig. Genes", textOutput("vb_total_sig_genes")))
        )
      ),
      div(class = "card-like",
        h4("Select Cell Type for Detailed Analysis"),
        uiOutput("celltype_selection_cards")
      ),
      div(class = "card-like",
        h4("Overview Heatmap"),
        plotOutput("overview_heatmap", height = "600px") %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB")
      ),
      div(class = "card-like",
        h4("Cell Type Totals"),
        uiOutput("overview_summary_dynamic")
      ),
      div(class = "card-like",
        h4("Cell Type Summary"),
        DTOutput("overview_summary_table") %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB")
      ),
      div(class = "card-like",
        h4("DESeq2 Parameters (from latest log)"),
        uiOutput("deseq2_params_ui")
      )
    )
  })
  
  output$deseq2_params_ui <- renderUI({
    df <- deseq2_params()
    if (is.null(df) || nrow(df) == 0) {
      return(tags$em("No parameter information detected at the start of the latest DESeq2 log."))
    }
    # Render as compact table
    DT::datatable(
      df,
      options = list(pageLength = 10, dom = 'tip', ordering = FALSE),
      rownames = FALSE,
      class = 'compact stripe'
    )
  })
  
  # No directory selected panel
  output$no_directory_panel <- renderUI({
    if (!is.null(input$result_dir) && nzchar(input$result_dir)) {
      return(NULL)
    }
    
    div(class = "card-like",
      h3("Welcome to DEA Results Explorer"),
      div(style = "text-align: center; padding: 40px;",
        tags$p("Please select a results directory from the sidebar to begin exploring your differential expression analysis results."),
        tags$p(class = "muted", "The app will scan the selected directory and show you an overview of all available cell types and contrasts.")
      )
    )
  })
  
  # Dynamic height for rotated barplot
  output$overview_summary_dynamic <- renderUI({
    data <- overview_data()
    n_ct <- length(data)
    # 28px per category + margins, min 300px, cap at 1200px
    height_px <- max(300, min(1200, n_ct * 28 + 80))
    plotOutput("overview_summary", height = paste0(height_px, "px")) %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB")
  })
  
  output$overview_heatmap <- renderPlot({
    make_overview_heatmap(overview_data())
  }, res = 110)
  
  output$overview_summary <- renderPlot({
    make_overview_summary_plot(overview_data())
  }, res = 110)
  
  output$overview_summary_table <- renderDT({
    data <- overview_data()
    if (length(data) == 0) return(datatable(data.frame()))
    
    summary_df <- data.frame(
      celltype = names(data),
      contrasts = sapply(data, function(x) length(x$contrasts)),
      comparisons = sapply(data, function(x) length(x$comparisons)),
      total_sig_genes = sapply(data, function(x) x$total_sig_genes),
      total_sig_paths = sapply(data, function(x) x$total_sig_paths),
      stringsAsFactors = FALSE
    )
    
    datatable(
      summary_df,
      options = list(pageLength = 10, dom = 'tip'),
      rownames = FALSE,
      selection = 'single'
    )
  })
  
  # Cell type selection cards
  output$celltype_selection_cards <- renderUI({
    data <- overview_data()
    if (length(data) == 0) return(tags$em("No cell types found."))
    
    # Order by major class using palette order, then by name
    celltypes <- names(data)
    major_classes <- vapply(celltypes, function(ct) {
      mc <- data[[ct]]$major_class
      if (is.null(mc) || !nzchar(mc)) "Other" else mc
    }, FUN.VALUE = character(1))
    class_order <- if (exists("MAJOR_CLASS_PALETTE", inherits = TRUE)) names(MAJOR_CLASS_PALETTE) else unique(c(major_classes, "Other"))
    order_idx <- match(major_classes, class_order)
    order_idx[is.na(order_idx)] <- length(class_order) + 1
    ordered_celltypes <- celltypes[order(order_idx, tolower(celltypes))]

    # For tiny bar normalization
    max_total <- max(sapply(data, function(x) x$total_sig_genes), na.rm = TRUE)
    if (!is.finite(max_total)) max_total <- 0

    # Compact cards with subtle color tint; larger name, minimal whitespace
    cards <- lapply(ordered_celltypes, function(ct) {
      ct_data <- data[[ct]]
      bg <- if (!is.null(ct_data$color)) color_alpha(ct_data$color, 0.08) else "#ffffff"
      br <- if (!is.null(ct_data$color)) color_alpha(ct_data$color, 0.35) else "#e3e6ea"
      pct <- if (max_total > 0) max(0, min(1, ct_data$total_sig_genes / max_total)) else 0
      bar_w <- sprintf("%.1f%%", pct * 100)
      tags$div(
        class = "card-like celltype-card",
        style = sprintf("margin-bottom: 6px; margin-right: 6px; display: inline-block; width: 200px; vertical-align: top; background:%s; border-color:%s; padding:8px;", bg, br),
        onclick = sprintf("Shiny.setInputValue('celltype_click', '%s', {priority: 'event'})", ct),
        tags$h5(style = "margin: 0 0 6px 0; font-size: 14px;", ct),
        tags$div(style = "font-size: 11px; line-height: 1.25; margin:0;",
          tags$div(tags$b("Contrasts:"), length(ct_data$contrasts)),
          tags$div(tags$b("Sig. genes:"), format(ct_data$total_sig_genes, big.mark = ","))
        ),
        tags$div(style = "margin-top:6px; background:#f1f3f5; height:6px; border-radius:4px; overflow:hidden;",
          tags$div(style = sprintf("width:%s; height:100%%; background:%s;", bar_w, br))
        )
      )
    })

    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 6px;",
      cards
    )
  })
  
  # Cell type selection handler - prioritize event; allow repeated same-value clicks
  observeEvent(input$celltype_click, {
    val <- input$celltype_click
    if (!is.null(val) && nzchar(val)) selected_celltype(val)
  }, ignoreInit = TRUE)
  
  # Selected paths for detailed view
  selected_result_path <- reactive({ resolve_selected_dir(input$result_dir) })
  selected_celltype_path <- reactive({ 
    ct <- selected_celltype()
    if (is.null(ct) || !nzchar(ct)) return(NA_character_)
    resolve_celltype_dir(input$result_dir, ct) 
  })
  
  output$selected_paths <- renderUI({
    res <- selected_result_path()
    ct <- selected_celltype_path()
    tagList(
      tags$p(tags$b("Results dir:"), ifelse(is.na(res), "—", res)),
      tags$p(tags$b("Cell type dir:"), ifelse(is.na(ct), "—", ct))
    )
  })
  
  # Discover files for selected cell type
  celltype_files <- reactive({ find_celltype_files(selected_celltype_path()) })
  
  # DEG
  deg_df <- reactive({ read_deg_table(celltype_files()$deg_tsv) })
  
  output$deg_info <- renderUI({
    files <- celltype_files()
    if (length(files$deg_tsv) == 0) return(tags$em("No DESeq2 TSV found in this cell type directory."))
    tags$p("Using file:", tags$code(files$deg_tsv[1]))
  })
  
  # Contrast selector based on available contrasts in the DE table (GLOBAL for all panels)
  output$contrast_selector <- renderUI({
    df <- deg_df()
    if (is.null(df) || !"contrast" %in% names(df)) return(NULL)
    contrasts <- sort(unique(df$contrast))
    if (length(contrasts) == 0) return(NULL)
    
    div(class = "card-like",
      h4("Select Contrast"),
      selectInput("contrast", "Choose contrast for analysis", choices = contrasts, selected = contrasts[1])
    )
  })

  # Render pseudobulk statistics
  output$pseudobulk_stats <- renderUI({
    stats <- pseudobulk_stats_reactive()
    if (is.null(stats)) return(NULL)
    
    # Create condition info text
    condition_info <- if (!is.null(stats$condition_column_used) && stats$condition_column_used != "Unknown") {
      div(style = "margin-bottom: 10px; font-size: 12px; color: #6c757d;",
        paste("Based on condition column:", stats$condition_column_used))
    } else {
      div(style = "margin-bottom: 10px; font-size: 12px; color: #6c757d;",
        "One-vs-all comparison mode")
    }
    
    div(
      class = "card-like pseudobulk-stats",
      h4("Sample Information"),
      condition_info,
      div(
        style = "display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 20px;",
        
        # Group 1 stats (left)
        div(
          class = "pseudobulk-box",
          h5(stats$group1_name, style = "color: #2E86AB;"),
          div(
            div(class = "pseudobulk-stat",
              div(class="pseudobulk-label", "Number of Samples"),
              div(class="pseudobulk-value", stats$group1_n_samples)
            ),
            div(class = "pseudobulk-stat",
              div(class="pseudobulk-label", "Total Cells"),
              div(class="pseudobulk-value", format(sum(stats$group1_n_cells_per_sample), big.mark=","))
            ),
            div(class = "pseudobulk-stat",
              div(class="pseudobulk-label", "Mean Cells per Sample"),
              div(class="pseudobulk-value", format(round(mean(stats$group1_n_cells_per_sample)), big.mark=","))
            ),
            div(class = "pseudobulk-stat",
              div(class="pseudobulk-label", "Median Cells per Sample"),
              div(class="pseudobulk-value", format(round(median(stats$group1_n_cells_per_sample)), big.mark=","))
            )
          )
        ),
        
        # Density plot (center)
        div(
          style = "display: flex; flex-direction: column; justify-content: center;",
          h5("Cell Count Distribution", style = "text-align: center; margin-bottom: 10px;"),
          plotOutput("pseudobulk_density_plot", height = "300px")
        ),
        
        # Group 2 stats (right)
        div(
          class = "pseudobulk-box",
          h5(stats$group2_name, style = "color: #A23B72;"),
          div(
            div(class = "pseudobulk-stat",
              div(class="pseudobulk-label", "Number of Samples"),
              div(class="pseudobulk-value", stats$group2_n_samples)
            ),
            div(class = "pseudobulk-stat",
                div(class="pseudobulk-label", "Total Cells"),
                div(class="pseudobulk-value", format(sum(stats$group2_n_cells_per_sample), big.mark=","))
            ),
            div(class = "pseudobulk-stat",
                div(class="pseudobulk-label", "Mean Cells per Sample"),
                div(class="pseudobulk-value", format(round(mean(stats$group2_n_cells_per_sample)), big.mark=","))
            ),
            div(class = "pseudobulk-stat",
                div(class="pseudobulk-label", "Median Cells per Sample"),
                div(class="pseudobulk-value", format(round(median(stats$group2_n_cells_per_sample)), big.mark=","))
            )
          )
        )
      )
    )
  })

  filtered_deg_df <- reactive({
    df <- deg_df()
    if (is.null(df)) return(NULL)
    if ("contrast" %in% names(df) && !is.null(input$contrast) && nzchar(input$contrast)) {
      df <- subset(df, contrast == input$contrast)
    }
    df
  })
  
  # Populate gene choices for highlight selectize
  observe({
    df <- filtered_deg_df()
    genes <- if (!is.null(df) && "gene" %in% names(df)) sort(unique(df$gene)) else character()
    updateSelectizeInput(session, "highlight_genes", choices = genes, server = TRUE)
  })
  
  # Overview value boxes using analysis metrics (matching thresholds)
  observe({
    df_deg <- deg_df()
    df_gsea <- gsea_df()
    n_contrasts <- if (!is.null(df_deg) && "contrast" %in% names(df_deg)) length(unique(df_deg$contrast)) else 0
    sig_genes <- if (!is.null(df_deg) && all(c("padj", "log2foldchange") %in% names(df_deg))) sum(df_deg$padj < PADJ_THRESH & abs(df_deg$log2foldchange) >= LFC_THRESH, na.rm = TRUE) else 0
    sig_paths <- if (!is.null(df_gsea) && "fdr" %in% names(df_gsea)) sum(df_gsea$fdr <= 0.25, na.rm = TRUE) else 0
    output$vb_contrasts <- renderText(format(n_contrasts, big.mark = ","))
    output$vb_sig_genes <- renderText(format(sig_genes, big.mark = ","))
    output$vb_sig_paths <- renderText(format(sig_paths, big.mark = ","))
  })
  
  # Density plot for pseudobulk cell distribution
  output$pseudobulk_density_plot <- renderPlot({
    stats <- pseudobulk_stats_reactive()
    if (is.null(stats)) return(NULL)
    
    # Create data frame for plotting
    plot_data <- data.frame(
      group = c(rep(stats$group1_name, length(stats$group1_n_cells_per_sample)),
                rep(stats$group2_name, length(stats$group2_n_cells_per_sample))),
      cells_per_sample = c(stats$group1_n_cells_per_sample, stats$group2_n_cells_per_sample)
    )
    
    # Create density plot
    library(ggplot2)
    ggplot(plot_data, aes(x = cells_per_sample, fill = group)) +
      geom_density(alpha = 0.7) +
      scale_fill_manual(values = c("#2E86AB", "#A23B72")) +
      labs(
        x = "Cells per Sample",
        y = "Density",
        fill = "Group"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10)
      ) +
      scale_x_continuous(labels = scales::comma_format())
  }, res = 110)
  
  output$volcano_plot <- renderPlot({
    hl <- input$highlight_genes
    hl_col <- if (!is.null(input$highlight_color_name) && nzchar(input$highlight_color_name)) input$highlight_color_name else "#FFD166"
    sig_col <- if (!is.null(input$signif_color_name) && nzchar(input$signif_color_name)) input$signif_color_name else "#E4572E"
    topn <- if (!is.null(input$volcano_label_topn)) input$volcano_label_topn else 10
    make_volcano(filtered_deg_df(), highlight_genes = hl, highlight_color = hl_col, signif_color = sig_col, label_top_n = topn)
  }, res = 110)
  
  output$deg_table <- renderDT({
    df <- filtered_deg_df()
    if (is.null(df)) return(datatable(data.frame()))
    # Order by adjusted p-value ascending if present
    if ("padj" %in% names(df)) {
      ord <- order(df$padj, na.last = TRUE)
      df <- df[ord, , drop = FALSE]
    }
    datatable(df, options = list(pageLength = 15, scrollX = TRUE), filter = "top", rownames = FALSE)
  })
  
  # GSEA
  gsea_df <- reactive({ read_gsea_results(celltype_files()$gsea_csvs) })
  
  output$gsea_info <- renderUI({
    files <- celltype_files()
    if (length(files$gsea_csvs) == 0) return(tags$em("No GSEA CSVs found in this cell type directory."))
    links_note <- if (HALLMARK_COUNT > 0) paste0(" (", HALLMARK_COUNT, " Hallmark pathways tested)") else ""
    tags$p("Merged ", length(files$gsea_csvs), " GSEA files", links_note, ".")
  })
  
  # GSEA contrast selector using the 'comparison' field
  output$gsea_contrast_selector <- renderUI({
    df <- gsea_df()
    if (is.null(df) || !"comparison" %in% names(df)) return(NULL)
    comps <- sort(unique(df$comparison))
    if (length(comps) == 0) return(NULL)
    selectInput("gsea_contrast", "Select comparison", choices = comps, selected = comps[1])
  })
  
  filtered_gsea_df <- reactive({
    df <- gsea_df()
    if (is.null(df)) return(NULL)
    if ("comparison" %in% names(df) && !is.null(input$gsea_contrast) && nzchar(input$gsea_contrast)) {
      df <- subset(df, comparison == input$gsea_contrast)
    }
    df
  })
  
  output$gsea_doc <- renderUI({
    ct <- selected_celltype()
    comp <- input$gsea_contrast
    if (is.null(ct) || !nzchar(ct) || is.null(comp) || !nzchar(comp)) return(NULL)
    tags$p(class = "muted",
      paste0("For ", tags$b(ct), " in contrast ", tags$b(comp), ": a positive normalized enrichment score (NES) means the pathway's genes are more up-weighted in the first group of the contrast (left side of 'A_vs_B'), while a negative NES means enrichment in the second group (right side). Magnitude reflects the strength of enrichment; FDR controls significance.")
    )
  })
  
  # GSEA enrichment PNG selection and display
  output$gsea_pathway_selector <- renderUI({
    df <- filtered_gsea_df()
    if (is.null(df) || !"pathway" %in% names(df)) return(NULL)
    # default: top by abs NES
    if ("nes" %in% names(df)) {
      df <- df[order(-abs(df$nes)), , drop = FALSE]
    }
    choices <- unique(df$pathway)
    selectInput("gsea_pathway", "Select pathway", choices = choices, selected = choices[1])
  })
  
  output$gsea_png <- renderUI({
    comp <- input$gsea_contrast
    pwy <- input$gsea_pathway
    ct_dir <- selected_celltype_path()
    if (is.null(comp) || !nzchar(comp) || is.null(pwy) || !nzchar(pwy) || is.na(ct_dir)) return(tags$em("Select a comparison and pathway to view the enrichment plot."))
    fname <- paste0(sanitize_pathway(pwy), ".png")
    fpath <- file.path(ct_dir, "gsea_plots", comp, fname)
    if (!file.exists(fpath)) return(tags$em("No enrichment PNG found for this pathway."))
    tags$img(src = base64enc::dataURI(file = fpath, mime = "image/png"), style = "max-width:100%; height:auto;")
  })
  
  output$gsea_plot <- renderPlot({
    topn <- if (!is.null(input$gsea_top_n)) input$gsea_top_n else 10
    make_gsea_plot(filtered_gsea_df(), top_n = topn)
  }, res = 110)
  
  output$gsea_table <- renderDT({
    df <- filtered_gsea_df()
    if (is.null(df)) return(datatable(data.frame()))

    # Add clickable pathway links if available
    if (!is.null(df$msigdb_url)) {
      df$pathway <- ifelse(!is.na(df$msigdb_url) & nzchar(df$msigdb_url),
                           paste0('<a href="', df$msigdb_url, '" target="_blank">', df$pathway, '</a>'),
                           df$pathway)
    }

    # Reorder columns: show key columns first
    preferred <- c("comparison", "pathway", "nes", "fdr")
    existing_pref <- intersect(preferred, colnames(df))
    rest <- setdiff(colnames(df), existing_pref)
    df <- df[, c(existing_pref, rest), drop = FALSE]

    # Identify wide/low-value columns to hide initially
    hide_cols <- intersect(c("source_file", "msigdb_url"), colnames(df))

    datatable(
      df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollY = '60vh',
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        dom = 'Btip',
        buttons = list('colvis', 'csv'),
        columnDefs = list(
          list(className = 'dt-nowrap', targets = seq_along(colnames(df)) - 1),
          list(visible = FALSE, targets = match(hide_cols, colnames(df)) - 1)
        )
      ),
      filter = "top",
      rownames = FALSE,
      escape = FALSE,
      class = 'nowrap stripe compact'
    )
  })
  
  # Summaries for Overview tab
  deg_summary <- reactive({
    df <- deg_df()
    if (is.null(df) || !all(c("contrast", "padj", "log2foldchange") %in% names(df))) return(data.frame())
    split_list <- split(df, df$contrast)
    res <- lapply(names(split_list), function(ct) {
      d <- split_list[[ct]]
      n_total <- nrow(d)
      n_sig <- sum(d$padj < PADJ_THRESH & abs(d$log2foldchange) >= LFC_THRESH, na.rm = TRUE)
      n_up <- sum(d$padj < PADJ_THRESH & d$log2foldchange >= LFC_THRESH, na.rm = TRUE)
      n_down <- sum(d$padj < PADJ_THRESH & d$log2foldchange <= -LFC_THRESH, na.rm = TRUE)
      data.frame(contrast = ct, genes_tested = n_total, sig_genes = n_sig, up = n_up, down = n_down, stringsAsFactors = FALSE)
    })
    do.call(rbind, res)
  })
  
  gsea_summary <- reactive({
    df <- gsea_df()
    if (is.null(df) || !"comparison" %in% names(df)) return(data.frame())
    split_list <- split(df, df$comparison)
    res <- lapply(names(split_list), function(cp) {
      d <- split_list[[cp]]
      n_total <- if (HALLMARK_COUNT > 0) HALLMARK_COUNT else nrow(d)
      n_sig <- if ("fdr" %in% names(d)) sum(d$fdr <= 0.25, na.rm = TRUE) else NA_integer_
      data.frame(comparison = cp, pathways_tested = n_total, sig_pathways = n_sig, stringsAsFactors = FALSE)
    })
    do.call(rbind, res)
  })
  
  output$deg_summary_table <- renderDT({
    df <- deg_summary()
    datatable(df, options = list(pageLength = 10, dom = 'tip'), rownames = FALSE)
  })
  
  output$gsea_summary_table <- renderDT({
    df <- gsea_summary()
    datatable(df, options = list(pageLength = 10, dom = 'tip'), rownames = FALSE)
  })
  
  # Advanced logs: find latest DESeq2 and GSEA log files under selected result dir
  latest_log_paths <- reactive({
    base <- selected_result_path()
    if (is.na(base)) return(list(deseq2 = NA_character_, gsea = NA_character_))
    # look for logs at the root of result directory (matching your pipelines)
    logs <- list.files(base, pattern = "(deseq2_analysis_.*\\.log|gsea_.*\\.log)$", full.names = TRUE, recursive = TRUE)
    if (length(logs) == 0) return(list(deseq2 = NA_character_, gsea = NA_character_))
    # pick latest by mtime separately
    finfo <- file.info(logs)
    deseq_idx <- grep("deseq2_analysis_.*\\.log$", basename(logs), ignore.case = TRUE)
    gsea_idx  <- grep("^gsea_.*\\.log$", basename(logs), ignore.case = TRUE)
    latest_deseq <- if (length(deseq_idx)) logs[deseq_idx[which.max(finfo$mtime[deseq_idx])]] else NA_character_
    latest_gsea  <- if (length(gsea_idx))  logs[gsea_idx[which.max(finfo$mtime[gsea_idx])]]   else NA_character_
    list(deseq2 = latest_deseq, gsea = latest_gsea)
  })
  
  # KPI and path outputs for overview (re-added)
  output$overview_path_info <- renderUI({
    res <- resolve_selected_dir(input$result_dir)
    tags$p(tags$b("Results dir:"), ifelse(is.na(res), "—", res))
  })
  
  output$vb_celltypes <- renderText({
    length(overview_data())
  })
  
  output$vb_total_contrasts <- renderText({
    data <- overview_data()
    if (length(data) == 0) return("0")
    all_contrasts <- unique(unlist(lapply(data, function(x) x$contrasts)))
    length(all_contrasts)
  })
  
  output$vb_total_sig_genes <- renderText({
    data <- overview_data()
    if (length(data) == 0) return("0")
    total <- sum(sapply(data, function(x) x$total_sig_genes))
    format(total, big.mark = ",")
  })
  
  output$deseq2_log_meta <- renderUI({
    p <- latest_log_paths()$deseq2
    if (is.na(p)) return(tags$em("No DESeq2 log found."))
    tags$p(tags$code(p))
  })
  
  output$gsea_log_meta <- renderUI({
    p <- latest_log_paths()$gsea
    if (is.na(p)) return(tags$em("No GSEA log found."))
    tags$p(tags$code(p))
  })
  
  output$deseq2_log <- renderText({
    p <- latest_log_paths()$deseq2
    if (is.na(p)) return("")
    paste(readLines(p, warn = FALSE), collapse = "\n")
  })
  
  output$gsea_log <- renderText({
    p <- latest_log_paths()$gsea
    if (is.na(p)) return("")
    paste(readLines(p, warn = FALSE), collapse = "\n")
  })
  
  # PNG preview
  output$volcano_png_ui <- renderUI({
    files <- celltype_files()
    if (length(files$volcano_png) == 0) return(tags$em("No volcano.png found."))
    tags$img(src = base64enc::dataURI(file = files$volcano_png[1], mime = "image/png"), style = "max-width:100%; height:auto;")
  })
  
  # Detailed panel
  output$detailed_panel <- renderUI({
    if (is.null(input$result_dir) || !nzchar(input$result_dir) || !nzchar(selected_celltype())) {
      return(NULL)
    }
    
    tagList(
      div(class = "card-like",
        uiOutput("selected_paths"),
        fluidRow(
          column(4, value_box("Contrasts", textOutput("vb_contrasts"))),
          column(4, value_box("Sig. genes", textOutput("vb_sig_genes"), subtitle = paste0("padj < ", PADJ_THRESH, ", |log2FC| >= ", LFC_THRESH))),
          column(4, value_box("Sig. pathways", textOutput("vb_sig_paths"), subtitle = "FDR <= 0.25"))
        ),
        div(style = "margin-top: 12px;",
          actionButton("back_to_overview", "← Back to Overview", class = "btn btn-outline-secondary")
        )
      ),
      # Global contrast selector for all panels
      uiOutput("contrast_selector"),
      # Pseudobulk stats (uses the global contrast selection)
      uiOutput("pseudobulk_stats"),
      tabsetPanel(type = "pills",
        tabPanel(
          title = "Overview",
          div(class = "card-like",
            h4("DE contrasts"),
            DTOutput("deg_summary_table") %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB")
          ),
          div(class = "card-like",
            h4("GSEA comparisons"),
            DTOutput("gsea_summary_table") %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB")
          )
        ),
        tabPanel(
          title = "DESeq2",
          div(class = "card-like",
            h4("Differential expression"),
            uiOutput("deg_info"),
            fluidRow(
              column(6,
                selectizeInput(
                  "highlight_genes",
                  label = "Highlight genes (type to search; press Enter to add)",
                  choices = NULL,
                  multiple = TRUE,
                  options = list(placeholder = 'Start typing gene symbols...', create = TRUE)
                )
              ),
              column(3,
                selectInput("highlight_color_name", "Highlight color", choices = VOLCANO_COLORS, selected = "#FFD166")
              ),
              column(3,
                selectInput("signif_color_name", "Significant color", choices = VOLCANO_COLORS, selected = "#E4572E")
              )
            ),
            sliderInput("volcano_label_topn", "Label top N significant genes", min = 0, max = 50, value = 10, step = 1),
            plotOutput("volcano_plot", height = "640px") %>% shinycssloaders::withSpinner(type = 7, color = "#E4572E")
          ),
          div(class = "card-like",
            h4("DE results table"),
            DTOutput("deg_table") %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB")
          )
        ),
        tabPanel(
          title = "GSEA",
          div(class = "card-like",
            h4("Enrichment (Hallmark)"),
            uiOutput("gsea_info"),
            uiOutput("gsea_contrast_selector"),
            sliderInput("gsea_top_n", "Top pathways to show", min = 5, max = 100, value = 10, step = 5),
            uiOutput("gsea_doc"),
            plotOutput("gsea_plot", height = "520px") %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB"),
            uiOutput("gsea_pathway_selector"),
            div(style = "margin-top: 12px;",
              h5("Enrichment plot (PNG)"),
              uiOutput("gsea_png")
            )
          ),
          div(class = "card-like",
            h4("GSEA table (merged)"),
            div(style = "overflow-x:auto;", DTOutput("gsea_table") %>% shinycssloaders::withSpinner(type = 7, color = "#2E86AB"))
          )
        ),
        tabPanel(
          title = "PNG Preview",
          div(class = "card-like",
            h4("Volcano (static PNG)"),
            uiOutput("volcano_png_ui")
          )
        ),
        tabPanel(
          title = "Logs (advanced)",
          div(class = "card-like",
            span(class = "muted", "Technical details for troubleshooting; not needed for normal use."),
            h4("Latest DESeq2 log"),
            uiOutput("deseq2_log_meta"),
            div(class = "logbox", verbatimTextOutput("deseq2_log")),
            h4(style = "margin-top:16px;", "Latest GSEA log"),
            uiOutput("gsea_log_meta"),
            div(class = "logbox", verbatimTextOutput("gsea_log"))
          )
        )
      )
    )
  })
} 