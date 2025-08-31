library(bslib)

ui <- page_fluid(
  theme = app_theme,
  tags$head(tags$style(HTML('
    .card-like { background: #fff; border: 1px solid #e3e6ea; border-radius: 8px; padding: 12px; margin-bottom: 12px; }
    .value-box { background: #f8f9fb; border: 1px solid #e3e6ea; border-radius: 8px; padding: 8px 12px; text-align: center; }
    .vb-title { font-size: 12px; color: #6c757d; }
    .vb-value { font-size: 20px; font-weight: 700; }
    .vb-subtitle { font-size: 11px; color: #6c757d; }
    .muted { color: #6c757d; font-size: 12px; }

    /* Compact celltype cards */
    .celltype-card { cursor: pointer; transition: transform 0.1s ease, box-shadow 0.1s ease; border: 2px solid transparent; }
    .celltype-card:hover { transform: translateY(-1px); box-shadow: 0 2px 10px rgba(0,0,0,0.08); }
    .celltype-card.selected { border-color: #2E86AB; }
    .celltype-chip { display:inline-block; padding: 2px 6px; border-radius: 999px; font-size: 11px; color: #fff; margin-left: 6px; }
  '))),
  layout_sidebar(
    sidebar = sidebar(
      width = 320,
      h4("Controls"),
      div(class = "card-like",
        div(class = "sidebar-label", "Results directory"),
        selectInput(
          inputId = "result_dir",
          label = NULL,
          choices = list_result_dirs(),
          selected = character(0),
          multiple = FALSE
        ),
        div(style = "margin-top:8px;",
          actionButton("refresh", "Refresh list", class = "btn btn-primary")
        ),
        div(style = "margin-top:8px; font-size:12px; color:#6c757d;",
          HTML(paste0("Significance criteria: <code>padj < ", PADJ_THRESH, "</code> and <code>|log2FC| >= ", LFC_THRESH, "</code>"))
        ),
        if (HALLMARK_COUNT > 0) div(style = "margin-top:6px; font-size:12px; color:#6c757d;",
          HTML(paste0("GSEA tests ", HALLMARK_COUNT, " Hallmark pathways (MSigDB H). ",
                      "See ", '<a href="https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H" target="_blank">collection H</a>.'))
        )
      )
    ),
    div(
      # Step 1: Overview when result directory is selected but no cell type
      uiOutput("no_directory_panel"),
      uiOutput("overview_panel"),
      # Step 2: Detailed view once a cell type has been selected
      uiOutput("detailed_panel")
    )
  )
) 