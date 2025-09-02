suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(ggplot2)
  library(DT)
  library(bslib)
  library(shinycssloaders)
  library(msigdbr)
  library(ggrepel)
  library(base64enc)
  library(future)
})

# Absolute path to the DEA output root
OUTPUT_ROOT <- "/nfs/data/COST_IBD/downstream_tasks/dea/output"

# Significance thresholds (aligned with deseq2.R)
PADJ_THRESH <- 0.05
LFC_THRESH  <- 1.0

# MSigDB Hallmark terms/count (constant across comparisons)
.HALLMARK_DF <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) NULL)
HALLMARK_TERMS <- if (!is.null(.HALLMARK_DF)) unique(.HALLMARK_DF$gs_name) else character()
HALLMARK_COUNT <- length(HALLMARK_TERMS)

# Theme
app_theme <- bs_theme(version = 5, bootswatch = "flatly")

# Fixed 10-color palette for UI choices
VOLCANO_COLORS <- c(
  "blue"="#1f77b4", "orange"="#ff7f0e", "green"="#2ca02c", "red"="#d62728", "purple"="#9467bd",
  "brown"="#8c564b", "pink"="#e377c2", "grey"="#7f7f7f", "yellow"="#bcbd22", "cyan"="#17becf"
)

# Major cell type mapping and palette
MAJOR_CLASS_PALETTE <- c(
  "T/NK/ILC" = "#1f77b4",
  "B/Plasma" = "#2ca02c",
  "Myeloid" = "#d62728",
  "Epithelial" = "#ff7f0e",
  "Endothelial" = "#9467bd",
  "Stromal" = "#8c564b",
  "Mast" = "#e377c2",
  "Neuronal/Glia" = "#17becf",
  "Other" = "#7f7f7f"
)

# Convert hex color to rgba string with alpha (0..1)
color_alpha <- function(hex, alpha = 1) {
  h <- gsub("#", "", hex)
  if (nchar(h) == 3) {
    h <- paste0(substr(h, 1, 1), substr(h, 1, 1), substr(h, 2, 2), substr(h, 2, 2), substr(h, 3, 3), substr(h, 3, 3))
  }
  r <- strtoi(substr(h, 1, 2), base = 16)
  g <- strtoi(substr(h, 3, 4), base = 16)
  b <- strtoi(substr(h, 5, 6), base = 16)
  sprintf("rgba(%d,%d,%d,%.3f)", r, g, b, max(0, min(1, alpha)))
}

major_celltype_of <- function(name) {
  n <- tolower(name)
  if (grepl("t cell|cd4|cd8|treg|nk|ilc", n)) return("T/NK/ILC")
  if (grepl("b cell|bcell|plasma|plasmablast", n)) return("B/Plasma")
  if (grepl("mono|macro|myeloid|neutro|dendritic|dc|mast", n)) {
    if (grepl("mast", n)) return("Mast")
    return("Myeloid")
  }
  if (grepl("epi|epith|enterocyte|goblet|tuft|paneth|hepat|keratino|basal|club|alveol|cilia", n)) return("Epithelial")
  if (grepl("endo|vascular", n)) return("Endothelial")
  if (grepl("fibro|stromal|myofibro|mesench|smooth|pericyte|stellate", n)) return("Stromal")
  if (grepl("neur|glia|schwann", n)) return("Neuronal/Glia")
  "Other"
}

color_for_major_celltype <- function(name) {
  cls <- major_celltype_of(name)
  col <- MAJOR_CLASS_PALETTE[[cls]]
  if (is.null(col) || length(col) == 0) MAJOR_CLASS_PALETTE[["Other"]] else col
}

# ---------------- Helpers ---------------- #
list_result_dirs <- function() {
  all_entries <- list.files(OUTPUT_ROOT, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  dirs <- all_entries[file.info(all_entries)$isdir]
  c('None', basename(dirs))
}

list_celltype_dirs <- function(result_basename) {
  if (is.null(result_basename) || !nzchar(result_basename)) return(character())
  base_dir <- file.path(OUTPUT_ROOT, result_basename)
  if (!dir.exists(base_dir)) return(character())
  entries <- list.files(base_dir, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  dirs <- entries[file.info(entries)$isdir]
  basename(dirs)
}

resolve_selected_dir <- function(selected_basename) {
  if (is.null(selected_basename) || is.na(selected_basename) || !nzchar(selected_basename)) return(NA_character_)
  candidate <- file.path(OUTPUT_ROOT, selected_basename)
  if (!dir.exists(candidate)) return(NA_character_)
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

resolve_celltype_dir <- function(result_basename, celltype_basename) {
  if (is.null(result_basename) || is.null(celltype_basename)) return(NA_character_)
  if (!nzchar(result_basename) || !nzchar(celltype_basename)) return(NA_character_)
  candidate <- file.path(OUTPUT_ROOT, result_basename, celltype_basename)
  if (!dir.exists(candidate)) return(NA_character_)
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

summarize_directory <- function(dir_path) {
  if (is.na(dir_path)) return(data.frame(metric = character(), value = character(), stringsAsFactors = FALSE))
  files <- list.files(dir_path, full.names = TRUE, recursive = TRUE, include.dirs = FALSE)
  if (length(files) == 0) {
    return(data.frame(metric = c("Total files", "Total size (MB)", "Subdirectories"), value = c(0, 0, 0), stringsAsFactors = FALSE))
  }
  finfo <- file.info(files)
  total_files <- nrow(finfo)
  total_size_mb <- round(sum(finfo$size, na.rm = TRUE) / (1024^2), 2)
  subdirs <- list.dirs(dir_path, recursive = TRUE, full.names = TRUE)
  subdirs <- setdiff(subdirs, dir_path)
  num_subdirs <- length(subdirs)
  data.frame(
    metric = c("Total files", "Total size (MB)", "Subdirectories"),
    value = c(total_files, total_size_mb, num_subdirs),
    stringsAsFactors = FALSE
  )
}

list_directory_files <- function(dir_path) {
  if (is.na(dir_path)) return(data.frame())
  files <- list.files(dir_path, full.names = TRUE, recursive = TRUE, include.dirs = FALSE)
  if (length(files) == 0) return(data.frame())
  finfo <- file.info(files)
  ext <- tools::file_ext(files)
  data.frame(
    path = files,
    name = basename(files),
    extension = ifelse(nzchar(ext), tolower(ext), ""),
    size_MB = round(finfo$size / (1024^2), 3),
    modified = finfo$mtime,
    stringsAsFactors = FALSE
  )
}

find_celltype_files <- function(celltype_dir) {
  if (is.na(celltype_dir)) return(list())
  volcano_png <- list.files(celltype_dir, pattern = "^volcano\\.png$", full.names = TRUE, ignore.case = TRUE)
  deg_tsv <- list.files(celltype_dir, pattern = "DESeq2.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)
  gsea_csvs <- list.files(celltype_dir, pattern = "gsea_results\\.csv$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  list(volcano_png = volcano_png, deg_tsv = deg_tsv, gsea_csvs = gsea_csvs)
}

read_deg_table <- function(path, max_rows = NULL) {
  if (length(path) == 0) return(NULL)
  existing <- path[file.exists(path)]
  if (length(existing) == 0) return(NULL)
  pieces <- lapply(existing, function(file) {
    args <- list(file = file, sep = "\t", header = TRUE, quote = "\"", na.strings = c("NA", ""), data.table = FALSE)
    if (!is.null(max_rows)) args$nrows <- max_rows
    dt <- tryCatch(do.call(fread, c(args, list(showProgress = FALSE))), error = function(e) NULL)
    if (is.null(dt)) return(NULL)
    df <- as.data.frame(dt)
    names(df) <- tolower(names(df))
    if (!"log2foldchange" %in% names(df)) {
      alt <- intersect(c("log2fc", "logfc", "log2_fold_change"), names(df))
      if (length(alt) > 0) names(df)[match(alt[1], names(df))] <- "log2foldchange"
    }
    if (!"padj" %in% names(df)) {
      alt <- intersect(c("fdr", "qvalue", "qval", "padjusted"), names(df))
      if (length(alt) > 0) names(df)[match(alt[1], names(df))] <- "padj"
    }
    if (!"pvalue" %in% names(df)) {
      alt <- intersect(c("pval", "p_value"), names(df))
      if (length(alt) > 0) names(df)[match(alt[1], names(df))] <- "pvalue"
    }
    if (!"gene" %in% names(df)) {
      alt <- intersect(c("symbol", "gene_symbol", "genes"), names(df))
      if (length(alt) > 0) names(df)[match(alt[1], names(df))] <- "gene"
    }
    # annotate source file to help disambiguate
    df$source_file <- basename(file)
    df
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) return(NULL)
  # align columns across pieces
  all_cols <- unique(unlist(lapply(pieces, names)))
  pieces <- lapply(pieces, function(d) {
    missing <- setdiff(all_cols, names(d))
    for (m in missing) d[[m]] <- NA
    d[, all_cols]
  })
  do.call(rbind, pieces)
}

make_volcano <- function(df, max_points = 50000, highlight_genes = NULL, highlight_color = "#FFD166", signif_color = "#E4572E", label_top_n = 10) {
  if (is.null(df)) return(NULL)
  req_cols <- c("log2foldchange", "padj")
  if (!all(req_cols %in% names(df))) return(NULL)
  plot_df <- df
  plot_df <- subset(plot_df, is.finite(log2foldchange) & is.finite(padj))
  plot_df$neglog10padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  plot_df$significant <- !is.na(plot_df$padj) & plot_df$padj < PADJ_THRESH & abs(plot_df$log2foldchange) >= LFC_THRESH
  if (!is.null(highlight_genes) && length(highlight_genes) > 0 && "gene" %in% names(plot_df)) {
    hl <- tolower(trimws(highlight_genes))
    plot_df$highlight <- tolower(plot_df$gene) %in% hl
  } else {
    plot_df$highlight <- FALSE
  }
  if (nrow(plot_df) > max_points) {
    set.seed(1)
    idx <- sample.int(nrow(plot_df), max_points)
    plot_df <- plot_df[idx, , drop = FALSE]
  }
  p <- ggplot(plot_df, aes(x = log2foldchange, y = neglog10padj, color = significant)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("FALSE" = "#9AA0A6", "TRUE" = signif_color)) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", color = "#B0BEC5") +
    geom_hline(yintercept = -log10(PADJ_THRESH), linetype = "dashed", color = "#B0BEC5") +
    labs(x = "log2 fold change", y = "-log10(padj)", color = paste0("sig: padj < ", PADJ_THRESH, " & |log2FC| >= ", LFC_THRESH)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  if (any(plot_df$highlight, na.rm = TRUE) && "gene" %in% names(plot_df)) {
    hl_df <- subset(plot_df, highlight)
    p <- p +
      geom_point(data = hl_df, aes(x = log2foldchange, y = neglog10padj), inherit.aes = FALSE,
                 color = highlight_color, size = 2.5, alpha = 0.95) +
      ggrepel::geom_text_repel(data = hl_df, aes(x = log2foldchange, y = neglog10padj, label = gene),
                               inherit.aes = FALSE, color = highlight_color, size = 3, max.overlaps = 100,
                               box.padding = 0.25, point.padding = 0.2, min.segment.length = 0)
  }
  # Label top-N significant genes by lowest padj (then highest |lfc|)
  if (label_top_n > 0 && "gene" %in% names(plot_df)) {
    sig_df <- subset(plot_df, significant & !is.na(padj))
    if (nrow(sig_df) > 0) {
      ord <- order(sig_df$padj, -abs(sig_df$log2foldchange))
      lab_df <- sig_df[head(ord, label_top_n), , drop = FALSE]
      p <- p + ggrepel::geom_text_repel(
        data = lab_df,
        aes(x = log2foldchange, y = neglog10padj, label = gene),
        size = 3, color = signif_color, max.overlaps = 100,
        box.padding = 0.25, point.padding = 0.2, min.segment.length = 0
      )
    }
  }
  p
}

msigdb_url <- function(id) {
  if (is.null(id) || is.na(id) || !nzchar(id)) return(NA_character_)
  paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/", toupper(gsub("[^A-Za-z0-9_]+", "_", id)), ".html")
}

read_gsea_results <- function(csv_paths) {
  if (length(csv_paths) == 0) return(NULL)
  pieces <- lapply(csv_paths, function(p) {
    dt <- tryCatch(fread(p, header = TRUE, data.table = FALSE, showProgress = FALSE), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    df <- as.data.frame(dt)
    names(df) <- tolower(names(df))
    pathway_col <- intersect(c("pathway", "term", "gs", "name", "id", "description"), names(df))
    if (length(pathway_col) == 0) return(NULL)
    names(df)[match(pathway_col[1], names(df))] <- "pathway"
    if (!"nes" %in% names(df)) {
      alt <- intersect(c("normalized_enrichment_score", "enrichment", "score"), names(df))
      if (length(alt) > 0) names(df)[match(alt[1], names(df))] <- "nes"
    }
    if (!"fdr" %in% names(df)) {
      alt <- intersect(c("padj", "qvalue", "qval", "fdr_qval", "fdr_bh", "p.adjust"), names(df))
      if (length(alt) > 0) names(df)[match(alt[1], names(df))] <- "fdr"
    }
    df$msigdb_url <- vapply(df$pathway, msigdb_url, FUN.VALUE = character(1))
    df$comparison <- basename(p)
    df$comparison <- sub("_gsea_results\\.csv$", "", df$comparison, ignore.case = TRUE)
    df$source_file <- p
    df
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) return(NULL)
  do.call(rbind, pieces)
}

make_gsea_plot <- function(df, fdr_cutoff = 0.25, top_n = 10) {
  if (is.null(df)) return(NULL)
  if (!all(c("pathway", "nes", "comparison") %in% names(df))) return(NULL)
  plot_df <- df
  if ("fdr" %in% names(df)) {
    plot_df <- subset(plot_df, is.na(fdr) | fdr <= fdr_cutoff)
  }
  plot_df$abs_nes <- abs(plot_df$nes)
  plot_df <- plot_df[order(plot_df$comparison, -plot_df$abs_nes), ]
  plot_df <- do.call(rbind, lapply(split(plot_df, plot_df$comparison), function(d) head(d, top_n)))
  if (nrow(plot_df) == 0) return(NULL)
  plot_df$pathway <- factor(plot_df$pathway)
  ggplot(plot_df, aes(x = reorder(pathway, nes), y = nes, fill = nes > 0)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ comparison, scales = "free_y") +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#E4572E"), guide = "none") +
    labs(x = "Pathway", y = "NES") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

value_box <- function(title, value, subtitle = NULL) {
  tags$div(class = "value-box",
    tags$div(class = "vb-title", title),
    tags$div(class = "vb-value", value),
    if (!is.null(subtitle)) tags$div(class = "vb-subtitle", subtitle)
  )
}

# ---------------- Overview Helpers ---------------- #
scan_result_overview <- function(result_basename, progress = NULL) {
  if (is.null(result_basename) || !nzchar(result_basename)) return(list())
  
  base_dir <- file.path(OUTPUT_ROOT, result_basename)
  if (!dir.exists(base_dir)) return(list())
  
  celltypes <- list_celltype_dirs(result_basename)
  if (length(celltypes) == 0) return(list())
  
  n <- length(celltypes)
  res_list <- vector("list", n)
  names(res_list) <- celltypes
  
  for (i in seq_along(celltypes)) {
    ct <- celltypes[i]
    if (!is.null(progress)) progress(value = (i - 1) / n, detail = paste0("Reading ", ct, " (", i, "/", n, ")"))
    ct_dir <- file.path(base_dir, ct)
    files <- find_celltype_files(ct_dir)
    
    # Read DEG data - full for accurate counts
    deg_data <- NULL
    if (length(files$deg_tsv) > 0) {
      deg_data <- read_deg_table(files$deg_tsv)
    }
    
    # Read GSEA data - full for accurate counts
    gsea_data <- NULL
    if (length(files$gsea_csvs) > 0) {
      gsea_data <- read_gsea_results(files$gsea_csvs)
    }
    
    # Calculate summary stats
    contrasts <- if (!is.null(deg_data) && "contrast" %in% names(deg_data)) unique(deg_data$contrast) else character()
    comparisons <- if (!is.null(gsea_data) && "comparison" %in% names(gsea_data)) unique(gsea_data$comparison) else character()
    
    # Count significant genes per contrast
    sig_genes_by_contrast <- if (!is.null(deg_data) && all(c("contrast", "padj", "log2foldchange") %in% names(deg_data))) {
      split_data <- split(deg_data, deg_data$contrast)
      sapply(split_data, function(d) {
        sum(d$padj < PADJ_THRESH & abs(d$log2foldchange) >= LFC_THRESH, na.rm = TRUE)
      })
    } else {
      setNames(integer(length(contrasts)), contrasts)
    }
    
    # Count significant pathways per comparison
    sig_paths_by_comparison <- if (!is.null(gsea_data) && "fdr" %in% names(gsea_data)) {
      split_data <- split(gsea_data, gsea_data$comparison)
      sapply(split_data, function(d) {
        sum(d$fdr <= 0.25, na.rm = TRUE)
      })
    } else {
      setNames(integer(length(comparisons)), comparisons)
    }
    
    res_list[[i]] <- list(
      celltype = ct,
      major_class = major_celltype_of(ct),
      color = color_for_major_celltype(ct),
      contrasts = contrasts,
      comparisons = comparisons,
      sig_genes_by_contrast = sig_genes_by_contrast,
      sig_paths_by_comparison = sig_paths_by_comparison,
      total_sig_genes = sum(sig_genes_by_contrast),
      total_sig_paths = sum(sig_paths_by_comparison),
      has_deg_data = !is.null(deg_data),
      has_gsea_data = !is.null(gsea_data)
    )
  }
  if (!is.null(progress)) progress(value = 1, detail = "Done")
  res_list
}

make_overview_heatmap <- function(overview_data) {
  if (length(overview_data) == 0) return(NULL)
  
  # Create data frame for heatmap
  heatmap_data <- do.call(rbind, lapply(names(overview_data), function(ct) {
    data <- overview_data[[ct]]
    if (length(data$contrasts) == 0) return(NULL)
    
    data.frame(
      celltype = ct,
      contrast = na.omit(data$contrasts),
      sig_genes = data$sig_genes_by_contrast,
      stringsAsFactors = FALSE
    )
  }))
  
  if (nrow(heatmap_data) == 0) return(NULL)
  
  # Create heatmap
  ggplot(heatmap_data, aes(x = contrast, y = celltype, fill = sig_genes)) +
    geom_tile() +
    geom_text(aes(label = sig_genes), color = "#263238", size = 2.4, fontface = "plain") +
    scale_fill_gradient2(
      low = "#f7f7f7", 
      mid = "#92c5de", 
      high = "#0571b0",
      midpoint = max(heatmap_data$sig_genes) / 2,
      name = "Significant\nGenes"
    ) +
    labs(x = "Contrast", y = "Cell Type", title = "Significant Genes by Cell Type and Contrast") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

make_overview_summary_plot <- function(overview_data) {
  if (length(overview_data) == 0) return(NULL)
  
  # Create summary data
  summary_data <- data.frame(
    celltype = names(overview_data),
    total_sig_genes = sapply(overview_data, function(x) x$total_sig_genes),
    total_sig_paths = sapply(overview_data, function(x) x$total_sig_paths),
    n_contrasts = sapply(overview_data, function(x) length(x$contrasts)),
    n_comparisons = sapply(overview_data, function(x) length(x$comparisons)),
    stringsAsFactors = FALSE
  )
  
  # Create bar plot
  ggplot(summary_data, aes(y = reorder(celltype, total_sig_genes), x = total_sig_genes)) +
    geom_col(fill = "#2E86AB", alpha = 0.8) +
    geom_text(aes(label = total_sig_genes), hjust = -0.2, size = 3, fontface = "bold") +
    labs(y = "Cell Type", x = "Total Significant Genes", 
         title = "Total Significant Genes by Cell Type") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )
} 