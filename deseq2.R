#!/usr/bin/env Rscript

cat("
############################################################################################################################################################
#           ## ###    ###                                 %%%%%%%  %%%  %%%%% %%%%%%%%%  %%%%%%   ########   #####    ####  #### ####  ###  #              #   
#        ###    ###  #### ####  ###   #####                  %      %%%   %    %%    %% %%     %%  ##    ##    ####     ###  ##   ## ###      ###          #   
#    ###        # #### ##  ##   ##   #     #   ########      %      %  %  %    %%%%%   %%      %%  #######    ##  ##    # ## ##   #####          ###       #   
#      ####     #  ##  ##   ##  #   #########                %      %   %%%    %%  %    %%     %%  ##    ##  #######    #  ####   ##  ##       ###         #   
#          ###  #      ##    ###     ##    ##                %      %    %%    %%        %%   %%   ##   ###  #     ##   #    ##   ##   ##   ###            #   
#               ##    ###     #        ####                %%%%%    %%    %    %%%%        %%%     #####    ###    ##  ###    #   ###    #                 #
#                         ######                                                                                                                           #
#  %                             #                             %%                         #  ##                          %%                           ##   #
# %%%*%%%               #######*###*###               %%%%%%%%#%%##%%%              ############*###               %%%% %%%#%%%#%%%              ########* #
#  ++ +++#%%#       ## ####*+++  ++ +++*##*        %% #%%*++++ ++  +++%%#        ## *##*++++  +  +++###           #%%%*  ++  ++++*+*%%#      ####*##*+ +   #
#         ++%%#    ####*+++  **+++**     +###    %%%%# ++  *+***       + %%#   ##### ++   #*++***    ++##*    %%%%# ++  *+**+++     ++*%%%   ####*         #
#            + ####++++       *+    **     ++ %%% ++++       *+**         + ####++++     *+    *+       + %%%%+**++        *           ++####++++          #
#             ######          *+    #*       %%%%%*         #+  **         #####         *               %%%%%#            *            #####+             #
#            # ***++          *+    *+     % +***+         #*******      ##+***+         *             %% *##*+            *         ####++++%%            #
#          ####*   %%%       #**###*+    %%%%%    ###     #**#  ##*#   #####    %%%       **###*+     %%%%    ###       ###**##      ####*   ***%%         #
#     #####+**++    ++%%%            %%%% ***+    ++####            ####+**++   ++*%%%      +++  %%%%%***+    +++###            #####++++      ###+%%% %%% #
# ### **#*+           +++%%%#%%%%%%%%*%%#+           +++##*########**##+           +++%%##%%% %%%#*##+           +++###*#######*+++++              +*+*%%* #
# +*++                    ++ ++++++++                    + ++++++++                   +++ +*+*++++                   ++ +*+++++                         +  #
############################################################################################################################################################ 
#   DESeq2 Differential Gene Expression Analysis                                                                                                           #
#   Author: Alexander Dietrich, Leon Hafner                                                                                                                #
#   Prepared for the EU COST Action MyeInfoBank                                                                                                            #
############################################################################################################################################################
")

library(argparse)

# ---------------- Helper functions ---------------- #

check_deseq_design <- function(se, formula_str) {
  meta <- as.data.frame(SummarizedExperiment::colData(se))
  vars <- all.vars(as.formula(formula_str))
  
  # 1. check constant variables
  no_var <- vars[vapply(meta[vars], function(x) length(unique(x)) == 1, logical(1))]
  if (length(no_var) > 0) {
    return(list(ok = FALSE,
                reason = paste("Variables with no variance:",
                               paste(no_var, collapse = ", "))))
  }
  
  # 2. try to build DESeqDataSet
  dds <- tryCatch(
    DESeqDataSet(se, design = as.formula(formula_str)),
    error = function(e) return(NULL)
  )
  if (is.null(dds)) {
    return(list(ok = FALSE, reason = "DESeqDataSet creation failed"))
  }
  
  # 3. check rank
  mm <- model.matrix(design(dds), data = colData(dds))
  qr_mm <- qr(mm)
  full_rank <- qr_mm$rank == ncol(mm)
  
  if (!full_rank) {
    return(list(ok = FALSE,
                reason = paste0("Model matrix not full rank (rank ",
                                qr_mm$rank, " vs ", ncol(mm), " columns)")))
  }
  
  list(ok = TRUE, reason = "✅ Design is valid", dds = dds)
}

save_volcano <- function(res_df, out_html, lfc_thresh = 1, padj_thresh = 0.05) {
  message("Saving DESeq2 result plot... ")
 
  res_df$significance <- "NS"
  res_df$significance[res_df$padj < padj_thresh & res_df$log2FoldChange > lfc_thresh] <- "Up"
  res_df$significance[res_df$padj < padj_thresh & res_df$log2FoldChange < -lfc_thresh] <- "Down"
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance, 
                          text = paste("Gene:", gene,
                                       "<br>log2FC:", round(log2FoldChange, 3),
                                       "<br>padj:", signif(padj, 3)))) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = 'dashed')+
    geom_vline(xintercept = -lfc_thresh, linetype = 'dashed')+
    geom_vline(xintercept = lfc_thresh, linetype = 'dashed')+
    geom_point(alpha = 0.8) +
    facet_wrap(~contrast)+
    theme_bw() +
    labs(x = "log2 Fold Change", y = "-log10 adjusted p-value")+
    # scale_color_manual(values = c(
    #   "Up"   = "#D55E00",   # orange-red
    #   "Down" = "#0072B2",   # blue
    #   "NS"   = "#999999"    # grey
    # ))+
    # scale_color_manual(values = c(
    #   "Up"   = "#E76F51",   # terracotta
    #   "Down" = "#2A9D8F",   # sea-green
    #   "NS"   = "#BDBDBD"    # light grey
    # ))+
    scale_color_manual(values = c(
      "Up"   = "#FF4E50",   # coral red
      "Down" = "#1E90FF",   # dodger blue
      "NS"   = "#C0C0C0"    # silver grey
    ))

  # create custom plot title with cell-type and contrast
  plot_title <- paste("Volcano Plot for", unique(obs[[cell_type_column]]), "by", unique(obs[[category_column]]))
  p <- p + ggtitle(plot_title)
  
  ggsave(p, filename = gsub(".html$", ".png", out_html), width = 14, height = 8)
  #htmlwidget <- plotly::ggplotly(p, tooltip = "text")
  #htmlwidgets::saveWidget(htmlwidget, out_html, selfcontained = TRUE)
}

# ---------------- Main script ---------------- #

snowparam <- BiocParallel::SnowParam(workers = 6, type = "SOCK")
BiocParallel::register(snowparam, default = TRUE)

parser <- ArgumentParser(description = "Pseudobulk DESeq2 analysis per cell type")
parser$add_argument("--input", required = TRUE, help = "Input .h5ad dataset")
parser$add_argument("--category_column", required = TRUE,
                    help = "Column to compare groups (e.g. condition)")
parser$add_argument("--cell_type_column", required = TRUE,
                    help = "Column with cell type labels")
parser$add_argument("--covariates", nargs = "*", default = NULL,
                    help = "Optional covariates")
parser$add_argument("--sample_column", nargs = "*", default = 'sample',
                    help = "Name of the column with sample identifiers for pseudobulking (default: sample)")      
parser$add_argument("--count_assay", nargs = "*", default = 'counts',
                    help = "Name of the assay with raw counts (default: counts)")   
parser$add_argument("--remove_groups", nargs = "*", default = NULL,
                    help = "Optional groups to remove from the category_column before analysis")                                     
parser$add_argument("--one_vs_all", action = "store_true", default = FALSE,
                    help = "If set, perform one-vs-all comparisons instead of all pairwise comparisons")
parser$add_argument("--output_directory", required = TRUE,
                    help = "Directory to save results")

args <- parser$parse_args()

category_column <- args$category_column
cell_type_column <- args$cell_type_column
covariates <- args$covariates
sample_column <- args$sample_column
count_assay <- args$count_assay
outdir <- args$output_directory
remove_groups <- args$remove_groups
one_vs_all <- args$one_vs_all

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message('Loading necessary packages ...')
suppressPackageStartupMessages({
  library(zellkonverter)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(scuttle)
  library(ggplot2)
  library(plotly)
})

# ---------------- Log file setup ---------------- #
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(outdir, paste0("deseq2_analysis_", timestamp, ".log"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Redirect output and messages to log
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")   # cat() and print() go here
sink(log_con, type = "message")  # message(), warning(), and errors go here
Sys.time()

message("\n================ Selected Parameters ================\n")
message("Input file:           ", args$input)
message("Category column:      ", args$category_column)
message("Cell type column:     ", args$cell_type_column)
message("Sample column:        ", args$sample_column)
message("Count assay:          ", args$count_assay)
message("Covariates:           ", ifelse(length(args$covariates) > 0,
                                         paste(args$covariates, collapse = ", "),
                                         "None"))
message("Remove groups:        ", ifelse(length(args$remove_groups) > 0,
                                         paste(args$remove_groups, collapse = ", "),
                                         "None"))
message("One-vs-all mode:      ", ifelse(args$one_vs_all, "Yes", "No"))                                         
message("Output directory:     ", args$output_directory)
message("====================================================\n")

# ---------------- actual DESeq2 starts here ---------------- #

message("Reading input: ", args$input)
adata <- suppressWarnings(zellkonverter::readH5AD(args$input))

if (!is.null(remove_groups)) {
  message("Removing groups from ", category_column, ": ", paste(remove_groups, collapse = ", "))
  adata <- adata[, !(adata[[category_column]] %in% remove_groups)]

  # remove empty factor levels from category_column
  if (is.factor(adata[[category_column]])) {
    adata[[category_column]] <- droplevels(adata[[category_column]])
  }
}
adata <- adata[, adata$dataset != "Devlin"]

if(count_assay %in% names(assays(adata))){
  message(paste0('Using ', count_assay, ' assay for DGE testing ...'))
  if(count_assay != 'counts'){
      assays(adata)$counts <- assays(adata)[[count_assay]]
  }
}

pseudobulk_group_cols <- c(category_column, cell_type_column, sample_column)
if (!is.null(covariates)) {
  pseudobulk_group_cols <- c(pseudobulk_group_cols, covariates)
}
pseudobulk_group_cols <- unlist(pseudobulk_group_cols)
min_cells_in_group <- 3
print(pseudobulk_group_cols)

obs <- colData(adata) |> data.frame(check.names = F, row.names = colnames(adata)) 
print(unique(obs[[cell_type_column]]))

celltypes <- unique(adata[[cell_type_column]])
celltypes <- celltypes[!is.na(celltypes)]

all_results <- lapply(celltypes, function(group){
  message("\n######## Starting DGE analysis for ", group)
  
  group_obs <- obs[obs[[cell_type_column]] == group & !is.na(obs[[cell_type_column]]), ]
  if (nrow(group_obs) < 4) {
    message("Skipping ", group, ": only ", nrow(group_obs), " cell(s).")
    return(NULL)
  }
  
  # Subset data
  sce_sub <- adata[, !is.na(adata[[cell_type_column]]) & adata[[cell_type_column]] == group]
  
  # Aggregate to pseudobulk
  pb <- scuttle::aggregateAcrossCells(
    sce_sub,
    ids = group_obs[, pseudobulk_group_cols, drop = FALSE] |> 
      unite("pb_id", everything(), remove = FALSE) |> pull("pb_id"),
    use.assay.type = "counts",
    statistics = "sum",
    store.number = "ncells"
  )
  colData(pb)$log2_cells_pseudobulk <- log2(as.numeric(colData(pb)$ncells))
  
  # Filter by min cells
  pb <- pb[, colData(pb)$ncells >= min_cells_in_group]
  if (ncol(pb) < 4) {
    message("Skipping ", group, ": only ", ncol(pb), " pseudobulk sample(s) after filtering.")
    return(NULL)
  }
  
  # Remove all-zero genes
  pb <- pb[which(rowSums(assays(pb)$counts) > 0), ]
  
  res_list <- list()
  
  if (!one_vs_all) {
    # -------- Pairwise mode (default) -------- #
    # Build design
    design_vars <- c("log2_cells_pseudobulk", category_column)
    if (!is.null(covariates)) {
      design_vars <- c(covariates, "log2_cells_pseudobulk", category_column)
    }
    formula_str <- as.formula(paste0("~", paste(design_vars, collapse = " + ")))
    
    diag <- check_deseq_design(pb, formula_str)
    if (!diag$ok) {
      message("❌ Skipping ", group, ": ", diag$reason)
      return(NULL)
    }
    
    dds <- diag$dds
    message("Running DESeq for ", group, " with ", ncol(dds), " samples.")
    dds <- estimateSizeFactors(dds)
    dds <- DESeq2::DESeq(dds, parallel = TRUE)
    
    possible_conditions <- levels(dds[[category_column]])
    contrasts <- combn(possible_conditions, 2, simplify = FALSE)
    
    for (c in contrasts) {
      contrast_vector <- c(category_column, c[1], c[2])
      contrast_name <- paste(c[1], "_vs_", c[2], sep = "")
      
      res <- results(dds, contrast = contrast_vector)
      res_df <- as.data.frame(res)
      res_df$contrast <- contrast_name
      res_df$gene <- rownames(res_df)
      res_df$celltype <- group
      res_list[[contrast_name]] <- res_df
    }
  } else {
    # -------- One-vs-all mode -------- #
    possible_conditions <- levels(pb[[category_column]])
    
    for (cond in possible_conditions) {
      message("One-vs-all: ", cond, " vs rest")
      
      # create new 2-level factor
      colData(pb)$tmp <- ifelse(colData(pb)[[category_column]] == cond, cond, "other")
      pb$tmp <- factor(pb$tmp)
      
      # build new design
      design_vars_ova <- c("log2_cells_pseudobulk", "tmp")
      if (!is.null(covariates)) {
        design_vars_ova <- c(covariates, "log2_cells_pseudobulk", "tmp")
      }
      formula_ova <- as.formula(paste0("~", paste(design_vars_ova, collapse = " + ")))
      
      diag_ova <- check_deseq_design(pb, formula_ova)
      if (!diag_ova$ok) {
        message("❌ Skipping one-vs-all for ", cond, ": ", diag_ova$reason)
        next
      }
      
      dds_ova <- diag_ova$dds
      dds_ova <- estimateSizeFactors(dds_ova)
      dds_ova <- DESeq2::DESeq(dds_ova, parallel = TRUE)
      
      res <- results(dds_ova, contrast = c("tmp", cond, "other"))
      res_df <- as.data.frame(res)
      res_df$contrast <- paste0(cond, "_vs_rest")
      res_df$gene <- rownames(res_df)
      res_df$celltype <- group
      res_list[[paste0(cond, "_vs_rest")]] <- res_df
    }
  }
  
  res_celltype <- bind_rows(res_list)
  
  # Save per-celltype results
  group_dir <- file.path(outdir, group)
  dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_tsv <- file.path(group_dir, paste0("DESeq2_results.tsv"))
  write.table(res_celltype, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  out_html <- file.path(group_dir, paste0("volcano.html"))
  save_volcano(res_celltype, out_html)
  
  return(res_celltype)
})


# Save combined results
combined <- bind_rows(all_results)
out_combined <- file.path(outdir, "DESeq2_all_results.tsv")
write.table(combined, out_combined, sep = "\t", quote = FALSE, row.names = FALSE)

message("\nAll done. Results written to: ", outdir)
Sys.time()

# ---------------- Restore console ---------------- #
sink(type = "message")
sink(type = "output")
close(log_con)
