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
#   GSEA on DESeq2 results                                                                                                                                 #
#   Author: Alexander Dietrich, Leon Hafner                                                                                                                #
#   Prepared for the EU COST Action MyeInfoBank                                                                                                            #
############################################################################################################################################################
")

# ---------------- Load packages ---------------- #
message("Loading necessary packages ...")
suppressPackageStartupMessages({
  library(argparse)
  library(clusterProfiler)
  library(msigdbr)
  library(ggplot2)
  library(tools)
})

# ---------------- Helper functions ---------------- #
sanitize <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

# ---------------- Argument parsing ---------------- #
parser <- ArgumentParser(description = "GSEA on DESeq2 results")
parser$add_argument(
  "--input_dir",
  required = TRUE,
  help = "Directory with DESeq2 results (tsv files)"
)

args <- parser$parse_args()
input_dir <- args$input_dir

# ---------------- Log file setup ---------------- #
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file  <- file.path(input_dir, paste0("gsea_", timestamp, ".log"))

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

Sys.time()

message("\n================ Selected Parameters ================\n")
message("Input directory:      ", input_dir)
message("====================================================\n")

# ---------------- Main script ---------------- #

message("Searching for DESeq2 result folders ...")
deseq_result_files <- list.files(
  input_dir,
  pattern = "DESeq2_results.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(deseq_result_files) == 0) {
  stop("No DESeq2 result files found in: ", input_dir)
}

message("Found ", length(deseq_result_files), " result file(s).")

message("Loading MSigDB Hallmark gene sets ...")
msigdbr_df <- msigdbr(species = "Homo sapiens", collection = "H")
term2gene  <- msigdbr_df[, c("gs_name", "gene_symbol")]
colnames(term2gene) <- c("term", "gene")

for (file in deseq_result_files) {
  message("\n######## Processing file: ", file)

  parent_dir <- basename(dirname(file))
  out_dir    <- dirname(file)

  res <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)

  required_cols <- c("stat", "contrast", "gene")
  missing_cols <- setdiff(required_cols, colnames(res))

  if (length(missing_cols) > 0) {
    message("Skipping ", file, ": missing required columns: ", paste(missing_cols, collapse = ", "))
    next
  }

  res <- res[!is.na(res$stat), ]

  if (nrow(res) == 0) {
    message("Skipping ", file, ": no rows with non-NA stat.")
    next
  }

  comparisons <- unique(res$contrast)
  comparisons <- comparisons[!is.na(comparisons)]

  message("Found ", length(comparisons), " comparison(s):")
  message(paste(comparisons, collapse = ", "))

  for (comp in comparisons) {
    message("\n---- Running GSEA for contrast: ", comp)

    comp_clean <- sanitize(comp)

    res_comp <- res[res$contrast == comp, ]

    gl <- res_comp$stat
    names(gl) <- res_comp$gene

    gl <- gl[!is.na(gl)]
    gl <- gl[!is.na(names(gl))]
    gl <- gl[names(gl) != ""]

    if (length(gl) < 10) {
      message("Skipping ", parent_dir, " / ", comp, ": fewer than 10 ranked genes.")
      next
    }

    # Deduplicate genes by keeping the gene with maximum absolute statistic
    if (anyDuplicated(names(gl))) {
      o <- order(abs(gl), decreasing = TRUE)
      gl <- gl[o][!duplicated(names(gl[o]))]
    }

    gl <- sort(gl, decreasing = TRUE)

    gseaRes <- tryCatch(
      {
        GSEA(
          gl,
          TERM2GENE = term2gene,
          pvalueCutoff = 0.05,
          verbose = FALSE,
          eps = 0
        )
      },
      error = function(e) {
        message("GSEA error for ", parent_dir, " / ", comp, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(gseaRes)) {
      message("No GSEA result object for ", parent_dir, " / ", comp)
      next
    }

    gsea_df <- as.data.frame(gseaRes)

    if (nrow(gsea_df) == 0) {
      message("No significant gene sets for ", parent_dir, " / ", comp)
      next
    }

    # Save GSEA results table
    out_csv <- file.path(out_dir, paste0(comp_clean, "_gsea_results.csv"))
    write.csv(gsea_df, out_csv, row.names = FALSE)
    message("Saved results table: ", out_csv)

    # Save simple pathway summary plot using ggplot2 only
    plot_df <- gsea_df
    plot_df <- plot_df[!is.na(plot_df$p.adjust) & !is.na(plot_df$NES), ]
    plot_df <- plot_df[order(plot_df$p.adjust), ]
    plot_df <- head(plot_df, 20)

    if (nrow(plot_df) > 0) {
      plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))

      out_summary_png <- file.path(out_dir, paste0(comp_clean, "_gsea_top_pathways.png"))

      p_summary <- ggplot(plot_df, aes(x = NES, y = Description)) +
        geom_point(aes(size = setSize, color = p.adjust)) +
        theme_bw() +
        labs(
          title = paste0("Top enriched Hallmark pathways: ", comp),
          x = "Normalized enrichment score (NES)",
          y = NULL,
          size = "Gene set size",
          color = "Adjusted p-value"
        )

      tryCatch(
        {
          ggsave(out_summary_png, plot = p_summary, width = 9, height = 7)
          message("Saved summary plot: ", out_summary_png)
        },
        error = function(e) {
          message("Could not save summary plot for ", parent_dir, " / ", comp, ": ", conditionMessage(e))
        }
      )
    } else {
      message("No valid rows for summary plot for ", parent_dir, " / ", comp)
    }

    message("Skipping enrichplot::gseaplot2 running-score plots for ", comp)
  }
}

message("\nAll done. Results written to: ", input_dir)
Sys.time()

# ---------------- Restore console ---------------- #
sink(type = "message")
sink(type = "output")
close(log_con)