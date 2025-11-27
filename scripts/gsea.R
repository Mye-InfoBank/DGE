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
  library(enrichplot)
  library(msigdbr)
  library(ggplot2)
  library(tools)
})

# ---------------- Helper functions ---------------- #
sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# ---------------- Argument parsing ---------------- #
parser <- ArgumentParser(description = "GSEA on DESeq2 results")
parser$add_argument("--input_dir", required = TRUE,
                    help = "Directory with DESeq2 results (tsv files)")

args <- parser$parse_args()

input_dir <- args$input_dir

# ---------------- Log file setup ---------------- #
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file  <- file.path(input_dir, paste0("gsea_", timestamp, ".log"))

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")   # cat(), print()
sink(log_con, type = "message")  # message(), warning(), error()
Sys.time()

message("\n================ Selected Parameters ================\n")
message("Input directory:      ", input_dir)
message("====================================================\n")

# ---------------- Main script ---------------- #

message("Searching for DESeq2 result folders ...")
deseq_result_files <- list.files(input_dir,
                                 pattern = "DESeq2_results.tsv$",
                                 full.names = TRUE,
                                 recursive = TRUE)
if (length(deseq_result_files) == 0) {
  stop("❌ No DESeq2 result folders found in: ", input_dir)
}
message("Found ", length(deseq_result_files), " result folders.")

message("Loading MSigDB Hallmark gene sets ...")
msigdbr_df <- msigdbr(species = "Homo sapiens", collection = "H")
term2gene  <- msigdbr_df[, c("gs_name", "gene_symbol")]
colnames(term2gene) <- c("term", "gene")

for (file in deseq_result_files) {
  message("\n######## Processing file: ", file)
  
  parent_dir <- basename(dirname(file))
  out_dir    <- dirname(file)
  
  res <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  res <- res[!is.na(res$stat), ]
  
  comparisons <- unique(res$contrast)
  message("Found ", length(comparisons), " comparisons: \n",
          paste(comparisons, collapse = ", "))
  
  for (comp in comparisons) {
    message("\n---- Running GSEA for contrast: ", comp)
    
    comp_clean <- gsub('/','|', comp)
    
    res_comp <- res[res$contrast == comp, ]
    gl <- res_comp$stat
    names(gl) <- res_comp$gene
    gl <- gl[!is.na(gl)]
    
    # Deduplicate by max |stat|
    if (anyDuplicated(names(gl))) {
      o <- order(abs(gl), decreasing = TRUE)
      gl <- gl[o][!duplicated(names(gl[o]))]
    }
    gl <- sort(gl, decreasing = TRUE)
    
    gseaRes <- tryCatch(
      GSEA(gl, TERM2GENE = term2gene,
           pvalueCutoff = 0.05, verbose = FALSE, eps = 0),
      error = function(e) {
        message("❌ GSEA error: ", conditionMessage(e))
        NULL
      }
    )
    gsea_df <- as.data.frame(gseaRes)
    if (nrow(gsea_df) == 0) {
      message("No significant gene sets for ", parent_dir, " / ", comp)
      next
    }
    
    # Save table
    out_csv <- file.path(out_dir, paste0(comp_clean, "_gsea_results.csv"))
    write.csv(gsea_df, out_csv, row.names = FALSE)
    message("Saved results table: ", out_csv)
    
    # Save plots
    gsea_plot_dir <- file.path(out_dir, "gsea_plots", comp_clean)
    dir.create(gsea_plot_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (i in seq_len(nrow(gsea_df))) {
      term <- gsea_df$ID[i]
      es <- gsea_df$enrichmentScore[i]
      p <- gseaplot2(gseaRes, geneSetID = term, title = paste0(term, " (ES: ", round(es, 4),")"), pvalue_table = TRUE)
      out_png <- file.path(gsea_plot_dir, paste0(sanitize(term), ".png"))
      ggsave(out_png, plot = p, width = 8, height = 6)
    }
    message("Saved ", nrow(gsea_df), " plots for ", comp)
  }
}

message("\n✅ All done. Results written to: ", input_dir)
Sys.time()

# ---------------- Restore console ---------------- #
sink(type = "message")
sink(type = "output")
close(log_con)