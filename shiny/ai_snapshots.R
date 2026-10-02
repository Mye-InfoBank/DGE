# Helpers for the "AI snapshot": a JSON export of the current view that the user
# downloads and uploads to an LLM of their choice.

# Data frame -> list of row records (JSON-safe: factors as strings, non-finite as NA)
records <- function(df) {
  if (is.null(df) || !nrow(df)) return(list())

  df[] <- lapply(df, function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.numeric(x)) x[!is.finite(x)] <- NA
    x
  })

  unname(split(df, seq_len(nrow(df))))
}

# Static context so an LLM can interpret the snapshot without further input
ai_snapshot_llm_context <- function() {
  list(
    purpose = paste(
      "This file is a snapshot of the COST IBD MyeInfoBank differential expression explorer",
      "(Shiny app). It contains what the user was looking at when it was exported."
    ),
    instructions_for_assistant = list(
      "Base your answers on this file. If something is not in it, say so instead of guessing.",
      "Start by briefly stating which cell type, contrast and comparison the snapshot shows.",
      "Treat top-N lists as truncated: absence from a list does not mean a gene or pathway is not significant.",
      "Distinguish statistical significance (padj, FDR) from effect size (log2FC, NES).",
      "Offer biological interpretation as hypotheses, not conclusions; mention caveats such as sample size.",
      "When citing genes or pathways, quote the numeric values (log2FC, padj, NES, FDR) from this file.",
      "If the user's question is vague, summarise the main findings and ask what they want to explore."
    ),
    study_context = list(
      project = "COST IBD MyeInfoBank (inflammatory bowel disease single-cell atlas)",
      unit_of_analysis = paste(
        "Pseudobulk: raw single-cell counts are summed per sample within each cell type",
        "(and per covariate level), so statistical power depends on the number of samples, not cells."
      ),
      cell_type_naming = "Cell type labels may be prefixed by tissue location or split group, e.g. 'Ileum/Plasma cells'.",
      contrast_naming = paste(
        "'A_vs_B' is a pairwise comparison of group A against reference group B.",
        "'A_vs_rest' compares group A against all other groups pooled.",
        "Group abbreviations are dataset-specific; do not assume their meaning, ask the user if unclear."
      )
    ),
    methods = list(
      differential_expression = paste(
        "DESeq2 (Wald test) per cell type on pseudobulk counts.",
        "Pseudobulk samples with fewer than 3 cells are dropped, cell types need at least 4 remaining pseudobulk samples,",
        "and genes with zero counts are removed.",
        "Design: ~ covariates + log2(cells per pseudobulk) + condition.",
        "log2FC values are NOT shrunk (no lfcShrink), so low-count genes can show inflated fold changes.",
        "padj is Benjamini-Hochberg with DESeq2's independent filtering (padj can be NA for filtered genes).",
        "'stat' is the Wald statistic. The significance cut-off in the app is padj < 0.05 and |log2FC| >= 1,",
        "applied after testing (the test itself is against log2FC = 0)."
      ),
      gsea = paste(
        "clusterProfiler GSEA on all tested genes ranked by the DESeq2 Wald statistic (duplicates reduced to max |stat|),",
        "against MSigDB Hallmark (H, human) gene sets, default gene-set size limits (10-500), eps = 0.",
        "The export used pvalueCutoff = 0.05 on the adjusted p-value, so ONLY significant pathways are saved:",
        "a pathway missing from the lists was either not significant or not testable, not 'unreported'.",
        "A contrast with no significant pathways has no GSEA entry at all."
      )
    ),
    field_glossary = list(
      baseMean = "Mean of normalised counts across samples (expression level).",
      log2foldchange = "log2 fold change of first vs second group in the contrast (positive = higher in first group).",
      pvalue = "Unadjusted Wald test p-value.",
      padj = "Benjamini-Hochberg adjusted p-value.",
      nes = "Normalised enrichment score (positive = enriched in first group).",
      fdr = "FDR-adjusted p-value of the GSEA result.",
      size = "Number of genes of the pathway present in the ranked list.",
      significant_genes = "Genes with padj < thresholds$padj and |log2FC| >= thresholds$abs_log2fc.",
      significant_pathways = "Pathways with fdr <= thresholds$gsea_fdr (all exported pathways already have fdr < 0.05)."
    ),
    truncation = list(
      top_up = "up to 50 genes, ordered by padj",
      top_down = "up to 50 genes, ordered by padj",
      top_by_padj = "up to 100 genes, lowest padj regardless of direction",
      top_positive_nes = "up to 25 pathways, highest NES",
      top_negative_nes = "up to 25 pathways, lowest NES"
    ),
    suggested_questions = list(
      "Summarise the main biology in this cell type and contrast.",
      "Which up- and down-regulated genes are most notable, and what do they suggest?",
      "How do the enriched pathways relate to the differentially expressed genes?",
      "What are the limitations of this result (sample size, truncation, thresholds)?"
    )
  )
}
