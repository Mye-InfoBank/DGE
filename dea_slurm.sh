#!/bin/bash
#SBATCH --job-name=dea_hnc
#SBATCH --output=dea_hnc_%A_%a.out
#SBATCH --error=dea_hnc_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-5
#SBATCH --partition=standard

INPUT="/nfs/data/COST_IBD/versions/HNC/v_Hamburg/atlas_scvi_extended_with_subatlases_COUNTS_final_annotations.h5ad"
SCRIPT="/nfs/data/COST_IBD/downstream_tasks/dea/deseq2.R"
CELLTYPE="final_annotation_combined"
SAMPLE="batch"
ASSAY="raw_counts"
OUTDIR="/nfs/data/COST_IBD/downstream_tasks/dea/output"

# Define categories and output subfolders
CATEGORIES=(HPV EBV classificatio sex site disease)
OUTPUTS=(hnc_hpv hnc_ebv hnc_classification hnc_sex hnc_site hnc_disease)

CATEGORY=${CATEGORIES[$SLURM_ARRAY_TASK_ID]}
OUTPUT=${OUTPUTS[$SLURM_ARRAY_TASK_ID]}

Rscript $SCRIPT --input $INPUT --category_column $CATEGORY --cell_type_column $CELLTYPE --sample_column $SAMPLE --count_assay $ASSAY --output_directory $OUTDIR/$OUTPUT
