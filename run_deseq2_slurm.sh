#!/bin/bash
#SBATCH --job-name=deseq2_hnc_site
#SBATCH --output=/nfs/data/COST_IBD/downstream_tasks/dea/logs/deseq2_hnc_site_%j.out
#SBATCH --error=/nfs/data/COST_IBD/downstream_tasks/dea/logs/deseq2_hnc_site_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --partition=exbio-cpu
#SBATCH --array=1-4

# Create logs directory if it doesn't exist
mkdir -p /nfs/data/COST_IBD/downstream_tasks/dea/logs

# Define input files for each array job
declare -a INPUT_FILES=(
    "/nfs/data/COST_IBD/versions/HNC/v_Hamburg/atlas_scvi_extended_with_subatlases_COUNTS_final_annotations_Fibroblast.h5ad"
    "/nfs/data/COST_IBD/versions/HNC/v_Hamburg/atlas_scvi_extended_with_subatlases_COUNTS_final_annotations_Myeloid.h5ad"
    "/nfs/data/COST_IBD/versions/HNC/v_Hamburg/atlas_scvi_extended_with_subatlases_COUNTS_final_annotations_Lymphocytes.h5ad"
    "/nfs/data/COST_IBD/versions/HNC/v_Hamburg/atlas_scvi_extended_with_subatlases_COUNTS_final_annotations_Endothelial.h5ad"
)

# Get the input file for this array task
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID-1]}

# Extract cell type name from filename for logging
CELL_TYPE=$(basename "$INPUT_FILE" | sed 's/.*_\([^_]*\)\.h5ad$/\1/')

echo "Starting DESeq2 analysis for cell type: $CELL_TYPE"
echo "Input file: $INPUT_FILE"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Date: $(date)"

# Run the DESeq2 analysis
/nfs/data/COST_IBD/downstream_tasks/dea/deseq2.R \
    --input "$INPUT_FILE" \
    --category_column site \
    --cell_type_column final_annotation_combined \
    --sample_column batch \
    --count_assay X \
    --one_vs_all \
    --output_directory /nfs/data/COST_IBD/downstream_tasks/dea/output/hnc_site_oneVSall \
    --remove_groups 'unknown' 'Larynx/Hypopharynx'

# Check exit status
if [ $? -eq 0 ]; then
    echo "DESeq2 analysis completed successfully for $CELL_TYPE at $(date)"
else
    echo "DESeq2 analysis failed for $CELL_TYPE at $(date)"
    exit 1
fi
