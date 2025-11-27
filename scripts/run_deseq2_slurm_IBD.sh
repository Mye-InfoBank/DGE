#!/bin/bash
#SBATCH --job-name=deseq2_ibd_condition
#SBATCH --output=/nfs/data/COST_IBD/downstream_tasks/dea/logs/deseq2_ibd_condition_%j.out
#SBATCH --error=/nfs/data/COST_IBD/downstream_tasks/dea/logs/deseq2_ibd_condition_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=150G
#SBATCH --partition=exbio-cpu

# Create logs directory if it doesn't exist
mkdir -p /nfs/data/COST_IBD/downstream_tasks/dea/logs


echo "Starting DESeq2 analysis for IBD condition"
echo "Input file: /nfs/data/COST_IBD/versions/IBD/07_00_00/build/results/finalized/merged_meta_Hamburg_annotated_noNA.h5ad"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"

# Run the DESeq2 analysis with IBD parameters
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/deseq2.R \
    --input "/nfs/data/COST_IBD/versions/IBD/07_00_00/build/results/finalized/merged_meta_Hamburg_annotated_noNA.h5ad" \
    --category_column condition \
    --cell_type_column annotation.large \
    --sample_column sample \
    --count_assay X \
    --output_directory /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_condition \
    --remove_groups 'NA' \
    --remove_cell_type 'NA' 

# Check exit status
if [ $? -eq 0 ]; then
    echo "DESeq2 analysis completed successfully at $(date)"
else
    echo "DESeq2 analysis failed at $(date)"
    exit 1
fi
