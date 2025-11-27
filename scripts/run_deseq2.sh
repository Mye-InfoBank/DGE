#/bin/bash

/nfs/data/COST_IBD/downstream_tasks/dea/deseq2.R \
 --input /nfs/data/COST_IBD/versions/HNC/v_Hamburg/atlas_scvi_extended_with_subatlases_COUNTS_final_annotations_Lymphocytes.h5ad \
 --category_column site \
 --cell_type_column final_annotation_combined \
 --sample_column batch \
 --count_assay X \
 --one_vs_all \
 --output_directory /nfs/data/COST_IBD/downstream_tasks/dea/output/hnc_site_oneVSall \
 --remove_groups 'unknown' 'Larynx/Hypopharynx'
