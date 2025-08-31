#/bin/bash

/nfs/data/COST_IBD/downstream_tasks/dea/deseq2.R --input /nfs/data/COST_IBD/versions/IBD/06_00_00/build/results/finalized/merged_annotated_Hamburg_noNA.h5ad --category_column Disease.location --cell_type_column annotation.large --sample_column batch --count_assay X --output_directory /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_location --remove_groups Unknown, Inapplicable 

#/nfs/data/COST_IBD/downstream_tasks/dea/deseq2.R --input /nfs/data/COST_IBD/versions/HNC/v_Hamburg/atlas_scvi_extended_with_subatlases_COUNTS_final_annotations_700k.h5ad --category_column site --cell_type_column final_annotation_combined --sample_column batch --count_assay X --output_directory /nfs/data/COST_IBD/downstream_tasks/dea/output/hnc_site
