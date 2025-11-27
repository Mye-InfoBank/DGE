#/bin/bash

/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_condition

/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Cecum
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Colon
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Left_Colon
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Rectum
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Right_Colon
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Sigmoid_Colon
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Small_Bowel
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Terminal_Ileum
/nfs/data/COST_IBD/downstream_tasks/dea/scripts/gsea.R --input_dir /nfs/data/COST_IBD/downstream_tasks/dea/output/ibd_region_inflammation/Transverse_Colon


echo "GSEA analysis completed at $(date)"
