import scanpy as sc
import numpy as np
import pandas as pd
from typing import List, Optional, Tuple
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:
    raise ImportError("Please install pyDESeq2: pip install pyDESeq2")

pseudobulk_group_cols=["sample", "annotation(old):fine"]
min_cells_in_group=3
pseudobulk_func="sum"


adata = sc.read_h5ad('/nfs/proj/COST_IBD/alex_dev/DGE/b_macro_subset_large.h5ad')


grouping_key = "pseudobulk_group_key"
while grouping_key in adata.obs.columns:
    grouping_key += "_dup"
adata.obs[grouping_key] = (
    adata.obs[pseudobulk_group_cols]
    .astype(str)
    .apply(lambda x: "--".join(x), axis=1)
)
print(f"Created grouping key '{grouping_key}' in adata.obs.")
print("Example grouping key values:")
print(adata.obs[grouping_key].head(5))


# 2. Convert to a cell-level counts DataFrame
print("Step 2: Creating cell_counts_df from adata.X...")
if not hasattr(adata.X, "toarray"):
    counts_matrix = adata.layers['raw']  # already dense
    print("Data is already dense.")
else:
    counts_matrix = adata.layers['raw'].toarray()  # convert sparse to dense
    print("Converted sparse matrix to dense.")

# Round only the non-integer values to the next highest integer
non_integer_mask = (counts_matrix % 1 != 0)
counts_matrix[non_integer_mask] = np.ceil(counts_matrix[non_integer_mask])

# Verify that all values are now integers
assert np.all(counts_matrix % 1 == 0), "There are still non-integer values!"
print("Non-integer values rounded up successfully.")

gene_names = adata.var_names.tolist()
cell_counts_df = pd.DataFrame(
    counts_matrix,
    columns=gene_names,
    index=adata.obs.index
)
cell_counts_df[grouping_key] = adata.obs[grouping_key].values
print(f"cell_counts_df shape: {cell_counts_df.shape}")
print("---------------------------------------------------")


# 3. Group and aggregate
print("Step 3: Grouping and aggregating...")
grouped = cell_counts_df.groupby(grouping_key)
group_sizes = grouped.size()
valid_groups = group_sizes[group_sizes >= min_cells_in_group].index
print(f"Number of total groups: {len(group_sizes)}")
print(f"Number of valid groups (>= {min_cells_in_group} cells): {len(valid_groups)}")
if len(valid_groups) == 0:
    raise ValueError("No groups found with at least min_cells_in_group cells.")
cell_counts_df = cell_counts_df[cell_counts_df[grouping_key].isin(valid_groups)]
grouped = cell_counts_df.groupby(grouping_key)
if not hasattr(pd.core.groupby.generic.DataFrameGroupBy, pseudobulk_func):
    raise ValueError(f"Invalid pseudobulk_func '{pseudobulk_func}'. "
                        "Must be a valid pandas groupby agg method.")
pseudobulk_df = grouped.aggregate(pseudobulk_func)
print("Aggregated pseudobulk dataframe shape:", pseudobulk_df.shape)
print("---------------------------------------------------")

# 3.1 Also store the actual number of cells in each group
n_cells_per_group = grouped.size()
if grouping_key in pseudobulk_df.columns:
    pseudobulk_df.drop(columns=grouping_key, inplace=True, errors="ignore")


# 4. Build pseudobulk-level metadata
print("Step 4: Building metadata for pseudobulk samples...")
meta_list = []
for idx in pseudobulk_df.index:
    split_vals = idx.split("--")
    meta_dict = {}
    for col, val in zip(pseudobulk_group_cols, split_vals):
        meta_dict[col] = val
    meta_list.append(meta_dict)
pseudo_metadata = pd.DataFrame(meta_list, index=pseudobulk_df.index)
print(f"pseudo_metadata shape: {pseudo_metadata.shape}")
print("pseudo_metadata head:")
print(pseudo_metadata.head(3))
print("---------------------------------------------------")

# 4.1 Add number of cells and log2 number of cells
pseudo_metadata["n_cells"] = n_cells_per_group.reindex(pseudo_metadata.index)
pseudo_metadata["log2_ncells"] = np.log2(pseudo_metadata["n_cells"].astype(float))
print("Added 'n_cells' and 'log2_ncells' to pseudo_metadata.")
print(pseudo_metadata[["n_cells", "log2_ncells"]].head(3))
print("---------------------------------------------------")


# 5 create DESeq2 object
dds = DeseqDataSet(
    counts=pseudobulk_df,
    metadata=pseudo_metadata,
    design_factors = 'annotation(old):fine',
    ref_level=['annotation(old):fine', 'B_cells']
)

pseudobulk_df.to_csv('pseudobulk_df.csv')
pseudo_metadata.to_csv('pseudo_metadata.csv')
