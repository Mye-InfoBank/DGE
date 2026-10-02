#!/usr/bin/env python3
"""
Simple differential gene expression (DGE) from an h5ad file.

Features:
- Reads an AnnData (.h5ad) file
- Optionally splits data by one metadata variable
- Runs one-vs-all DGE for each group of a specified metadata variable
- Outputs a single table with top N genes per group (default: 100)

Example:
python dge_h5ad.py \
  --input data.h5ad \
  --groupby condition \
  --split-by tissue \
  --output dge_top100.tsv
"""

from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import List

# Suppress anndata and scanpy deprecation warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import numpy as np 
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one-vs-all DGE on h5ad data and export top genes as a table."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input .h5ad file",
    )
    parser.add_argument(
        "--groupby",
        default=None,
        help="obs column used for one-vs-all differential expression",
    )
    parser.add_argument(
        "--groupby-template",
        default=None,
        help="Template for groupby column name (e.g. 'scanvi-{split_value}-2.0_leiden'). Used if --groupby is not provided.",
    )
    parser.add_argument(
        "--split-by",
        default=None,
        help="Optional obs column to split dataset before running DGE",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path for consolidated top DE genes",
    )
    parser.add_argument(
        "--n-top",
        type=int,
        default=100,
        help="Top N genes to keep per group (default: 100)",
    )
    parser.add_argument(
        "--method",
        choices=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"],
        default="wilcoxon",
        help="Method for scanpy.tl.rank_genes_groups (default: wilcoxon)",
    )
    parser.add_argument(
        "--layer",
        default=None,
        help="Optional layer name to use for DGE",
    )
    parser.add_argument(
        "--use-raw",
        action="store_true",
        help="Use adata.raw for DGE",
    )
    parser.add_argument(
        "--min-cells-per-group",
        type=int,
        default=3,
        help="Minimum cells required per group in a subset (default: 3)",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Number of parallel jobs for DGE computation (default: -1, use all CPUs)",
    )
    return parser.parse_args()


def ensure_categorical(adata: sc.AnnData, col: str) -> None:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in adata.obs")
    if not pd.api.types.is_categorical_dtype(adata.obs[col]):
        adata.obs[col] = adata.obs[col].astype("category")


def get_split_values(adata: sc.AnnData, split_by: str | None) -> List[str]:
    if split_by is None:
        return ["__ALL__"]
    ensure_categorical(adata, split_by)
    values = adata.obs[split_by].cat.categories.tolist()
    values = [v for v in values if pd.notna(v)]
    return values


def subset_adata(adata: sc.AnnData, split_by: str | None, split_value: str) -> sc.AnnData:
    if split_by is None:
        return adata.copy()
    mask = adata.obs[split_by] == split_value
    return adata[mask].copy()


def _compute_dge_for_group(
    adata: sc.AnnData,
    group: str,
    groupby: str,
    n_top: int,
    method: str,
    layer: str | None,
    use_raw: bool,
) -> pd.DataFrame:
    """Compute DGE for a single group (one-vs-all) in parallel."""
    group_mask = adata.obs[groupby] == group
    other_mask = ~group_mask

    group_data = adata[group_mask]
    other_data = adata[other_mask]

    if group_data.n_obs < 2 or other_data.n_obs < 2:
        return pd.DataFrame()

    adata_temp = adata.copy()
    adata_temp.obs["_temp_group"] = "other"
    adata_temp.obs.loc[group_mask, "_temp_group"] = group
    adata_temp = adata_temp[adata_temp.obs["_temp_group"].isin([group, "other"])].copy()

    sc.tl.rank_genes_groups(
        adata_temp,
        groupby="_temp_group",
        method=method,
        use_raw=use_raw,
        layer=layer,
        rankby_abs=True,
    )

    df_group = sc.get.rank_genes_groups_df(adata_temp, group=group)
    if df_group.empty:
        return pd.DataFrame()

    df_group = df_group.head(n_top).copy()
    df_group.insert(0, "rank", np.arange(1, len(df_group) + 1))
    df_group.insert(0, "group", group)

    return df_group


def run_dge_one_subset(
    adata_subset: sc.AnnData,
    groupby: str,
    n_top: int,
    method: str,
    layer: str | None,
    use_raw: bool,
    min_cells_per_group: int,
    logger: logging.Logger,
    n_jobs: int = -1,
) -> pd.DataFrame:
    ensure_categorical(adata_subset, groupby)

    counts = adata_subset.obs[groupby].value_counts()
    valid_groups = counts[counts >= min_cells_per_group].index.tolist()

    if len(valid_groups) < 2:
        return pd.DataFrame()

    adata_subset = adata_subset[adata_subset.obs[groupby].isin(valid_groups)].copy()
    adata_subset.obs[groupby] = adata_subset.obs[groupby].astype("category")

    logger.info(
        f"    Running parallel DGE on {len(valid_groups)} groups with method={method} "
        f"(n_jobs={n_jobs})"
    )

    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(_compute_dge_for_group)(
            adata_subset,
            group,
            groupby,
            n_top,
            method,
            layer,
            use_raw,
        )
        for group in valid_groups
    )

    rows = []
    for group, df_result in zip(valid_groups, results):
        if not df_result.empty:
            rows.append(df_result)
        else:
            logger.debug(f"      Group '{group}': no results")

    if not rows:
        return pd.DataFrame()

    result = pd.concat(rows, ignore_index=True)
    return result


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(__name__)

    logger.info("Starting DGE analysis")
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if args.groupby is None and args.groupby_template is None:
        raise ValueError("Must provide either --groupby or --groupby-template")

    logger.info(f"Loading h5ad file: {input_path}")
    adata = sc.read_h5ad(input_path)
    logger.info(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

    if args.split_by is not None:
        ensure_categorical(adata, args.split_by)
    elif args.groupby is not None:
        ensure_categorical(adata, args.groupby)

    split_values = get_split_values(adata, args.split_by)
    logger.info(f"Processing {len(split_values)} split(s)")

    all_results = []

    for i, split_value in enumerate(split_values, 1):
        logger.info(f"[{i}/{len(split_values)}] Processing split: {split_value}")
        adata_subset = subset_adata(adata, args.split_by, split_value)

        if adata_subset.n_obs == 0:
            logger.warning(f"No cells in split '{split_value}', skipping")
            continue

        logger.info(f"  Subset size: {adata_subset.n_obs} cells")

        groupby_col = args.groupby
        if args.groupby_template is not None:
            groupby_col = args.groupby_template.format(split_value=split_value)

        logger.info(f"  Grouping by: {groupby_col}")
        ensure_categorical(adata_subset, groupby_col)

        df = run_dge_one_subset(
            adata_subset=adata_subset,
            groupby=groupby_col,
            n_top=args.n_top,
            method=args.method,
            layer=args.layer,
            use_raw=args.use_raw,
            min_cells_per_group=args.min_cells_per_group,
            logger=logger,
            n_jobs=args.n_jobs,
        )

        if df.empty:
            logger.warning(f"  No valid DGE results for split '{split_value}'")
            continue

        logger.info(f"  Found {len(df)} top genes across groups")
        df.insert(0, "split_by", args.split_by if args.split_by is not None else "none")
        df.insert(1, "split_value", split_value if args.split_by is not None else "all")

        # Prepare output with preferred column order
        preferred_cols = [
            "split_by",
            "split_value",
            "group",
            "rank",
            "names",
            "scores",
            "logfoldchanges",
            "pvals",
            "pvals_adj",
            "pct_nz_group",
            "pct_nz_reference",
        ]
        remaining_cols = [c for c in df.columns if c not in preferred_cols]
        df = df[[c for c in preferred_cols if c in df.columns] + remaining_cols]

        # Create split-specific output filename
        if args.split_by is not None:
            output_base = output_path.stem
            output_ext = "".join(output_path.suffixes)
            split_value_safe = str(split_value).replace("/", "_").replace(" ", "_")
            split_output_path = output_path.parent / f"{output_base}_{split_value_safe}{output_ext}"
        else:
            split_output_path = output_path

        split_output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(split_output_path, sep="\t", index=False)
        logger.info(f"  Saved {len(df)} rows to {split_output_path}")
        all_results.append(split_output_path)

    if not all_results:
        logger.error("No valid DGE results produced. Check metadata columns and group sizes.")
        sys.exit(1)

    logger.info(f"Completed DGE analysis. Generated {len(all_results)} output file(s).")


if __name__ == "__main__":
    main()
