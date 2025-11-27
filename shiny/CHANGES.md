# Shiny App Changes - November 27, 2025

## Support for Nested Directory Structures

### Summary
Modified the Shiny app to handle both flat and nested directory structures when reading DEA results.

### Changes Made

#### 1. `global.R` - Updated `list_celltype_dirs()` function
- **Previous behavior**: Only looked for directories directly under the result folder
- **New behavior**: 
  - Detects if directories contain DESeq2 result files (TSV files)
  - If a directory doesn't contain results directly, checks subdirectories
  - Returns paths like "Subfolder/CellType" for nested structures
  - Maintains backward compatibility with flat structures

#### 2. `global.R` - Updated `scan_result_overview()` function
- Enhanced to handle nested paths correctly
- Uses basename of cell type for classification (e.g., "B cell" from "Colon/B cell")
- Displays full path for user interface

#### 3. `server.R` - Updated file path resolution
- Changed `pseudobulk_metadata()` to use `selected_celltype_path()` instead of manually constructing paths
- Changed `celltype_params()` to use `selected_celltype_path()` instead of manually constructing paths
- Ensures consistent path handling throughout the application

#### 4. `README.md` - Documentation
- Added section explaining supported directory structures
- Documented both flat and nested patterns
- Provided examples of each structure

### Testing Recommendations

Test with:
1. **Flat structure**: `output/ibd_condition/` 
   - Direct cell type folders
2. **Nested structure**: `output/ibd_region_inflammation/`
   - Region folders containing cell type folders

### Example Usage

For a nested structure like:
```
output/
  ibd_region_inflammation/
    Colon/
      B cell/
        DESeq2_results_Colon.tsv
        volcano.png
        ...
      T cell/
        ...
    Cecum/
      B cell/
        ...
```

The app will:
- Detect "Colon/B cell", "Colon/T cell", "Cecum/B cell", etc.
- Display them in the UI with their full path names
- Allow selection and proper file loading

### Backward Compatibility
✅ Fully backward compatible with existing flat directory structures
✅ No changes required to existing data organization
✅ Automatic detection of structure type
