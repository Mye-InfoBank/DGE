# COST MyeInfoBank Differential Expression Analysis (DEA) 

A comprehensive package for performing differential expression analysis and visualizing results for the COST IBD project.

## 📦 Package Contents

- **Analysis Scripts**: DESeq2 and GSEA analysis pipelines
- **Shiny Application**: Interactive web interface for results exploration
- **Docker Support**: Containerized deployment for easy setup
- **Documentation**: Comprehensive guides and examples

## 🚀 Quick Start

### Option 1: Docker (Recommended)

```bash
# Clone the repository
git clone <your-repo-url>
cd cost-ibd-dea

# Start the Shiny application
cd shiny
docker-compose up -d

# Access the app at http://localhost:3838
```

### Option 2: Local Development

```bash
# Install R dependencies
Rscript install_packages.R

# Run DESeq2 analysis
./run_deseq2.sh

# Run GSEA analysis
./run_gsea.sh

# Start Shiny app
cd shiny
R -e "shiny::runApp()"
```

## 📁 Directory Structure

```
cost-ibd-dea/
├── deseq2.R              # DESeq2 differential expression analysis
├── gsea.R                # GSEA pathway enrichment analysis
├── run_deseq2.sh         # SLURM script for DESeq2
├── run_gsea.sh           # SLURM script for GSEA
├── dea_slurm.sh          # Combined SLURM submission script
├── install_packages.R    # R package installation script
├── output/               # Analysis results (not in repo)
│   ├── celltype1/
│   │   └── condition1_vs_condition2/
│   │       ├── DESeq2_results.tsv
│   │       ├── gsea_results.tsv
│   │       └── gsea_plots/
│   └── logs/
└── shiny/                # Shiny application
    ├── app.R
    ├── global.R
    ├── ui.R
    ├── server.R
    ├── Dockerfile
    ├── docker-compose.yml
    └── README.md
```

## 🔬 Analysis Scripts

### DESeq2 Analysis (`deseq2.R`)

Performs differential expression analysis using DESeq2:

```r
# Key features:
- Multiple cell type support
- Flexible contrast definitions
- Quality control metrics
- Results export in multiple formats
```

**Usage:**
```bash
Rscript deseq2.R --help
```

### GSEA Analysis (`gsea.R`)

Performs Gene Set Enrichment Analysis:

```r
# Key features:
- MSigDB pathway analysis
- Multiple gene set collections
- Enrichment plot generation
- Results visualization
```

**Usage:**
```bash
Rscript gsea.R --help
```

### SLURM Scripts

- `run_deseq2.sh`: Submit DESeq2 analysis to SLURM cluster
- `run_gsea.sh`: Submit GSEA analysis to SLURM cluster
- `dea_slurm.sh`: Combined submission script

## 🖥️ Shiny Application

### Features

- **Interactive Results Browser**: Explore DEA results from multiple cell types
- **Volcano Plots**: Customizable plots with gene labeling
- **GSEA Visualization**: Interactive pathway analysis
- **Data Tables**: Sortable and searchable results
- **Progress Tracking**: Real-time loading indicators

### UI Components

#### Overview Page
- Cell type selection cards (color-coded by major cell type)
- Summary statistics and KPIs
- DESeq2 parameter display
- Interactive heatmaps

#### Detailed Analysis
- Volcano plots with configurable options
- GSEA results with pathway selection
- Data tables with sorting and filtering

### Docker Deployment

```bash
cd shiny
docker-compose up -d
```

## ⚙️ Configuration

### Environment Variables

```bash
# R settings
export R_MAX_MEM_SIZE=4G
export R_NUM_THREADS=4

# Analysis parameters
export DESEQ2_ALPHA=0.05
export GSEA_PVALUE=0.05
```

### Data Requirements

The analysis expects input data in the following format:

```
input_data/
├── counts_matrix.tsv     # Gene expression counts
├── metadata.tsv         # Sample metadata
└── gene_annotations.tsv # Gene information
```

## 🐳 Docker Support

### Shiny App Container

```dockerfile
# Features:
- Ubuntu 22.04 base
- R 4.3+ with all required packages
- Shiny application server
- Data volume mounting support
```

### Build and Run

```bash
# Build image
docker build -t cost-ibd-dea-shiny ./shiny

# Run container
docker run -d \
  --name dea-shiny-app \
  -p 3838:3838 \
  -v /path/to/results:/data:ro \
  cost-ibd-dea-shiny
```

## 📊 Expected Output Structure

```
output/
├── celltype1/
│   ├── condition1_vs_condition2/
│   │   ├── DESeq2_results.tsv
│   │   ├── gsea_results.tsv
│   │   ├── volcano_plot.png
│   │   └── gsea_plots/
│   │       ├── HALLMARK_ADIPOGENESIS.png
│   │       └── ...
│   └── condition3_vs_condition4/
├── celltype2/
└── logs/
    ├── DESeq2_20240101_120000.log
    └── GSEA_20240101_140000.log
```

## 🔧 Customization

### Analysis Parameters

Modify analysis parameters in the R scripts:

```r
# DESeq2 parameters
alpha <- 0.05
lfc_threshold <- 0.5

# GSEA parameters
pvalue_cutoff <- 0.05
min_gene_set_size <- 15
max_gene_set_size <- 500
```

### Color Schemes

Update color palettes in `shiny/global.R`:

```r
MAJOR_CLASS_PALETTE <- c(
  "T/NK" = "#FF6B6B",
  "Myeloid" = "#4ECDC4",
  "B" = "#45B7D1",
  "Stromal" = "#96CEB4",
  "Other" = "#FFEAA7"
)
```

## 🐛 Troubleshooting

### Common Issues

1. **Package Installation**
   ```bash
   # Reinstall packages
   Rscript install_packages.R
   ```

2. **Memory Issues**
   ```bash
   # Increase memory limit
   export R_MAX_MEM_SIZE=8G
   ```

3. **Docker Issues**
   ```bash
   # Check logs
   docker logs dea-shiny-app
   
   # Rebuild image
   docker build --no-cache -t cost-ibd-dea-shiny ./shiny
   ```

### Logs and Debugging

```bash
# Analysis logs
tail -f output/logs/DESeq2_*.log
tail -f output/logs/GSEA_*.log

# Shiny app logs
docker logs -f dea-shiny-app
```

## 📝 Usage Examples

### Running Complete Analysis

```bash
# 1. Submit DESeq2 analysis
sbatch run_deseq2.sh

# 2. Submit GSEA analysis (after DESeq2 completes)
sbatch run_gsea.sh

# 3. Start Shiny app to explore results
cd shiny
docker-compose up -d
```

### Custom Analysis

```r
# Load functions
source("deseq2.R")
source("gsea.R")

# Run custom analysis
run_deseq2_analysis(
  counts_file = "my_counts.tsv",
  metadata_file = "my_metadata.tsv",
  output_dir = "my_results"
)
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📞 Support

For issues and questions:
- Create an issue in the GitHub repository
- Check the troubleshooting section
- Review the application logs

## 📄 License

[Add your license information here]

---

**Note**: This package is designed specifically for the COST IBD project differential expression analysis. Ensure your data follows the expected structure for optimal functionality. 