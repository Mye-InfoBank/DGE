# DEA Shiny App Docker Setup

This directory contains the Docker configuration for the Differential Expression Analysis (DEA) Shiny application.

## Prerequisites

- Docker
- Docker Compose (optional, for easier deployment)

## Quick Start

### Using Docker Compose (Recommended)

1. **Build and run the application:**
   ```bash
   docker-compose up --build
   ```

2. **Access the application:**
   Open your browser and navigate to `http://localhost:3838`

3. **Stop the application:**
   ```bash
   docker-compose down
   ```

### Using Docker directly

1. **Build the Docker image:**
   ```bash
   docker build -t dea-shiny-app .
   ```

2. **Run the container:**
   ```bash
   docker run -d \
     --name dea-shiny-app \
     -p 3838:3838 \
     -v /nfs/data/COST_IBD/downstream_tasks/dea/output:/data:ro \
     dea-shiny-app
   ```

3. **Access the application:**
   Open your browser and navigate to `http://localhost:3838`

4. **Stop and remove the container:**
   ```bash
   docker stop dea-shiny-app
   docker rm dea-shiny-app
   ```

## Configuration

### Data Mounting

The application expects DEA output data to be mounted at `/data` inside the container. The default configuration mounts:
- `/nfs/data/COST_IBD/downstream_tasks/dea/output:/data:ro`

You can modify the volume mapping in `docker-compose.yml` or the `docker run` command to point to your actual data directory.

#### Supported Directory Structures

The application supports two directory structures:

1. **Flat structure** (cell types directly under result directory):
   ```
   output/
     ibd_condition/
       B cell/
       T cell/
       ...
   ```

2. **Nested structure** (subfolders containing cell types):
   ```
   output/
     ibd_region_inflammation/
       Colon/
         B cell/
         T cell/
       Cecum/
         B cell/
         ...
   ```

The app automatically detects the structure and displays cell types appropriately. For nested structures, cell types are shown as "Subfolder/Cell Type" (e.g., "Colon/B cell").

### Environment Variables

You can customize the following environment variables:

- `R_MAX_MEM_SIZE`: Maximum memory for R (default: 4G)
- `R_MAX_NUM_THREADS`: Maximum number of threads for R (default: 4)

### Port Configuration

The application runs on port 3838 by default. You can change this by modifying the port mapping in `docker-compose.yml` or the `docker run` command.

## Included Packages

The Docker image includes all necessary R packages for the DEA Shiny app:

### Core Shiny Packages
- `shiny`, `bslib`, `shinycssloaders`, `DT`

### Data Processing
- `data.table`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`

### Visualization
- `ggplot2`, `ggrepel`, `plotly`, `viridis`, `RColorBrewer`, `scales`
- `gridExtra`, `cowplot`, `patchwork`

### Bioconductor Packages
- `msigdbr`, `clusterProfiler`, `enrichplot`
- `org.Hs.eg.db`, `AnnotationDbi`
- Various single-cell analysis packages

### Additional Utilities
- `base64enc`, `future`, `jsonlite`, `yaml`, `fs`
- Development tools: `roxygen2`, `devtools`, `testthat`

## Troubleshooting

### Build Issues

If you encounter build issues:

1. **Memory issues during package installation:**
   ```bash
   docker build --memory=8g --memory-swap=8g -t dea-shiny-app .
   ```

2. **Network issues:**
   Ensure you have a stable internet connection for downloading R packages.

### Runtime Issues

1. **Application not accessible:**
   - Check if the container is running: `docker ps`
   - Check container logs: `docker logs dea-shiny-app`
   - Verify port mapping: `docker port dea-shiny-app`

2. **Data not loading:**
   - Verify the data directory is correctly mounted
   - Check file permissions on the host system
   - Ensure the data directory contains the expected structure

3. **Performance issues:**
   - Increase memory allocation: `R_MAX_MEM_SIZE=8G`
   - Increase thread count: `R_MAX_NUM_THREADS=8`
   - Consider using a more powerful host system

## Development

### Modifying the Application

1. Make changes to the R files (`app.R`, `global.R`, `ui.R`, `server.R`)
2. Rebuild the Docker image: `docker-compose build`
3. Restart the container: `docker-compose up -d`

### Adding New Packages

1. Edit `install_packages.R` to add new package names
2. Rebuild the Docker image: `docker-compose build --no-cache`

## Security Notes

- The application runs as the `shiny` user inside the container
- Data volumes are mounted as read-only (`:ro`) for security
- The container includes only necessary packages to minimize attack surface

## Support

For issues related to:
- Docker setup: Check Docker logs and configuration
- Application functionality: Review the Shiny app code and R package documentation
- Data access: Verify file permissions and directory structure 