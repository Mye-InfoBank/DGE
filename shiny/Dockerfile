# Use the official R base image with Ubuntu 22.04
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV R_VERSION=4.3.2
ENV CRAN_MIRROR=https://cloud.r-project.org/

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    gnupg2 \
    software-properties-common \
    apt-transport-https \
    ca-certificates \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    libgsl-dev \
    libblas-dev \
    liblapack-dev \
    libopenblas-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    libreadline-dev \
    zlib1g-dev \
    libgfortran5 \
    libquadmath0 \
    gfortran \
    gcc \
    g++ \
    make \
    cmake \
    pkg-config \
    git \
    && rm -rf /var/lib/apt/lists/*

# Add R repository and install R (modern GPG key method)
RUN wget -O- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor > /usr/share/keyrings/cran-archive-keyring.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" > /etc/apt/sources.list.d/r-project.list \
    && apt-get update \
    && apt-get install -y r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
RUN mkdir -p /app

# Copy R package installation script
COPY install_packages.R /tmp/install_packages.R

# Install R packages
RUN Rscript /tmp/install_packages.R

# Copy app files
COPY app.R /app/
COPY global.R /app/
COPY ui.R /app/
COPY server.R /app/

# Create data directory for mounting
RUN mkdir -p /data

# Set permissions
RUN chown -R root:root /app

# Expose port
EXPOSE 3838

# Create startup script
RUN echo '#!/bin/bash\n\
# Update OUTPUT_ROOT to use mounted data directory\n\
sed -i "s|OUTPUT_ROOT <- \".*\"|OUTPUT_ROOT <- \"/data\"|g" /app/global.R\n\
\n\
# Start R Shiny app\n\
cd /app\n\
R -e "shiny::runApp(host = \"0.0.0.0\", port = 3838)"' > /start.sh \
    && chmod +x /start.sh

# Set working directory
WORKDIR /app

# Start R Shiny app
CMD ["/start.sh"] 