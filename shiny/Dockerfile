FROM rocker/shiny:4.4.2

# Use Posit Package Manager for pre-built Ubuntu binaries
RUN echo 'options(repos = c(PPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))' \
    >> /usr/local/lib/R/etc/Rprofile.site

# Install system libraries needed by some packages
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
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
    libuv1-dev \
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

# Copy R package installation script
COPY install_packages.R /tmp/install_packages.R

# Install R packages
RUN Rscript /tmp/install_packages.R && Rscript -e "library(shiny); cat('shiny OK\n')"

# Copy app files
#COPY app.R /app/
#COPY global.R /app/
#COPY ui.R /app/
#COPY server.R /app/
#COPY ai_snapshots /app/ai_snapshots/
# Copy all R files used by the app
COPY *.R /app/

# Copy snapshot directory
COPY ai_snapshots /app/ai_snapshots/

RUN mkdir -p /data && chown -R root:root /app

EXPOSE 3838

RUN echo '#!/bin/bash\n\
sed -i "s|OUTPUT_ROOT <- \".*\"|OUTPUT_ROOT <- \"/data\"|g" /app/global.R\n\
cd /app\n\
R -e "shiny::runApp(host = \"0.0.0.0\", port = 3838)"' > /start.sh \
    && chmod +x /start.sh

WORKDIR /app
CMD ["/start.sh"]
