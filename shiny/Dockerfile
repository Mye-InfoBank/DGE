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
    && rm -rf /var/lib/apt/lists/*

# Copy R package installation script
COPY install_packages.R /tmp/install_packages.R

# Install R packages — fail the build if any package is missing
RUN Rscript /tmp/install_packages.R && \
    Rscript -e "stopifnot(requireNamespace('shiny'), requireNamespace('DT'), requireNamespace('bslib'))"

# Copy app files
COPY app.R /app/
COPY global.R /app/
COPY ui.R /app/
COPY server.R /app/

RUN mkdir -p /data && chown -R root:root /app

EXPOSE 3838

RUN echo '#!/bin/bash\n\
sed -i "s|OUTPUT_ROOT <- \".*\"|OUTPUT_ROOT <- \"/data\"|g" /app/global.R\n\
cd /app\n\
R -e "shiny::runApp(host = \"0.0.0.0\", port = 3838)"' > /start.sh \
    && chmod +x /start.sh

WORKDIR /app
CMD ["/start.sh"]
