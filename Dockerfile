FROM rocker/r-ver:4.3.2

# System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install FarmCPUpp R dependencies
RUN R -e "install.packages(c('remotes', 'Rcpp', 'data.table', 'optparse', 'vcfR', 'bigmemory', 'ggplot2', 'qqman'), repos='https://cloud.r-project.org/')"

# Install FarmCPUpp
RUN R -e "remotes::install_github('amkusmec/FarmCPUpp')"

# Set working directory and copy your app
WORKDIR /usr/src/app
COPY . /usr/src/app

CMD ["R"]