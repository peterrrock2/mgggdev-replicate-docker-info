# Use r-base 4.3.2 as the base image
FROM r-base:4.3.2

ENTRYPOINT ["/bin/bash"]

# Install system dependencies for Rust, Julia, and Vim
RUN apt-get update && apt-get install -y \
    time \
    git \
    wget \
    build-essential \
    libssl-dev \
    pkg-config \
    curl \
    vim \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Rust and Cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install the latest version of Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.2-linux-x86_64.tar.gz \
    && tar zxvf julia-1.10.2-linux-x86_64.tar.gz \
    && rm julia-1.10.2-linux-x86_64.tar.gz \ 
    && mv julia-1.10.2 /opt \
    && ln -s /opt/julia-1.10.2/bin/julia /usr/local/bin/julia


# Install base R dependencies
RUN R -e "install.packages(c('argparser', 'dplyr', 'ggplot2', 'remotes', 'sf'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
# Install redist. This is split to make cached rebuilds easer for the dev
RUN R -e "remotes::install_version('redist', version='4.2.0', repos='http://cran.rstudio.com/')"

RUN cargo install --git https://github.com/peterrrock2/msms_parser.git --tag=v0.1.0
RUN cargo install --git https://github.com/peterrrock2/smc_parser.git --tag=v0.1.0
RUN cargo install binary-ensemble --version 0.1.2

# Install the necessary runner files every time
ARG CACHEBUST=1
COPY ./home /home

# # Set up forest recom
RUN cd /home/forest \
    && julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.build();'
RUN julia -e 'using Pkg; Pkg.add("ArgParse"); Pkg.develop(PackageSpec(path="/home/forest"))'

RUN cd /home/frcw \
    && cargo install --path .