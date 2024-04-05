# Use r-base 4.3.2 as the base image
FROM r-base:4.3.2

ENTRYPOINT ["/bin/bash"]

# Install system dependencies for Rust, Julia, and Vim
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    libssl-dev \
    pkg-config \
    curl \
    vim \  
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Rust and Cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install the latest version of Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.2-linux-x86_64.tar.gz \
    && tar zxvf julia-1.10.2-linux-x86_64.tar.gz \ 
    && mv julia-1.10.2 /opt \
    && ln -s /opt/julia-1.10.2/bin/julia /usr/local/bin/julia

# Verify installations
RUN cargo --version && julia --version && vim --version
