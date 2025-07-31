# Federated Genotype Imputation Pipeline - Docker Container
# GA4GH Hackathon 2025 - African Genomics Team

FROM ubuntu:22.04

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    cmake \
    git \
    autoconf \
    automake \
    libtool \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    bc \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install bcftools (latest version)
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 \
    && tar -xjf bcftools-1.20.tar.bz2 \
    && cd bcftools-1.20 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf bcftools-1.20*

# Install htslib (for bgzip, tabix)
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 \
    && tar -xjf htslib-1.20.tar.bz2 \
    && cd htslib-1.20 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf htslib-1.20*

# Install savvy library (required by Minimac4)
RUN git clone https://github.com/statgen/savvy.git \
    && cd savvy \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make \
    && make install \
    && cd ../.. \
    && rm -rf savvy

# Install Minimac4
RUN git clone https://github.com/statgen/Minimac4.git \
    && cd Minimac4 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && cp minimac4 /usr/local/bin/ \
    && cd ../.. \
    && rm -rf Minimac4

# Update library paths
RUN ldconfig

# Add /usr/local/bin to PATH
ENV PATH="/usr/local/bin:${PATH}"

# Verify installations
RUN bcftools --version && \
    minimac4 --help && \
    bgzip -h && \
    tabix -h

# Create working directories
RUN mkdir -p /data /output

# Set default working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]

# Metadata
LABEL maintainer="GA4GH Hackathon 2025 - African Genomics Team"
LABEL description="Federated Genotype Imputation Pipeline Container"
LABEL version="1.0.0"
LABEL tools="bcftools-1.20,minimac4-latest,htslib-1.20" 