#!/bin/bash

# Script to rebuild mamana/minimac4-all with multi-platform support
# GA4GH Hackathon 2025 - African Genomics Team

set -e

echo "🔧 Updating mamana/minimac4-all Docker image for multi-platform support..."

# Step 1: Set up Docker buildx for multi-platform builds
echo "📦 Setting up Docker buildx..."
docker buildx create --name multiplatform-builder --use || true
docker buildx inspect --bootstrap

# Step 2: Create improved Dockerfile with multi-platform support
cat > Dockerfile.multiplatform << 'EOF'
# Multi-platform Dockerfile for mamana/minimac4-all
# Supports both AMD64 and ARM64 architectures

FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:$PATH"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    cmake \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libhts-dev \
    tabix \
    bc \
    unzip \
    autoconf \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install htslib (architecture-aware build)
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 \
    && tar -xjf htslib-1.20.tar.bz2 \
    && cd htslib-1.20 \
    && autoreconf -i \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf htslib-1.20*

# Install bcftools (architecture-aware build)
RUN wget https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 \
    && tar -xjf bcftools-1.20.tar.bz2 \
    && cd bcftools-1.20 \
    && autoreconf -i \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf bcftools-1.20*

# Install Minimac4 (build from source for cross-platform compatibility)
RUN git clone https://github.com/statgen/Minimac4.git \
    && cd Minimac4 \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j$(nproc) \
    && cp minimac4 /usr/local/bin/ \
    && cd ../.. \
    && rm -rf Minimac4

# Clean up
RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/*

# Verify installations
RUN bcftools --version \
    && minimac4 --help \
    && which bcftools \
    && which minimac4

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]

# Metadata
LABEL maintainer="GA4GH Hackathon 2025 - African Genomics Team"
LABEL description="Multi-platform imputation tools: bcftools + Minimac4"
LABEL version="1.0.0-multiplatform"
LABEL platforms="linux/amd64,linux/arm64"
EOF

# Step 3: Build multi-platform image
echo "🏗️ Building multi-platform Docker image..."
docker buildx build \
    --platform linux/amd64,linux/arm64 \
    --file Dockerfile.multiplatform \
    --tag mamana/minimac4-all:multiplatform \
    --tag mamana/minimac4-all:latest \
    --push \
    .

echo "✅ Multi-platform build completed!"

# Step 4: Test the images on current platform
echo "🧪 Testing the updated image..."
docker run --rm mamana/minimac4-all:latest which bcftools
docker run --rm mamana/minimac4-all:latest which minimac4
docker run --rm mamana/minimac4-all:latest bcftools --version
docker run --rm mamana/minimac4-all:latest minimac4 --help | head -5

echo "🎉 Docker image mamana/minimac4-all updated successfully!"
echo "📋 Image now supports:"
echo "   ✅ linux/amd64 (Intel/AMD x86_64)"
echo "   ✅ linux/arm64 (Apple Silicon M1/M2/M3)"
echo ""
echo "🚀 Usage:"
echo "   docker pull mamana/minimac4-all:latest"
echo "   # No more --platform flags needed!"

# Step 5: Clean up
rm -f Dockerfile.multiplatform

EOF 