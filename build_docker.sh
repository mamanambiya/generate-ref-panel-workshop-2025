#!/bin/bash

# Federated Genotype Imputation Pipeline - Docker Build Script
# GA4GH Hackathon 2025 - African Genomics Team

set -e

# Configuration
IMAGE_NAME="federated-imputation"
TAG="latest"
FULL_TAG="${IMAGE_NAME}:${TAG}"

echo "🚀 Building Federated Genotype Imputation Pipeline Docker Image"
echo "=================================================="
echo "Image: ${FULL_TAG}"
echo "Build context: $(pwd)"
echo ""

# Build the Docker image
echo "📦 Building Docker image..."
docker build -t ${FULL_TAG} .

echo ""
echo "✅ Docker build completed successfully!"
echo ""
echo "📋 Image Details:"
docker images ${IMAGE_NAME}

echo ""
echo "🧪 Testing container tools..."
echo "Testing bcftools:"
docker run --rm ${FULL_TAG} bcftools --version | head -2

echo ""
echo "Testing minimac4:"
docker run --rm ${FULL_TAG} minimac4 --help | head -5

echo ""
echo "Testing bgzip:"
docker run --rm ${FULL_TAG} bgzip -h | head -3

echo ""
echo "🎉 Container is ready for use!"
echo ""
echo "📖 Usage Examples:"
echo "  # Run interactively:"
echo "  docker run -it --rm ${FULL_TAG}"
echo ""
echo "  # Mount data directory:"
echo "  docker run -it --rm -v \$(pwd):/data ${FULL_TAG}"
echo ""
echo "  # For GitHub Container Registry:"
echo "  ./push_to_github.sh"
echo "  # Or manually:"
echo "  docker tag ${FULL_TAG} ghcr.io/USERNAME/${IMAGE_NAME}:${TAG}"
echo "  docker push ghcr.io/USERNAME/${IMAGE_NAME}:${TAG}"
echo ""
echo "  # Update WDL file to use:"
echo "  docker: \"${FULL_TAG}\""
echo "  # or"
echo "  docker: \"ghcr.io/USERNAME/${IMAGE_NAME}:${TAG}\"" 