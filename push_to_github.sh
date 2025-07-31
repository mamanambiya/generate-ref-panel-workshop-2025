#!/bin/bash

# GitHub Container Registry Push Script  
# GA4GH Hackathon 2025 - African Genomics Team

set -e

# Configuration - UPDATE THESE VALUES
GITHUB_USERNAME="mamanambiya"  # Replace with your GitHub username/organization
IMAGE_NAME="federated-imputation"
TAG="latest"

LOCAL_TAG="${IMAGE_NAME}:${TAG}"
GHCR_TAG="ghcr.io/${GITHUB_USERNAME}/${IMAGE_NAME}:${TAG}"

echo "🚀 Pushing to GitHub Container Registry"
echo "======================================="
echo "Local image: ${LOCAL_TAG}"
echo "GHCR image: ${GHCR_TAG}"
echo "Registry: ghcr.io"
echo ""

# Check if local image exists
if ! docker images --format "table {{.Repository}}:{{.Tag}}" | grep -q "^${LOCAL_TAG}$"; then
    echo "❌ Local image ${LOCAL_TAG} not found!"
    echo "🔧 Please run ./build_docker.sh first"
    exit 1
fi

# Check if user is logged into GitHub Container Registry
echo "🔐 Checking GitHub Container Registry login..."
if ! docker info | grep -q "Username.*${GITHUB_USERNAME}"; then
    echo "❌ Not logged into GitHub Container Registry!"
    echo ""
    echo "🔧 Please login using your GitHub Personal Access Token:"
    echo "  1. Create PAT with 'write:packages' scope at: https://github.com/settings/tokens"
    echo "  2. Login: echo \$GITHUB_PAT | docker login ghcr.io -u ${GITHUB_USERNAME} --password-stdin"
    echo "  3. Or interactive: docker login ghcr.io -u ${GITHUB_USERNAME}"
    echo ""
    read -p "Press Enter after logging in, or Ctrl+C to exit..."
fi

# Tag image for GitHub Container Registry
echo "🏷️  Tagging image for GitHub Container Registry..."
docker tag ${LOCAL_TAG} ${GHCR_TAG}

# Push to GitHub Container Registry
echo "📤 Pushing image to GitHub Container Registry..."
docker push ${GHCR_TAG}

echo ""
echo "✅ Successfully pushed to GitHub Container Registry!"
echo ""
echo "📋 Image Details:"
echo "  Registry: ghcr.io"
echo "  Username: ${GITHUB_USERNAME}"
echo "  Image: ${IMAGE_NAME}"
echo "  Tag: ${TAG}"
echo "  Full path: ${GHCR_TAG}"
echo ""
echo "🔧 Update your WDL file to use:"
echo "  docker: \"${GHCR_TAG}\""
echo ""
echo "📖 Verify the push:"
echo "  Visit: https://github.com/${GITHUB_USERNAME}?tab=packages"
echo "  Or check: docker pull ${GHCR_TAG}"
echo ""
echo "🌐 Make package public (optional):"
echo "  1. Go to: https://github.com/users/${GITHUB_USERNAME}/packages/container/${IMAGE_NAME}/settings"
echo "  2. Scroll to 'Danger Zone' → 'Change package visibility'"
echo "  3. Select 'Public' to allow public access"
echo ""
echo "🎉 Ready for deployment anywhere!" 