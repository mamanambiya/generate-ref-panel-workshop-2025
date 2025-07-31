# GitHub Container Registry Guide

## Overview

This guide explains how to use GitHub Container Registry (ghcr.io) with the Federated Genotype Imputation Pipeline for the GA4GH Hackathon 2025.

## Why GitHub Container Registry?

- **Free for public repositories**
- **Integrated with GitHub** - packages appear alongside your code
- **No vendor lock-in** - works with any Docker-compatible platform
- **Better for open science** - easier sharing and collaboration
- **Azure compatible** - works seamlessly with Azure Batch

## Pre-built Container

### Using the Public Container

The pipeline is available as a pre-built container:
```
ghcr.io/mamanambiya/federated-imputation:latest
```

**Usage with WDL:**
```bash
java -jar cromwell.jar run federated_imputation_pipeline_github.wdl -i test_input.json
```

### Container Contents
- **bcftools 1.20** - VCF manipulation and statistics
- **htslib 1.20** - bgzip, tabix for file compression/indexing  
- **minimac4 (latest)** - MSAV format conversion
- **Ubuntu 22.04** base with all genomics dependencies

## Building Your Own Container

### 1. Setup GitHub Personal Access Token

Create a Personal Access Token with `write:packages` scope:
1. Go to: https://github.com/settings/tokens
2. Click "Generate new token (classic)"
3. Select scopes: `write:packages`, `read:packages`
4. Copy the token

### 2. Login to GitHub Container Registry

```bash
# Export your token
export GITHUB_PAT=your_token_here

# Login to GitHub Container Registry
echo $GITHUB_PAT | docker login ghcr.io -u YOUR_USERNAME --password-stdin
```

### 3. Build and Push Container

```bash
# Build locally
./build_docker.sh

# Push to GitHub Container Registry
./push_to_github.sh
```

### 4. Update Configuration

Edit `push_to_github.sh` to use your GitHub username:
```bash
GITHUB_USERNAME="your-username"  # Replace with your username
```

## Usage Examples

### Local Testing
```bash
# Test with GitHub Container Registry version
java -jar cromwell.jar run federated_imputation_pipeline_github.wdl -i test_input.json
```

### Azure Batch Deployment
```bash
# Azure deployment with GitHub container
java -jar cromwell.jar run federated_imputation_pipeline_github.wdl \
  -i test_input_azure.json \
  --backend Azure
```

### Terra/AnVIL
1. Upload `federated_imputation_pipeline_github.wdl` to your workspace
2. Configure inputs pointing to your data
3. Run workflow - Terra will automatically pull the container

### DNAstack Workbench
1. Import `federated_imputation_pipeline_github.wdl`
2. Configure input JSON
3. Execute - DNAstack handles container pulling

## Container Visibility

### Making Container Public
To allow public access (recommended for open science):
1. Go to: https://github.com/users/YOUR_USERNAME/packages/container/federated-imputation/settings
2. Scroll to "Danger Zone" → "Change package visibility"
3. Select "Public"

### Private Container Access
For private containers, users need authentication:
```bash
docker login ghcr.io -u USERNAME
docker pull ghcr.io/USERNAME/federated-imputation:latest
```

## Troubleshooting

### Authentication Issues
```bash
# Check login status
docker info | grep Username

# Re-login if needed
docker logout ghcr.io
echo $GITHUB_PAT | docker login ghcr.io -u USERNAME --password-stdin
```

### Container Not Found
- Verify the container name: `ghcr.io/USERNAME/federated-imputation:latest`
- Check if the container is public or if you need authentication
- Ensure you have read access to the repository

### Push Permission Denied
- Verify your Personal Access Token has `write:packages` scope
- Check that you're logged in as the correct user
- Ensure you have push permissions to the repository

## Workflow File Comparison

| File | Container Source | Use Case |
|------|------------------|----------|
| `federated_imputation_pipeline.wdl` | Remote (mamana/minimac4-all) | Original/legacy |
| `federated_imputation_pipeline_github.wdl` | GitHub Registry (ghcr.io) | **Recommended** |
| `federated_imputation_pipeline_local.wdl` | Local build | Development |

## Benefits for GA4GH Hackathon 2025

- **Easy sharing** - participants can easily access the same container
- **Version control** - container builds are linked to code versions  
- **Reproducibility** - exact same environment across all deployments
- **Collaboration** - teams can contribute to container improvements
- **Open science** - transparent and accessible to the community

## Support

For GitHub Container Registry issues:
- GitHub Documentation: https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry
- Container Issues: Check the repository's Issues page
- Technical Support: GA4GH Hackathon 2025 Discord/Slack channels

Ready for federated genomics research! 🧬 