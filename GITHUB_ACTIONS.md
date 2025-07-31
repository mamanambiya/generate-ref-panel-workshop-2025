# GitHub Actions CI/CD

## Overview

This repository includes automated GitHub Actions workflows for building, testing, and deploying the Federated Genotype Imputation Pipeline for the GA4GH Hackathon 2025.

## Workflows

### 1. Build and Push Docker Image (`docker-build.yml`)

**Triggers:**
- Push to `main` or `develop` branches 
- Pull requests to `main` or `develop` branches
- Manual workflow dispatch
- New releases

**Features:**
- ✅ **Multi-platform builds**: linux/amd64, linux/arm64
- ✅ **GitHub Container Registry**: Automatic push to ghcr.io
- ✅ **Security scanning**: Trivy vulnerability scanning
- ✅ **Container testing**: Validates bcftools, minimac4, bgzip functionality
- ✅ **Automatic tagging**: latest, branch names, GA4GH-specific tags

**Container Output:**
```
ghcr.io/mamanambiya/federated-imputation:latest
ghcr.io/mamanambiya/federated-imputation:ga4gh-hackathon-2025
```

### 2. Test WDL Workflows (`test-wdl.yml`)

**Triggers:**
- Push/PR with changes to `*.wdl` or `*.json` files
- Manual workflow dispatch

**Features:**
- ✅ **WDL syntax validation**: All workflow files validated with Cromwell
- ✅ **JSON input validation**: Test input files checked for valid JSON
- ✅ **Multi-workflow testing**: Tests all WDL variants
- ✅ **Java 11 + Cromwell 85**: Latest stable versions

**Validated Files:**
- `federated_imputation_pipeline.wdl`
- `federated_imputation_pipeline_github.wdl` 
- `federated_imputation_pipeline_local.wdl`
- `federated_imputation_pipeline_clean.wdl`
- `test_input.json`
- `test_input_azure.json`

## Build Status Badges

The README includes status badges that show current build status:

[![Build and Push Docker Image](https://github.com/mamanambiya/generate-ref-panel-workshop-2025/actions/workflows/docker-build.yml/badge.svg)](https://github.com/mamanambiya/generate-ref-panel-workshop-2025/actions/workflows/docker-build.yml)
[![Test WDL Workflows](https://github.com/mamanambiya/generate-ref-panel-workshop-2025/actions/workflows/test-wdl.yml/badge.svg)](https://github.com/mamanambiya/generate-ref-panel-workshop-2025/actions/workflows/test-wdl.yml)

## Using Automated Builds

### For GA4GH Hackathon Participants

**Quick Start:**
```bash
# Pull the latest automated container
docker pull ghcr.io/mamanambiya/federated-imputation:latest

# Use with GitHub Registry WDL
java -jar cromwell.jar run federated_imputation_pipeline_github.wdl -i test_input.json
```

**Benefits:**
- 🚀 **Always up-to-date**: Latest builds on every commit
- 🔒 **Security scanned**: Automatic vulnerability detection  
- 🌍 **Multi-platform**: Works on Apple Silicon and x86_64
- ✅ **Tested**: All tools validated before publishing
- 📦 **Version tracked**: Git commit linked to container

### For Developers

**Development Workflow:**
1. Make changes to Dockerfile or WDL files
2. Push to `develop` branch
3. GitHub Actions automatically:
   - Builds new container
   - Validates all WDL files  
   - Runs security scans
   - Tests genomics tools
4. Check status badges for build results
5. Container available at `ghcr.io/mamanambiya/federated-imputation:develop`

## Manual Triggering

Both workflows support manual triggering:

1. Go to **Actions** tab in GitHub repository
2. Select the workflow (Docker Build or WDL Test)
3. Click **"Run workflow"**
4. Choose branch and click **"Run workflow"**

## Secrets and Permissions

**Required Permissions:**
- `contents: read` - Repository access
- `packages: write` - Push to GitHub Container Registry

**Automatic Secrets:**
- `GITHUB_TOKEN` - Automatically provided by GitHub
- No additional secrets needed!

## Troubleshooting

### Build Failures

**Common Issues:**
- **Docker build timeout**: Large genomics tools take time to compile
- **Platform compatibility**: ARM64 builds may take longer
- **Network issues**: Downloads of bcftools/minimac4 sources

**Solutions:**
- Check Actions logs for specific error messages
- Re-run failed jobs (often transient issues)
- Build locally to test changes before pushing

### Container Access Issues

**Making Container Public:**
1. Go to GitHub profile → Packages
2. Find `federated-imputation` package
3. Package settings → Change visibility → Public

### WDL Validation Errors

**Common Issues:**
- Syntax errors in WDL files
- Invalid JSON in input files
- Missing required parameters

**Solutions:**
- Test locally with Cromwell before pushing
- Validate JSON files with `python -m json.tool`
- Check WDL syntax with `womtool validate`

## Integration with GA4GH Hackathon 2025

**Participant Benefits:**
- ✅ **No setup required**: Pre-built containers ready to use
- ✅ **Reliable builds**: Automated testing ensures quality
- ✅ **Platform support**: Works on all common architectures
- ✅ **Version control**: Reproducible research with tagged versions
- ✅ **Security**: Regular vulnerability scans

**Workflow Integration:**
- Use `federated_imputation_pipeline_github.wdl` for automatic container pulling
- Reference `ghcr.io/mamanambiya/federated-imputation:latest` in custom workflows
- Fork repository to customize for specific research needs

## Support

For GitHub Actions issues:
- Check workflow logs in Actions tab
- Review this documentation
- Open issue in repository
- Contact GA4GH Hackathon organizers

Ready for automated genomics research! 🧬🚀 