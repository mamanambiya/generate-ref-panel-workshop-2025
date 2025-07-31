# Cloud Deployment Guide

## Prerequisites

- Azure Batch account configured (for Azure deployment)
- Azure Storage account for input/output files (for Azure deployment)  
- WDL execution engine (e.g., Cromwell on Azure, Terra, DNAstack)
- GitHub Container Registry access (for custom container builds)

## Container Registry Options

### GitHub Container Registry (Recommended)
Use pre-built container from GitHub Container Registry:
```
ghcr.io/mamanambiya/federated-imputation:latest
```

### Custom Container Build
To build and push your own container:
```bash
./build_docker.sh
./push_to_github.sh
```

## Step-by-Step Deployment

### 1. Upload Files to Azure Storage

Upload these files to your Azure Blob Storage container:
```
federated_imputation_pipeline.wdl              # Remote container (original)
federated_imputation_pipeline_github.wdl       # GitHub Container Registry
federated_imputation_pipeline_local.wdl        # Local container build
test_input_azure.json                          # Azure configuration
test_chr22_region_bgzip.vcf.gz                 # Test data
test_chr22_region_bgzip.vcf.gz.csi             # Test data index
```

Choose the appropriate WDL file based on your container preference.

### 2. Configure Input Paths

Edit `test_input_azure.json` to match your Azure environment:

```json
{
  "FederatedImputationPipeline.input_vcf": "$AZ_BATCH_TASK_WORKING_DIR/test_chr22_region_bgzip.vcf.gz"
}
```

Or use absolute Azure paths:
```json
{
  "FederatedImputationPipeline.input_vcf": "/mnt/batch/tasks/workitems/your-job/job-1/test_chr22_region_bgzip.vcf.gz"
}
```

### 3. Submit to Azure Batch

Using Cromwell on Azure:
```bash
# Option 1: GitHub Container Registry (Recommended)
java -jar cromwell.jar run \
  federated_imputation_pipeline_github.wdl \
  -i test_input_azure.json \
  --backend Azure

# Option 2: Original remote container  
java -jar cromwell.jar run \
  federated_imputation_pipeline.wdl \
  -i test_input_azure.json \
  --backend Azure

# Option 3: Custom local build
java -jar cromwell.jar run \
  federated_imputation_pipeline_local.wdl \
  -i test_input_azure.json \
  --backend Azure
```

### 4. Common Azure Issues

**Exit Code 10 Error:**
- File not found in Azure working directory
- Incorrect file paths in input JSON
- Docker container not accessible

**Solutions:**
1. Ensure all input files are uploaded to the correct Azure location
2. Use Azure environment variables in paths: `$AZ_BATCH_TASK_WORKING_DIR`
3. Verify Docker container `mamana/minimac4-all` is accessible
4. Check Azure Batch pool configuration

### 5. Output Files

Successful execution will generate:
- `test_federated_panel.msav` - Main output for federated imputation
- Various statistics and summary files
- All outputs will be in the Azure Batch task working directory

### 6. Download Results

Download output files from Azure Storage to complete the process.

## Azure-Specific Configuration

The WDL has been optimized for Azure by:
- Removing Apple Silicon-specific Docker args
- Using standard Docker runtime configuration
- Supporting Azure Batch file path patterns
- Providing Azure-specific input templates

## Support

For Azure-specific issues:
- Check Azure Batch task logs
- Verify file accessibility in Azure environment
- Ensure proper Azure permissions and networking 