#!/bin/bash

## Federated Imputation Pipeline Execution Script
## GA4GH Hackathon 2025 - African Genomics Team

set -euo pipefail

echo "🧬 Federated Imputation Pipeline"
echo "GA4GH Hackathon 2025 - African Genomics Team"
echo ""

# Check if workflow file exists
if [[ ! -f "workflows/federated_imputation_pipeline.wdl" ]]; then
    echo "❌ Workflow file not found: workflows/federated_imputation_pipeline.wdl"
    exit 1
fi

# Default inputs
INPUTS=${1:-"inputs/test_local.json"}

if [[ ! -f "$INPUTS" ]]; then
    echo "❌ Input file not found: $INPUTS"
    echo "Usage: $0 [input_file.json]"
    exit 1
fi

echo "Workflow: workflows/federated_imputation_pipeline.wdl"
echo "Inputs: $INPUTS"
echo ""

# Check for Cromwell
if [[ ! -f "cromwell.jar" ]]; then
    echo "Downloading Cromwell..."
    wget -O cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
fi

# Pull container
echo "Checking container..."
docker pull mamana/imputation:minimac4-4.1.6

# Create output directory
mkdir -p outputs/$(date +%Y%m%d_%H%M%S)

echo "🚀 Starting workflow execution..."
java -jar cromwell.jar run workflows/federated_imputation_pipeline.wdl -i "$INPUTS"

echo "✅ Workflow completed!"
