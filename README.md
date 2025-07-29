# Federated Genotype Imputation Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![WDL](https://img.shields.io/badge/WDL-1.0-blue.svg)](https://github.com/openwdl/wdl)
[![Docker](https://img.shields.io/badge/Docker-Ready-brightgreen.svg)](https://hub.docker.com/r/mamana/imputation)

## GA4GH Hackathon 2025 - African Genomics Team

A production-ready WDL workflow for processing VCF files into MSAV format to enable federated genotype imputation networks while preserving data sovereignty.

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start) 
- [Workflow Architecture](#workflow-architecture)
- [Installation](#installation)
- [Usage](#usage)
- [Testing](#testing)
- [Deployment](#deployment)
- [Contributing](#contributing)

## Overview

### What This Pipeline Does
1. **Extract genomic regions** from VCF files using bcftools
2. **Apply quality control** filters (MAF, call rate, variant type)
3. **Convert to MSAV format** using Minimac4 for federated imputation
4. **Generate comprehensive statistics** and validation reports

### Key Features
- **Data Sovereignty**: Raw genomic data never leaves your institution
- **African Genomics Focus**: Optimized for H3Africa and AWI-Gen reference panels
- **GA4GH Compliant**: Standard WDL format for cross-platform execution
- **Production Ready**: Comprehensive error handling and validation
- **Scalable**: Containerized approach supports large datasets

## Quick Start

### 1. Prerequisites
- Docker or Singularity
- WDL execution engine (Cromwell, DNAstack Workbench, Terra, etc.)
- Input VCF file (bgzip compressed and indexed)

### 2. Run Test Pipeline
```bash
# Pull the container
docker pull mamana/imputation:minimac4-4.1.6

# Run test workflow
java -jar cromwell.jar run workflows/federated_imputation_pipeline.wdl -i inputs/test_local.json
```

### 3. Expected Output
- `federated_panel.msav` - Ready for federated imputation networks
- Quality-controlled VCF files with statistics
- Comprehensive processing reports

## Workflow Architecture

### Modular Design
```
workflows/
├── federated_imputation_pipeline.wdl    # Main orchestration workflow
│
tasks/
├── extract_region.wdl                   # Genomic region extraction
├── quality_control.wdl                  # VCF filtering and QC
└── minimac_conversion.wdl               # MSAV format conversion
```

### Pipeline Flow
```mermaid
graph LR
    A[Input VCF] --> B[Extract Region]
    B --> C[Quality Control]
    C --> D[MSAV Conversion]
    D --> E[Federated Reference Panel]
```

### Data Flow & Sovereignty
- **Input**: VCF files remain in your secure environment
- **Processing**: All computation happens within your infrastructure
- **Output**: MSAV files enable collaboration without data sharing
- **Result**: Federated imputation while maintaining data sovereignty

## Installation

### Container Requirements
The pipeline uses the validated container: `mamana/imputation:minimac4-4.1.6`

**Contains:**
- Minimac4 v4.1.6 (MSAV conversion)
- bcftools 1.20 (VCF processing)
- All dependencies (bc, bgzip, tabix)

### Local Setup
```bash
# Clone repository
git clone https://github.com/mamanambiya/generate-ref-panel-workshop-2025.git
cd generate-ref-panel-workshop-2025

# Pull container
docker pull mamana/imputation:minimac4-4.1.6

# Download Cromwell (if needed)
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
```

## Usage

### Input Configuration

Create an input JSON file based on your requirements:

**For Production (Cloud)**
```json
{
  "FederatedImputationPipeline.input_vcf": "gs://your-bucket/h3africa_chr22.vcf.gz",
  "FederatedImputationPipeline.chromosome": "22",
  "FederatedImputationPipeline.start_position": 16000000,
  "FederatedImputationPipeline.end_position": 51304566,
  "FederatedImputationPipeline.output_prefix": "h3africa_chr22_federated"
}
```

**For Local Testing**
```json
{
  "FederatedImputationPipeline.input_vcf": "/data/test_chr22.vcf.gz",
  "FederatedImputationPipeline.chromosome": "22",
  "FederatedImputationPipeline.start_position": 16000000,
  "FederatedImputationPipeline.end_position": 16200000,
  "FederatedImputationPipeline.output_prefix": "test_panel"
}
```

### Execution Commands

**Local Execution with Cromwell**
```bash
java -jar cromwell-85.jar run \
  workflows/federated_imputation_pipeline.wdl \
  -i inputs/production.json
```

**DNAstack Workbench**
1. Upload `workflows/federated_imputation_pipeline.wdl`
2. Upload input JSON configuration
3. Execute through web interface

**Terra/AnVIL**
1. Import workflow to workspace
2. Configure inputs via UI
3. Run workflow

### Parameter Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_vcf` | File | Required | Input VCF file (bgzip compressed) |
| `chromosome` | String | Required | Target chromosome (e.g., "22") |
| `start_position` | Int | Required | Start genomic position |
| `end_position` | Int | Required | End genomic position |
| `maf_threshold` | Float | 0.01 | Minimum allele frequency |
| `call_rate_threshold` | Float | 0.95 | Minimum call rate |
| `remove_indels` | Boolean | true | Keep SNPs only |
| `output_prefix` | String | "federated_panel" | Output file prefix |
| `compression_level` | Int | 5 | MSAV compression (1-9) |

## Testing

### Automated CI/CD Testing
The pipeline includes comprehensive GitHub Actions workflows for automated testing:

#### Unit Tests - Individual Component Testing
- **Triggers**: Every push/PR to `main`/`develop` branches
- **Tests**: 4 individual component unit tests
  - `ExtractRegion` task validation
  - `QualityControl` task validation  
  - `MinimacConversion` task validation
  - Complete pipeline integration testing
- **Matrix Strategy**: Parallel execution for faster feedback
- **Artifacts**: Test reports, validation results, logs
- **Status**: Automated via GitHub Actions on push/PR

#### Integration Tests - End-to-End Pipeline Testing  
- **Triggers**: Push to `main`, daily at 2 AM UTC, manual dispatch
- **Tests**: Complete workflow execution with validation
  - Main pipeline execution
  - Enhanced test pipeline with validation tasks
- **Outputs**: MSAV files, comprehensive reports, performance metrics
- **Status**: Automated via GitHub Actions daily and on main branch changes

#### Local Unit Testing
```bash
# Run comprehensive unit test suite locally
cd tests/
chmod +x run_unit_tests.sh
./run_unit_tests.sh

# Run individual unit tests
java -jar cromwell.jar run tests/unit/test_extract_region.wdl -i tests/inputs/unit_test_config.json
```

#### Manual Integration Testing
```bash
# Run complete pipeline test
java -jar cromwell.jar run tests/test_pipeline.wdl -i inputs/test_local.json

# Run main pipeline
java -jar cromwell.jar run workflows/federated_imputation_pipeline.wdl -i inputs/test_local.json
```

#### Test Results & Reporting
- **Automatic PR Comments**: Test results posted as PR comments
- **Detailed Reports**: Validation reports with pass/fail status
- **Artifacts Storage**: 30-90 day retention for test outputs
- **Success Metrics**: Real-time success rates and statistics

### Manual Validation
```bash
# Validate MSAV output
file output/federated_panel.msav
# Expected: Zstandard compressed data

# Check variant counts
bcftools view -H output/quality_controlled.vcf.gz | wc -l
```

### Test Data
- `test_chr22_region_bgzip.vcf.gz` - 100 variants from chromosome 22
- Expected output: ~18 variants after QC, ~5KB MSAV file

### Test Coverage
- **Unit Tests**: Individual task validation (4 tests)
- **Integration Tests**: End-to-end pipeline validation (2 tests)  
- **Performance Tests**: Resource usage and timing validation
- **Cross-platform**: Docker platform compatibility testing
- **Quality Gates**: Automated pass/fail determination

## Deployment

### GitHub Actions CI/CD Pipeline
The repository includes automated CI/CD workflows for quality assurance:

#### Continuous Integration
```yaml
# .github/workflows/unit-tests.yml
on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main, develop ]
```

- **Automated Testing**: Every commit triggers comprehensive unit tests
- **Quality Gates**: PRs blocked if tests fail
- **Parallel Execution**: Matrix strategy for faster feedback
- **Cross-platform**: Tests run on `ubuntu-latest` with platform compatibility

#### Continuous Deployment Readiness
```yaml
# .github/workflows/integration-tests.yml  
on:
  push:
    branches: [ main ]
  schedule:
    - cron: '0 2 * * *'  # Daily monitoring
```

- **Production Validation**: Main branch changes trigger integration tests
- **Daily Monitoring**: Automated health checks at 2 AM UTC
- **Release Readiness**: Validates complete pipeline before deployment

#### Manual Triggers
```yaml
workflow_dispatch:
  inputs:
    test_scope:
      type: choice
      options: [all, extract_region, quality_control, minimac_conversion]
```

- **On-demand Testing**: Run specific test suites manually
- **Debug Mode**: Individual component testing for troubleshooting
- **Performance Testing**: Optional stress testing configurations

### DNAstack Workbench
1. **Upload Files**:
   - `workflows/federated_imputation_pipeline.wdl`
   - `inputs/production.json` (configured for your data)

2. **Configure Container**: Ensure `mamana/imputation:minimac4-4.1.6` is available

3. **Execute**: Run workflow through web interface

### Cloud Platforms

**Google Cloud (Terra)**
- Import workflow to Terra workspace
- Configure with Google Cloud Storage paths
- Execute with Cromwell backend

**AWS (Amazon Genomics)**
- Deploy using AWS Batch backend
- Configure S3 storage paths
- Scale with EC2 instances

**Azure (Microsoft Genomics)**
- Use Azure Batch for execution
- Configure Azure Blob Storage
- Integrate with Azure AD

### HPC Clusters
```bash
# Singularity setup
singularity pull docker://mamana/imputation:minimac4-4.1.6

# SLURM submission script
sbatch --job-name=federated-imputation run_pipeline.sh
```

## Production Use Cases

### African Genomics Research
- **H3Africa Consortium**: Process reference panels while maintaining data sovereignty
- **AWI-Gen Study**: Enable cross-site imputation without data transfer
- **Population Studies**: Create federated imputation networks across institutions

### International Collaboration
- **Multi-site Studies**: Collaborative research without data sharing
- **Reference Panel Development**: Contribute to global resources while preserving privacy
- **Regulatory Compliance**: Meet data protection requirements while enabling research

## Performance

### Benchmarks
- **Small regions** (100K variants): 5-10 minutes
- **Chromosome arm** (500K variants): 15-30 minutes  
- **Whole chromosome** (2M+ variants): 1-2 hours

### Resource Requirements
- **CPU**: 2-4 cores recommended
- **Memory**: 4-8 GB depending on dataset size
- **Storage**: 2x input file size for temporary files

## Contributing

### Development Guidelines
1. Follow WDL best practices for task organization
2. Include comprehensive parameter metadata
3. Add validation and error handling
4. Update documentation for new features

### Testing Requirements
- All tasks must include unit tests
- Integration tests for full pipeline
- Validation with multiple VCF formats

### Pull Request Process
1. Create feature branch
2. Add tests for new functionality
3. Update documentation
4. Submit PR with clear description

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Support

### Documentation
- [WDL Specification](https://github.com/openwdl/wdl)
- [GA4GH Standards](https://www.ga4gh.org/)
- [Minimac4 Documentation](https://genome.sph.umich.edu/wiki/Minimac4)

### Contact
- **GA4GH Hackathon 2025 Team**
- **African Genomics Initiative**
- **Email**: [team@afrigenomics.org](mailto:team@afrigenomics.org)

---

## Ready for African Genomics Research!

This pipeline enables federated genotype imputation across African institutions while preserving data sovereignty and promoting collaborative genomics research.

**Container**: `mamana/imputation:minimac4-4.1.6`  
**Workflow**: Production-ready WDL pipeline  
**Output**: MSAV files for federated networks  
**Status**: **READY FOR DEPLOYMENT**
