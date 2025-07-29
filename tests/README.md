# Federated Imputation Pipeline - Testing Framework

## GA4GH Hackathon 2025 - African Genomics Team

This directory contains a comprehensive unit testing framework for the Federated Genotype Imputation Pipeline. The tests validate individual components and the complete pipeline across multiple scenarios.

---

## Directory Structure

```
tests/
├── README.md                           # This documentation
├── run_unit_tests.sh                   # Main test runner script
├── test_pipeline.wdl                   # Integration test workflow
├── test_chr22_region_bgzip.vcf.gz      # Test data (100 variants)
├── test_chr22_region_bgzip.vcf.gz.csi  # Index file
│
├── unit/                               # Individual component tests
│   ├── test_extract_region.wdl        # ExtractRegion task test
│   ├── test_quality_control.wdl       # QualityControl task test  
│   ├── test_minimac_conversion.wdl     # MinimacConversion task test
│   └── test_complete_pipeline.wdl     # Multi-scenario pipeline test
│
└── inputs/                            # Test input configurations
    ├── unit_test_config.json          # ExtractRegion test config
    ├── qc_unit_test_config.json       # QualityControl test config
    ├── minimac_unit_test_config.json  # MinimacConversion test config
    └── complete_pipeline_test_config.json # Complete pipeline config
```

---

## Test Categories

### 1. Individual Component Tests

**ExtractRegion Task Test** (`test_extract_region.wdl`)
- Tests genomic region extraction functionality
- Validates VCF file parsing and bcftools integration
- Checks output file generation and statistics
- Verifies genomic coordinate accuracy

**QualityControl Task Test** (`test_quality_control.wdl`)
- Tests VCF filtering and quality control
- Validates MAF and call rate filtering
- Checks indel removal functionality
- Verifies QC statistics and reporting

**MinimacConversion Task Test** (`test_minimac_conversion.wdl`)
- Tests MSAV format conversion using Minimac4
- Validates file format and compression
- Checks conversion statistics and metadata
- Verifies output file integrity

### 2. Integration Tests

**Complete Pipeline Test** (`test_complete_pipeline.wdl`)
- Tests the full federated imputation pipeline
- Runs multiple filtering scenarios:
  - Basic configuration (standard filtering)
  - Strict filtering (high thresholds)
  - Relaxed filtering (low thresholds)
- Validates cross-scenario consistency
- Comprehensive end-to-end testing

**Integration Test** (`test_pipeline.wdl`)
- Tests main pipeline workflow
- Includes additional validation tasks
- Generates comprehensive test reports
- Validates production readiness

---

## Quick Start

### Run All Unit Tests
```bash
# Navigate to tests directory
cd tests/

# Make test runner executable
chmod +x run_unit_tests.sh

# Run complete test suite
./run_unit_tests.sh
```

### Run Individual Tests
```bash
# Extract Region test
java -jar cromwell.jar run unit/test_extract_region.wdl -i inputs/unit_test_config.json

# Quality Control test
java -jar cromwell.jar run unit/test_quality_control.wdl -i inputs/qc_unit_test_config.json

# Minimac Conversion test
java -jar cromwell.jar run unit/test_minimac_conversion.wdl -i inputs/minimac_unit_test_config.json

# Complete Pipeline test
java -jar cromwell.jar run unit/test_complete_pipeline.wdl -i inputs/complete_pipeline_test_config.json
```

---

## Test Output

### Test Runner Output
```
Federated Imputation Pipeline - Unit Test Suite
GA4GH Hackathon 2025 - African Genomics Team
==================================================

Checking prerequisites...
All prerequisites met

Running extract_region_test...
Description: Test genomic region extraction functionality
extract_region_test PASSED

Running quality_control_test...
Description: Test VCF quality control and filtering
quality_control_test PASSED

Running minimac_conversion_test...
Description: Test MSAV format conversion using Minimac4
minimac_conversion_test PASSED

Running complete_pipeline_test...
Description: Test complete pipeline with basic, strict, and relaxed filtering
complete_pipeline_test PASSED

==================================================
Unit Test Suite Complete

Results Summary:
  Total Tests: 4
  Passed: 4
  Failed: 0
  Success Rate: 100.0%

All unit tests passed - Pipeline ready for production!
```

### Individual Test Reports

Each test generates detailed validation reports:

**ExtractRegion Validation Report**
```
ExtractRegion Unit Test Validation
=================================
Test Date: [timestamp]

File Existence Tests:
✅ Extracted VCF file exists
✅ Summary file exists
✅ Stats file exists

Content Validation Tests:
Variant count: 21
✅ VCF contains variants
✅ Chromosome matches expected (22)
✅ Positions within expected range (16000000-16200000)
   Actual range: 16000000-16200000
✅ Summary contains variant count
✅ Summary contains correct region

Overall Test Result:
TEST PASSED: ExtractRegion working correctly
```

---

## Test Coverage

### Unit Test Coverage
- **ExtractRegion**: File parsing, region extraction, statistics generation
- **QualityControl**: MAF filtering, call rate filtering, indel removal, reporting
- **MinimacConversion**: MSAV conversion, file validation, metadata generation
- **Pipeline Integration**: Multi-task workflow, data flow, error propagation

### Validation Scope
- Input file validation
- Output file generation
- Content verification
- Format compliance
- Error handling
- Performance metrics

### Cross-Platform Testing
- Docker platform compatibility (`linux/amd64`)
- Container execution validation
- Resource allocation testing
- Dependency verification

---

## Configuration

### Test Data
- **File**: `test_chr22_region_bgzip.vcf.gz`
- **Size**: 100 variants from chromosome 22
- **Range**: 16000000-16990000 (chr22)
- **Format**: bgzip compressed VCF with CSI index

### Input Configurations
All test configurations use absolute paths and realistic parameters:

```json
{
  "test_vcf": "/absolute/path/to/test_chr22_region_bgzip.vcf.gz",
  "test_chromosome": "22",
  "test_start": 16000000,
  "test_end": 16200000,
  "test_prefix": "unit_test_output"
}
```

### Expected Results
- **ExtractRegion**: ~21 variants in extracted region
- **QualityControl**: 10-18 variants after filtering (depends on parameters)
- **MinimacConversion**: Valid MSAV file (Zstandard compressed)
- **Pipeline**: Complete MSAV panel ready for federated imputation

---

## Troubleshooting

### Common Issues

**Cromwell Not Found**
```bash
# Download Cromwell if missing
curl -L -o cromwell.jar \
  "https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar"
```

**Docker Platform Issues**
```bash
# Pull with platform specification
docker pull --platform linux/amd64 mamana/imputation:minimac4-4.1.6
```

**Permission Issues**
```bash
# Make test runner executable
chmod +x run_unit_tests.sh
```

**Path Issues**
- Ensure test data files exist: `tests/test_chr22_region_bgzip.vcf.gz`
- Update paths in input JSON files for your environment
- Use absolute paths for Cromwell execution

### Debug Mode
```bash
# Run with verbose output
java -jar cromwell.jar run unit/test_extract_region.wdl \
  -i inputs/unit_test_config.json \
  --options verbose_options.json
```

---

## Performance

### Expected Execution Times
- **ExtractRegion**: 30-60 seconds
- **QualityControl**: 45-90 seconds  
- **MinimacConversion**: 60-120 seconds
- **Complete Pipeline**: 3-5 minutes total

### Resource Usage
- **Memory**: 2-4 GB per test
- **CPU**: 1-2 cores
- **Storage**: 100-500 MB temporary files

### Optimization Tips
- Use local test data to avoid network delays
- Ensure sufficient Docker resources
- Run tests on SSD storage for better I/O performance

---

## Quality Assurance

### Validation Criteria
- All output files must be generated
- File formats must be valid (VCF, MSAV)
- Content must match expected patterns
- No critical errors in execution logs
- Resource usage within expected limits

### Acceptance Thresholds
- Test success rate: 100%
- Execution time: < 2x expected duration
- Memory usage: < 150% of allocated resources
- Output file sizes: Within 10% of expected

### Continuous Integration
Tests are automatically run via GitHub Actions:
- On every push to main/develop branches
- On all pull requests
- Daily monitoring at 2 AM UTC
- Manual triggers for debugging

---

## Contributing

### Adding New Tests
1. Create new WDL workflow in `unit/` directory
2. Add corresponding input JSON in `inputs/` directory
3. Update `run_unit_tests.sh` to include new test
4. Update this documentation

### Test Standards
- Use descriptive test names
- Include comprehensive validation
- Generate detailed reports
- Handle error cases gracefully
- Document expected behavior

### Best Practices
- Test edge cases and boundary conditions
- Validate both success and failure scenarios
- Use realistic test data
- Minimize test execution time
- Provide clear failure messages

---

## Support

### Documentation
- **Main README**: `../README.md`
- **Contributing Guide**: `../.github/CONTRIBUTING.md`
- **WDL Documentation**: [OpenWDL Specification](https://github.com/openwdl/wdl)

### Getting Help
- Review test output logs in `test_results/` directory
- Check Cromwell execution logs for detailed errors
- Validate input file paths and permissions
- Ensure Docker containers are accessible

### Contact
- **Team**: GA4GH Hackathon 2025 - African Genomics Team
- **Repository**: [generate-ref-panel-workshop-2025](https://github.com/mamanambiya/generate-ref-panel-workshop-2025)
- **Issues**: Use GitHub Issues for bug reports and feature requests

---

**The testing framework ensures production readiness for federated genotype imputation across African genomics institutions while maintaining data sovereignty.** 