# Federated Imputation Pipeline - Testing Framework

## 🧬 GA4GH Hackathon 2025 - African Genomics Team

This directory contains a comprehensive unit testing framework for the Federated Genotype Imputation Pipeline. The tests validate individual components and the complete pipeline across multiple scenarios.

---

## 📁 **Directory Structure**

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

## 🧪 **Test Categories**

### **1. Individual Component Tests**

#### **ExtractRegion Task Test**
- **File**: `unit/test_extract_region.wdl`
- **Purpose**: Validates genomic region extraction functionality
- **Validations**:
  - ✅ File existence (VCF, summary, stats)
  - ✅ Variant count accuracy
  - ✅ Chromosome consistency  
  - ✅ Position range validation
  - ✅ Summary content verification

#### **QualityControl Task Test**
- **File**: `unit/test_quality_control.wdl`
- **Purpose**: Tests VCF filtering and quality control
- **Validations**:
  - ✅ MAF threshold filtering
  - ✅ Call rate filtering
  - ✅ Indel removal (when requested)
  - ✅ Variant count consistency
  - ✅ Summary report accuracy

#### **MinimacConversion Task Test**
- **File**: `unit/test_minimac_conversion.wdl`
- **Purpose**: Validates MSAV format conversion
- **Validations**:
  - ✅ MSAV file creation
  - ✅ Zstandard compression format
  - ✅ File size reasonableness
  - ✅ Filename accuracy
  - ✅ Conversion summary content

### **2. Integration Tests**

#### **Complete Pipeline Test**
- **File**: `unit/test_complete_pipeline.wdl`
- **Purpose**: Tests complete pipeline with multiple scenarios
- **Test Scenarios**:
  - **Basic**: Standard filtering (MAF≥0.01, CR≥0.95)
  - **Strict**: Stringent filtering (MAF≥0.10, CR≥0.99)
  - **Relaxed**: Lenient filtering (MAF≥0.001, CR≥0.80)
- **Cross-Scenario Validations**:
  - ✅ Format consistency across outputs
  - ✅ Filtering effect validation
  - ✅ Performance comparison
  - ✅ File size analysis

---

## 🚀 **Running Tests**

### **Quick Start**

```bash
# Make test runner executable
chmod +x tests/run_unit_tests.sh

# Run all unit tests
./tests/run_unit_tests.sh
```

### **Individual Test Execution**

```bash
# Test ExtractRegion task only
java -jar cromwell.jar run tests/unit/test_extract_region.wdl \
    -i tests/inputs/unit_test_config.json

# Test QualityControl task only  
java -jar cromwell.jar run tests/unit/test_quality_control.wdl \
    -i tests/inputs/qc_unit_test_config.json

# Test MinimacConversion task only
java -jar cromwell.jar run tests/unit/test_minimac_conversion.wdl \
    -i tests/inputs/minimac_unit_test_config.json

# Test complete pipeline with multiple scenarios
java -jar cromwell.jar run tests/unit/test_complete_pipeline.wdl \
    -i tests/inputs/complete_pipeline_test_config.json
```

---

## 📊 **Test Output & Reports**

### **Automated Test Results**
- **Location**: `test_results/unit_tests/`
- **Main Report**: `UNIT_TEST_REPORT.md`
- **Individual Logs**: `*_output.log`
- **Validation Reports**: `*_validation_report.txt`

### **Test Result Interpretation**

#### **Success Indicators**
- ✅ All validation checks pass
- ✅ Files created with expected formats
- ✅ Variant counts within expected ranges
- ✅ No error messages in logs

#### **Failure Indicators**
- ❌ Missing output files
- ❌ Incorrect file formats
- ❌ Variant count inconsistencies  
- ❌ Error messages in logs

---

## 🎯 **Test Coverage**

### **Functional Coverage**
- [x] **Region Extraction**: All genomic coordinates and chromosome handling
- [x] **Quality Filtering**: MAF, call rate, and variant type filters
- [x] **Format Conversion**: VCF to MSAV transformation
- [x] **Pipeline Integration**: Complete workflow execution
- [x] **Parameter Handling**: Multiple configuration scenarios

### **Data Coverage**
- [x] **Test Dataset**: 100 variants from chr22:16000000-16990000
- [x] **Multiple Scenarios**: Basic, strict, and relaxed filtering
- [x] **Edge Cases**: Empty regions, strict filters, format variations

### **Platform Coverage**
- [x] **Docker Integration**: Containerized execution
- [x] **Cross-Platform**: Apple Silicon (ARM64) with linux/amd64 emulation
- [x] **WDL Compliance**: Standards-compliant workflow execution

---

## 🔧 **Customizing Tests**

### **Adding New Test Scenarios**

1. **Create New WDL Test**:
   ```wdl
   version 1.0
   
   import "../../tasks/your_task.wdl" as YourTask
   
   workflow TestYourTask {
       # Define test inputs and validation logic
   }
   ```

2. **Add Input Configuration**:
   ```json
   {
     "TestYourTask.input_parameter": "test_value"
   }
   ```

3. **Update Test Runner**:
   Add new test case to `run_unit_tests.sh`

### **Modifying Test Parameters**

Edit input JSON files in `tests/inputs/` to change:
- Test regions and coordinates
- Filtering thresholds
- Output prefixes
- Compression levels

### **Adding Custom Validation**

Enhance validation tasks in test WDL files to check:
- Additional file properties
- Content-specific validations
- Performance metrics
- Error condition handling

---

## 🐛 **Troubleshooting**

### **Common Issues**

#### **Docker Platform Errors**
```bash
# Error: platform mismatch
# Solution: Tests include --platform linux/amd64 flag
```

#### **Missing Test Data**
```bash
# Error: test_chr22_region_bgzip.vcf.gz not found
# Solution: Ensure test data is in tests/ directory
```

#### **Cromwell Download Issues**
```bash
# Error: cromwell.jar not found
# Solution: Run main pipeline test first or download manually
curl -L -o cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
```

#### **Permission Issues**
```bash
# Error: Permission denied
# Solution: Make test runner executable
chmod +x tests/run_unit_tests.sh
```

### **Log Analysis**

Check these files for debugging:
- `test_results/unit_tests/*_output.log` - Full execution logs
- `test_results/unit_tests/*_validation_report.txt` - Detailed validation results
- `cromwell-executions/*/call-*/execution/stderr` - Task-specific errors

---

## 📈 **Performance Expectations**

### **Test Execution Times**
- **ExtractRegion**: ~10 seconds
- **QualityControl**: ~15 seconds  
- **MinimacConversion**: ~10 seconds
- **Complete Pipeline**: ~2-3 minutes (3 scenarios)
- **Total Suite**: ~4-5 minutes

### **Resource Usage**
- **Memory**: 2-4 GB per task
- **CPU**: 1-2 cores per task
- **Disk**: ~100 MB temporary files
- **Network**: Container download (first run only)

---

## ✅ **Quality Assurance**

### **Test Reliability**
- **Deterministic**: Same inputs always produce same outputs
- **Isolated**: Each test runs independently
- **Comprehensive**: Covers all major functionality
- **Fast**: Complete suite runs in under 5 minutes

### **Production Readiness Criteria**
- ✅ All unit tests pass
- ✅ All integration tests pass  
- ✅ Performance within acceptable ranges
- ✅ Error handling validated
- ✅ Documentation complete

---

## 🤝 **Contributing**

### **Adding New Tests**
1. Create test WDL workflow
2. Add input configuration
3. Update test runner script
4. Document new test case
5. Verify all tests still pass

### **Reporting Issues**
Include in bug reports:
- Test execution logs
- Validation reports
- System information
- Steps to reproduce

---

## 📞 **Support**

- **Team**: GA4GH Hackathon 2025 - African Genomics Initiative
- **Documentation**: Complete README files in each directory  
- **Container**: `mamana/imputation:minimac4-4.1.6`
- **WDL Version**: 1.0

**Status**: ✅ **COMPREHENSIVE TESTING FRAMEWORK READY** 