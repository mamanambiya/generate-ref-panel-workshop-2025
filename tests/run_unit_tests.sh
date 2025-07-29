#!/bin/bash

## Federated Imputation Pipeline - Unit Test Runner
## GA4GH Hackathon 2025 - African Genomics Team

set -euo pipefail

echo "🧬 Federated Imputation Pipeline - Unit Test Suite"
echo "GA4GH Hackathon 2025 - African Genomics Team"
echo "=================================================="
echo ""

# Configuration
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$TEST_DIR")"
CROMWELL_JAR="$ROOT_DIR/cromwell.jar"
RESULTS_DIR="$ROOT_DIR/test_results/unit_tests"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test results tracking
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Create results directory
mkdir -p "$RESULTS_DIR"
echo "Test results will be saved to: $RESULTS_DIR"
echo ""

# Check prerequisites
check_prerequisites() {
    echo "🔍 Checking prerequisites..."
    
    # Check if Cromwell exists
    if [[ ! -f "$CROMWELL_JAR" ]]; then
        echo -e "${RED}❌ Cromwell not found at $CROMWELL_JAR${NC}"
        echo "Run the main pipeline test first to download Cromwell"
        exit 1
    fi
    
    # Check if Docker is running
    if ! docker info >/dev/null 2>&1; then
        echo -e "${RED}❌ Docker is not running${NC}"
        exit 1
    fi
    
    # Check if test data exists
    if [[ ! -f "$ROOT_DIR/tests/test_chr22_region_bgzip.vcf.gz" ]]; then
        echo -e "${RED}❌ Test data not found${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}✅ All prerequisites met${NC}"
    echo ""
}

# Run individual test
run_test() {
    local test_name="$1"
    local workflow_file="$2"
    local input_file="$3"
    local description="$4"
    
    echo -e "${BLUE}🧪 Running $test_name${NC}"
    echo "Description: $description"
    echo "Workflow: $workflow_file"
    echo "Inputs: $input_file"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    # Run the test
    if java -jar "$CROMWELL_JAR" run "$workflow_file" -i "$input_file" > "$RESULTS_DIR/${test_name}_output.log" 2>&1; then
        echo -e "${GREEN}✅ $test_name PASSED${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        
        # Copy validation report if it exists
        local validation_file=$(find cromwell-executions -name "*validation_report.txt" -type f | head -1)
        if [[ -n "$validation_file" ]]; then
            cp "$validation_file" "$RESULTS_DIR/${test_name}_validation_report.txt"
        fi
        
    else
        echo -e "${RED}❌ $test_name FAILED${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        
        # Show last few lines of error
        echo "Error details:"
        tail -20 "$RESULTS_DIR/${test_name}_output.log" | head -10
    fi
    
    echo ""
}

# Run all unit tests
run_unit_tests() {
    echo "🚀 Starting unit test execution..."
    echo ""
    
    # Test 1: ExtractRegion task
    run_test \
        "extract_region_test" \
        "$TEST_DIR/unit/test_extract_region.wdl" \
        "$TEST_DIR/inputs/unit_test_config.json" \
        "Test genomic region extraction functionality"
    
    # Test 2: QualityControl task
    run_test \
        "quality_control_test" \
        "$TEST_DIR/unit/test_quality_control.wdl" \
        "$TEST_DIR/inputs/qc_unit_test_config.json" \
        "Test VCF quality control and filtering"
    
    # Test 3: MinimacConversion task
    run_test \
        "minimac_conversion_test" \
        "$TEST_DIR/unit/test_minimac_conversion.wdl" \
        "$TEST_DIR/inputs/minimac_unit_test_config.json" \
        "Test MSAV format conversion using Minimac4"
    
    # Test 4: Complete pipeline with multiple scenarios
    run_test \
        "complete_pipeline_test" \
        "$TEST_DIR/unit/test_complete_pipeline.wdl" \
        "$TEST_DIR/inputs/complete_pipeline_test_config.json" \
        "Test complete pipeline with basic, strict, and relaxed filtering"
}

# Generate test report
generate_report() {
    local report_file="$RESULTS_DIR/UNIT_TEST_REPORT.md"
    
    echo "📝 Generating test report..."
    
    cat > "$report_file" << EOF
# Federated Imputation Pipeline - Unit Test Report

## 🧬 GA4GH Hackathon 2025 - African Genomics Team

**Test Date**: $(date)  
**Test Environment**: $(uname -s) $(uname -m)  
**Docker Version**: $(docker --version)

---

## 📊 Test Summary

| Metric | Value |
|--------|-------|
| **Total Tests** | $TOTAL_TESTS |
| **Passed Tests** | $PASSED_TESTS |
| **Failed Tests** | $FAILED_TESTS |
| **Success Rate** | $(echo "scale=1; $PASSED_TESTS * 100 / $TOTAL_TESTS" | bc)% |

---

## 🧪 Test Results

### Individual Component Tests

EOF

    # Add results for each test
    if [[ -f "$RESULTS_DIR/extract_region_test_validation_report.txt" ]]; then
        echo "#### ✅ ExtractRegion Task Test" >> "$report_file"
        echo '```' >> "$report_file"
        cat "$RESULTS_DIR/extract_region_test_validation_report.txt" >> "$report_file"
        echo '```' >> "$report_file"
        echo "" >> "$report_file"
    fi
    
    if [[ -f "$RESULTS_DIR/quality_control_test_validation_report.txt" ]]; then
        echo "#### ✅ QualityControl Task Test" >> "$report_file"
        echo '```' >> "$report_file"
        cat "$RESULTS_DIR/quality_control_test_validation_report.txt" >> "$report_file"
        echo '```' >> "$report_file"
        echo "" >> "$report_file"
    fi
    
    if [[ -f "$RESULTS_DIR/minimac_conversion_test_validation_report.txt" ]]; then
        echo "#### ✅ MinimacConversion Task Test" >> "$report_file"
        echo '```' >> "$report_file"
        cat "$RESULTS_DIR/minimac_conversion_test_validation_report.txt" >> "$report_file"
        echo '```' >> "$report_file"
        echo "" >> "$report_file"
    fi
    
    if [[ -f "$RESULTS_DIR/complete_pipeline_test_validation_report.txt" ]]; then
        echo "#### ✅ Complete Pipeline Test" >> "$report_file"
        echo '```' >> "$report_file"
        cat "$RESULTS_DIR/complete_pipeline_test_validation_report.txt" >> "$report_file"
        echo '```' >> "$report_file"
        echo "" >> "$report_file"
    fi
    
    cat >> "$report_file" << EOF

---

## 📁 Test Artifacts

- **Test Logs**: \`test_results/unit_tests/*_output.log\`
- **Validation Reports**: \`test_results/unit_tests/*_validation_report.txt\`
- **Generated Files**: Check cromwell-executions directories

---

## 🎯 Test Coverage

The unit test suite validates:

1. **Individual Task Functionality**
   - ✅ Region extraction accuracy
   - ✅ Quality control filtering
   - ✅ MSAV format conversion
   - ✅ File format validation
   - ✅ Parameter handling

2. **Pipeline Integration**
   - ✅ Task interdependencies
   - ✅ Data flow between stages
   - ✅ Multiple filtering scenarios
   - ✅ Error handling

3. **Output Validation**
   - ✅ File existence and format
   - ✅ Content accuracy
   - ✅ Statistical summaries
   - ✅ Cross-scenario consistency

---

## 🚀 Production Readiness

EOF

    if [[ $FAILED_TESTS -eq 0 ]]; then
        cat >> "$report_file" << EOF
**Status**: ✅ **ALL TESTS PASSED**

The federated imputation pipeline has successfully passed all unit tests and is ready for production deployment.

### Key Validations Completed:
- ✅ All individual tasks working correctly
- ✅ Complete pipeline integration verified
- ✅ Multiple filtering scenarios tested
- ✅ MSAV output format validated
- ✅ Error handling verified

**Recommendation**: ✅ **APPROVED FOR GA4GH HACKATHON 2025**
EOF
    else
        cat >> "$report_file" << EOF
**Status**: ❌ **SOME TESTS FAILED**

$FAILED_TESTS out of $TOTAL_TESTS tests failed. Review the individual test reports above for details.

**Recommendation**: 🔄 **REQUIRES FIXES BEFORE PRODUCTION**
EOF
    fi
    
    echo "Report generated: $report_file"
}

# Cleanup function
cleanup() {
    echo ""
    echo "🧹 Cleaning up..."
    # Remove large cromwell execution directories if tests passed
    if [[ $FAILED_TESTS -eq 0 ]]; then
        find cromwell-executions -name "*Test*" -type d -exec rm -rf {} + 2>/dev/null || true
        echo "Cleaned up successful test executions"
    fi
}

# Main execution
main() {
    # Change to root directory
    cd "$ROOT_DIR"
    
    # Run tests
    check_prerequisites
    run_unit_tests
    generate_report
    cleanup
    
    # Final summary
    echo "=================================================="
    echo "🎯 Unit Test Suite Complete"
    echo ""
    echo "Results Summary:"
    echo -e "  Total Tests: $TOTAL_TESTS"
    echo -e "  Passed: ${GREEN}$PASSED_TESTS${NC}"
    echo -e "  Failed: ${RED}$FAILED_TESTS${NC}"
    echo -e "  Success Rate: $(echo "scale=1; $PASSED_TESTS * 100 / $TOTAL_TESTS" | bc)%"
    echo ""
    echo "Detailed results available in: $RESULTS_DIR"
    
    if [[ $FAILED_TESTS -eq 0 ]]; then
        echo -e "${GREEN}🎉 ALL TESTS PASSED - Pipeline ready for production!${NC}"
        exit 0
    else
        echo -e "${RED}❌ Some tests failed - Review reports for details${NC}"
        exit 1
    fi
}

# Handle script interruption
trap cleanup EXIT

# Run main function
main "$@" 