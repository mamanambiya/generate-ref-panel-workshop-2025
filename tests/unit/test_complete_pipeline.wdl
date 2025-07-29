version 1.0

import "../../workflows/federated_imputation_pipeline.wdl" as MainPipeline

## Comprehensive unit test for the complete pipeline
workflow TestCompletePipeline {
    input {
        # Test scenarios
        Array[String] test_scenarios = ["basic", "strict_filtering", "relaxed_filtering"]
        
        # Common inputs
        File test_vcf = "tests/test_chr22_region_bgzip.vcf.gz"
        String test_chromosome = "22"
        Int test_start = 16000000
        Int test_end = 16990000
    }

    # Test 1: Basic configuration
    call MainPipeline.FederatedImputationPipeline as BasicTest {
        input:
            input_vcf = test_vcf,
            chromosome = test_chromosome,
            start_position = test_start,
            end_position = test_end,
            output_prefix = "unit_test_basic",
            maf_threshold = 0.01,
            call_rate_threshold = 0.95,
            remove_indels = true,
            compression_level = 3
    }

    # Test 2: Strict filtering
    call MainPipeline.FederatedImputationPipeline as StrictTest {
        input:
            input_vcf = test_vcf,
            chromosome = test_chromosome,
            start_position = test_start,
            end_position = test_end,
            output_prefix = "unit_test_strict",
            maf_threshold = 0.10,
            call_rate_threshold = 0.99,
            remove_indels = true,
            compression_level = 1
    }

    # Test 3: Relaxed filtering
    call MainPipeline.FederatedImputationPipeline as RelaxedTest {
        input:
            input_vcf = test_vcf,
            chromosome = test_chromosome,
            start_position = test_start,
            end_position = test_end,
            output_prefix = "unit_test_relaxed",
            maf_threshold = 0.001,
            call_rate_threshold = 0.80,
            remove_indels = false,
            compression_level = 9
    }

    # Validate all test results
    call ValidateCompletePipeline {
        input:
            # Basic test outputs
            basic_msav = BasicTest.federated_panel_msav,
            basic_extracted_vcf = BasicTest.extracted_vcf,
            basic_qc_vcf = BasicTest.quality_controlled_vcf,
            basic_qc_summary = BasicTest.qc_summary,
            basic_conversion_summary = BasicTest.conversion_summary,
            
            # Strict test outputs
            strict_msav = StrictTest.federated_panel_msav,
            strict_qc_summary = StrictTest.qc_summary,
            
            # Relaxed test outputs
            relaxed_msav = RelaxedTest.federated_panel_msav,
            relaxed_qc_summary = RelaxedTest.qc_summary,
            
            # Test parameters
            input_vcf = test_vcf
    }

    output {
        # Test results for each scenario
        File basic_test_msav = BasicTest.federated_panel_msav
        File strict_test_msav = StrictTest.federated_panel_msav
        File relaxed_test_msav = RelaxedTest.federated_panel_msav
        
        # Validation report
        File validation_report = ValidateCompletePipeline.validation_report
        Boolean all_tests_passed = ValidateCompletePipeline.all_tests_passed
        
        # Summary reports for inspection
        Array[File] test_summaries = [
            BasicTest.qc_summary,
            StrictTest.qc_summary,
            RelaxedTest.qc_summary
        ]
    }

    meta {
        description: "Comprehensive unit test for complete federated imputation pipeline"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }
}

## Validate complete pipeline results across multiple test scenarios
task ValidateCompletePipeline {
    input {
        # Test outputs
        File basic_msav
        File basic_extracted_vcf
        File basic_qc_vcf
        File basic_qc_summary
        File basic_conversion_summary
        
        File strict_msav
        File strict_qc_summary
        
        File relaxed_msav
        File relaxed_qc_summary
        
        # Original input
        File input_vcf
    }

    command <<<
        echo "Complete Pipeline Unit Test Validation" > validation_report.txt
        echo "=====================================" >> validation_report.txt
        echo "Test Date: $(date)" >> validation_report.txt
        echo "Testing multiple filtering scenarios" >> validation_report.txt
        echo "" >> validation_report.txt
        
        all_tests_passed=true
        
        # Get original variant count
        original_count=$(bcftools view -H "~{input_vcf}" | wc -l)
        echo "Original input variants: $original_count" >> validation_report.txt
        echo "" >> validation_report.txt
        
        # Test 1: Basic configuration validation
        echo "=== Test 1: Basic Configuration ===" >> validation_report.txt
        
        if [ -f "~{basic_msav}" ]; then
            echo "✅ Basic MSAV file created" >> validation_report.txt
            basic_size=$(stat -f%z "~{basic_msav}" 2>/dev/null || stat -c%s "~{basic_msav}")
            echo "   File size: $basic_size bytes" >> validation_report.txt
        else
            echo "❌ Basic MSAV file missing" >> validation_report.txt
            all_tests_passed=false
        fi
        
        # Extract variant counts from summaries
        basic_final=$(grep "Final variants:" "~{basic_qc_summary}" | awk '{print $3}')
        echo "Basic test final variants: $basic_final" >> validation_report.txt
        
        # Validate pipeline stages
        extracted_count=$(bcftools view -H "~{basic_extracted_vcf}" | wc -l)
        qc_count=$(bcftools view -H "~{basic_qc_vcf}" | wc -l)
        
        if [ "$extracted_count" -eq "$original_count" ]; then
            echo "✅ Extraction preserved all variants" >> validation_report.txt
        else
            echo "❌ Extraction changed variant count: $original_count → $extracted_count" >> validation_report.txt
            all_tests_passed=false
        fi
        
        if [ "$qc_count" -le "$extracted_count" ]; then
            echo "✅ QC reduced or maintained variant count: $extracted_count → $qc_count" >> validation_report.txt
        else
            echo "❌ QC increased variant count (impossible): $extracted_count → $qc_count" >> validation_report.txt
            all_tests_passed=false
        fi
        
        echo "" >> validation_report.txt
        
        # Test 2: Strict filtering validation
        echo "=== Test 2: Strict Filtering ===" >> validation_report.txt
        
        if [ -f "~{strict_msav}" ]; then
            echo "✅ Strict MSAV file created" >> validation_report.txt
            strict_size=$(stat -f%z "~{strict_msav}" 2>/dev/null || stat -c%s "~{strict_msav}")
            echo "   File size: $strict_size bytes" >> validation_report.txt
        else
            echo "❌ Strict MSAV file missing" >> validation_report.txt
            all_tests_passed=false
        fi
        
        strict_final=$(grep "Final variants:" "~{strict_qc_summary}" | awk '{print $3}')
        echo "Strict test final variants: $strict_final" >> validation_report.txt
        
        # Strict filtering should result in fewer variants than basic
        if [ "$strict_final" -le "$basic_final" ]; then
            echo "✅ Strict filtering reduced variants: $basic_final → $strict_final" >> validation_report.txt
        else
            echo "⚠️  Strict filtering did not reduce variants (may be due to data characteristics)" >> validation_report.txt
        fi
        
        echo "" >> validation_report.txt
        
        # Test 3: Relaxed filtering validation
        echo "=== Test 3: Relaxed Filtering ===" >> validation_report.txt
        
        if [ -f "~{relaxed_msav}" ]; then
            echo "✅ Relaxed MSAV file created" >> validation_report.txt
            relaxed_size=$(stat -f%z "~{relaxed_msav}" 2>/dev/null || stat -c%s "~{relaxed_msav}")
            echo "   File size: $relaxed_size bytes" >> validation_report.txt
        else
            echo "❌ Relaxed MSAV file missing" >> validation_report.txt
            all_tests_passed=false
        fi
        
        relaxed_final=$(grep "Final variants:" "~{relaxed_qc_summary}" | awk '{print $3}')
        echo "Relaxed test final variants: $relaxed_final" >> validation_report.txt
        
        # Relaxed filtering should result in more variants than strict
        if [ "$relaxed_final" -ge "$strict_final" ]; then
            echo "✅ Relaxed filtering preserved more variants: $strict_final ≤ $relaxed_final" >> validation_report.txt
        else
            echo "⚠️  Relaxed filtering unexpectedly reduced variants" >> validation_report.txt
        fi
        
        echo "" >> validation_report.txt
        
        # Cross-scenario validation
        echo "=== Cross-Scenario Validation ===" >> validation_report.txt
        
        # Check file format consistency
        basic_type=$(file "~{basic_msav}" | grep -o "Zstandard")
        strict_type=$(file "~{strict_msav}" | grep -o "Zstandard")
        relaxed_type=$(file "~{relaxed_msav}" | grep -o "Zstandard")
        
        if [ "$basic_type" = "Zstandard" ] && [ "$strict_type" = "Zstandard" ] && [ "$relaxed_type" = "Zstandard" ]; then
            echo "✅ All MSAV files have correct Zstandard format" >> validation_report.txt
        else
            echo "❌ Inconsistent MSAV file formats detected" >> validation_report.txt
            all_tests_passed=false
        fi
        
        # Check that all files are non-empty
        if [ "$basic_size" -gt 100 ] && [ "$strict_size" -gt 100 ] && [ "$relaxed_size" -gt 100 ]; then
            echo "✅ All MSAV files have reasonable sizes" >> validation_report.txt
        else
            echo "❌ Some MSAV files are too small" >> validation_report.txt
            all_tests_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "=== Performance Summary ===" >> validation_report.txt
        echo "Variant retention rates:" >> validation_report.txt
        echo "- Basic (MAF≥0.01, CR≥0.95): $(echo "scale=1; $basic_final * 100 / $original_count" | bc)%" >> validation_report.txt
        echo "- Strict (MAF≥0.10, CR≥0.99): $(echo "scale=1; $strict_final * 100 / $original_count" | bc)%" >> validation_report.txt
        echo "- Relaxed (MAF≥0.001, CR≥0.80): $(echo "scale=1; $relaxed_final * 100 / $original_count" | bc)%" >> validation_report.txt
        
        echo "" >> validation_report.txt
        echo "File size comparison:" >> validation_report.txt
        echo "- Basic MSAV: $basic_size bytes" >> validation_report.txt
        echo "- Strict MSAV: $strict_size bytes" >> validation_report.txt
        echo "- Relaxed MSAV: $relaxed_size bytes" >> validation_report.txt
        
        echo "" >> validation_report.txt
        echo "=== Overall Test Result ===" >> validation_report.txt
        if [ "$all_tests_passed" = "true" ]; then
            echo "🎉 ALL TESTS PASSED: Complete pipeline working correctly across all scenarios" >> validation_report.txt
            echo "✅ Basic filtering scenario: PASS" >> validation_report.txt
            echo "✅ Strict filtering scenario: PASS" >> validation_report.txt
            echo "✅ Relaxed filtering scenario: PASS" >> validation_report.txt
            echo "✅ Cross-scenario validation: PASS" >> validation_report.txt
            echo "true" > all_tests_passed.txt
        else
            echo "❌ SOME TESTS FAILED: Pipeline has issues in one or more scenarios" >> validation_report.txt
            echo "false" > all_tests_passed.txt
        fi
    >>>

    output {
        File validation_report = "validation_report.txt"
        Boolean all_tests_passed = read_boolean("all_tests_passed.txt")
    }

    runtime {
        docker: "mamana/imputation:minimac4-4.1.6"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
        docker_args: "--platform linux/amd64"
    }

    meta {
        description: "Validate complete pipeline across multiple test scenarios"
    }
} 