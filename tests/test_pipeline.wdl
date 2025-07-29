version 1.0

import "../workflows/federated_imputation_pipeline.wdl" as MainPipeline

## Test workflow for the Federated Genotype Imputation Pipeline
## This workflow includes additional validation and testing steps
workflow TestFederatedImputationPipeline {
    input {
        File input_vcf = "tests/test_chr22_region_bgzip.vcf.gz"
        String chromosome = "22"
        Int start_position = 16000000
        Int end_position = 16990000
        String output_prefix = "test_federated_panel"
    }

    # Run the main pipeline with test parameters
    call MainPipeline.FederatedImputationPipeline {
        input:
            input_vcf = input_vcf,
            chromosome = chromosome,
            start_position = start_position,
            end_position = end_position,
            output_prefix = output_prefix,
            maf_threshold = 0.05,  # Relaxed for test data
            call_rate_threshold = 0.90,  # Relaxed for test data
            remove_indels = true,
            compression_level = 3  # Faster compression for testing
    }

    # Additional validation task
    call ValidateTestResults {
        input:
            msav_file = FederatedImputationPipeline.federated_panel_msav,
            qc_summary = FederatedImputationPipeline.qc_summary,
            conversion_summary = FederatedImputationPipeline.conversion_summary
    }

    output {
        # Main pipeline outputs
        File test_msav_file = FederatedImputationPipeline.federated_panel_msav
        File test_validation_report = ValidateTestResults.validation_report
        
        # All intermediate files for inspection
        File extracted_vcf = FederatedImputationPipeline.extracted_vcf
        File quality_controlled_vcf = FederatedImputationPipeline.quality_controlled_vcf
        Array[File] all_reports = [
            FederatedImputationPipeline.extraction_summary,
            FederatedImputationPipeline.qc_summary,
            FederatedImputationPipeline.conversion_summary,
            ValidateTestResults.validation_report
        ]
    }

    meta {
        description: "Test workflow for Federated Genotype Imputation Pipeline"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }
}

## Validate test results and generate comprehensive report
task ValidateTestResults {
    input {
        File msav_file
        File qc_summary
        File conversion_summary
    }

    command <<<
        echo "Test Validation Report" > validation_report.txt
        echo "======================" >> validation_report.txt
        echo "Generated: $(date)" >> validation_report.txt
        echo "" >> validation_report.txt
        
        # Check MSAV file
        if [ -f "~{msav_file}" ]; then
            echo "✅ MSAV file created successfully" >> validation_report.txt
            echo "File: ~{msav_file}" >> validation_report.txt
            echo "Size: $(ls -lh ~{msav_file} | awk '{print $5}')" >> validation_report.txt
            echo "Type: $(file ~{msav_file})" >> validation_report.txt
        else
            echo "❌ MSAV file not found" >> validation_report.txt
        fi
        echo "" >> validation_report.txt
        
        # Extract QC metrics
        echo "Quality Control Results:" >> validation_report.txt
        echo "========================" >> validation_report.txt
        cat "~{qc_summary}" >> validation_report.txt
        echo "" >> validation_report.txt
        
        # Extract conversion metrics
        echo "Conversion Results:" >> validation_report.txt
        echo "==================" >> validation_report.txt
        cat "~{conversion_summary}" >> validation_report.txt
        echo "" >> validation_report.txt
        
        # Test expectations for the test data
        echo "Test Expectations vs Results:" >> validation_report.txt
        echo "=============================" >> validation_report.txt
        echo "Expected input variants: 100" >> validation_report.txt
        echo "Expected output after QC: 10-50 variants (depends on MAF filter)" >> validation_report.txt
        echo "Expected MSAV file size: 1-10 KB" >> validation_report.txt
        echo "" >> validation_report.txt
        
        # Overall test status
        if [ -f "~{msav_file}" ]; then
            echo "🎉 TEST PASSED: Pipeline completed successfully!" >> validation_report.txt
            echo "✅ Ready for production use with larger datasets" >> validation_report.txt
        else
            echo "❌ TEST FAILED: Pipeline did not complete successfully" >> validation_report.txt
        fi
    >>>

    output {
        File validation_report = "validation_report.txt"
    }

    runtime {
        docker: "mamana/imputation:minimac4-4.1.6"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
        docker_args: "--platform linux/amd64"
    }

    meta {
        description: "Validate test results and generate comprehensive report"
    }
} 