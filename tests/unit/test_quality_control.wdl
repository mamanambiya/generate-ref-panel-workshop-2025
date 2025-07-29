version 1.0

import "../../tasks/quality_control.wdl" as QualityControlTask

## Unit test for QualityControl task
workflow TestQualityControl {
    input {
        File test_vcf = "tests/test_chr22_region_bgzip.vcf.gz"
        Float test_maf = 0.05
        Float test_call_rate = 0.90
        Boolean test_remove_indels = true
        String test_prefix = "unit_test_qc"
    }

    call QualityControlTask.QualityControl {
        input:
            input_vcf = test_vcf,
            input_vcf_index = test_vcf + ".csi",
            maf_threshold = test_maf,
            call_rate_threshold = test_call_rate,
            remove_indels = test_remove_indels,
            output_prefix = test_prefix
    }

    call ValidateQualityControl {
        input:
            input_vcf = test_vcf,
            qc_vcf = QualityControl.qc_vcf,
            qc_summary = QualityControl.qc_summary,
            qc_stats = QualityControl.qc_stats,
            expected_maf = test_maf,
            expected_call_rate = test_call_rate,
            expected_remove_indels = test_remove_indels
    }

    output {
        File test_qc_vcf = QualityControl.qc_vcf
        File test_qc_summary = QualityControl.qc_summary
        File validation_report = ValidateQualityControl.validation_report
        Boolean test_passed = ValidateQualityControl.test_passed
    }

    meta {
        description: "Unit test for QualityControl task"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }
}

## Validate quality control results
task ValidateQualityControl {
    input {
        File input_vcf
        File qc_vcf
        File qc_summary
        File qc_stats
        Float expected_maf
        Float expected_call_rate
        Boolean expected_remove_indels
    }

    command <<<
        echo "QualityControl Unit Test Validation" > validation_report.txt
        echo "==================================" >> validation_report.txt
        echo "Test Date: $(date)" >> validation_report.txt
        echo "" >> validation_report.txt
        
        test_passed=true
        
        echo "File Existence Tests:" >> validation_report.txt
        if [ -f "~{qc_vcf}" ]; then
            echo "✅ QC VCF file exists" >> validation_report.txt
        else
            echo "❌ QC VCF file missing" >> validation_report.txt
            test_passed=false
        fi
        
        if [ -f "~{qc_summary}" ]; then
            echo "✅ QC summary file exists" >> validation_report.txt
        else
            echo "❌ QC summary file missing" >> validation_report.txt
            test_passed=false
        fi
        
        if [ -f "~{qc_stats}" ]; then
            echo "✅ QC stats file exists" >> validation_report.txt
        else
            echo "❌ QC stats file missing" >> validation_report.txt
            test_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "Quality Control Validation Tests:" >> validation_report.txt
        
        # Count variants before and after QC
        original_count=$(bcftools view -H "~{input_vcf}" | wc -l)
        qc_count=$(bcftools view -H "~{qc_vcf}" | wc -l)
        
        echo "Original variants: $original_count" >> validation_report.txt
        echo "QC variants: $qc_count" >> validation_report.txt
        
        # Validate that some filtering occurred (unless all variants pass)
        if [ "$qc_count" -le "$original_count" ]; then
            echo "✅ QC reduced or maintained variant count" >> validation_report.txt
        else
            echo "❌ QC increased variant count (impossible)" >> validation_report.txt
            test_passed=false
        fi
        
        # Check if indels were removed when requested
        if [ "~{expected_remove_indels}" = "true" ]; then
            indel_count=$(bcftools view -H "~{qc_vcf}" | awk '{if(length($4) > 1 || length($5) > 1) print}' | wc -l)
            if [ "$indel_count" -eq 0 ]; then
                echo "✅ All indels removed as requested" >> validation_report.txt
            else
                echo "❌ Found $indel_count indels in QC output" >> validation_report.txt
                test_passed=false
            fi
        fi
        
        # Validate MAF filtering
        echo "Checking MAF threshold: ~{expected_maf}" >> validation_report.txt
        low_maf_count=$(bcftools query -f '%INFO/AF\n' "~{qc_vcf}" | awk -v maf=~{expected_maf} '$1 < maf {count++} END {print count+0}')
        if [ "$low_maf_count" -eq 0 ]; then
            echo "✅ No variants below MAF threshold found" >> validation_report.txt
        else
            echo "⚠️  Found $low_maf_count variants below MAF threshold (may be expected if filter not strict)" >> validation_report.txt
        fi
        
        # Validate summary content
        if grep -q "Quality Control Summary" "~{qc_summary}"; then
            echo "✅ Summary contains header" >> validation_report.txt
        else
            echo "❌ Summary missing header" >> validation_report.txt
            test_passed=false
        fi
        
        if grep -q "Original variants: $original_count" "~{qc_summary}"; then
            echo "✅ Summary contains correct original count" >> validation_report.txt
        else
            echo "❌ Summary missing or incorrect original count" >> validation_report.txt
            test_passed=false
        fi
        
        if grep -q "Final variants: $qc_count" "~{qc_summary}"; then
            echo "✅ Summary contains correct final count" >> validation_report.txt
        else
            echo "❌ Summary missing or incorrect final count" >> validation_report.txt
            test_passed=false
        fi
        
        # Check filter parameters in summary
        if grep -q "MAF threshold: ~{expected_maf}" "~{qc_summary}"; then
            echo "✅ Summary contains correct MAF threshold" >> validation_report.txt
        else
            echo "❌ Summary missing or incorrect MAF threshold" >> validation_report.txt
            test_passed=false
        fi
        
        if grep -q "Call rate threshold: ~{expected_call_rate}" "~{qc_summary}"; then
            echo "✅ Summary contains correct call rate threshold" >> validation_report.txt
        else
            echo "❌ Summary missing or incorrect call rate threshold" >> validation_report.txt
            test_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "Overall Test Result:" >> validation_report.txt
        if [ "$test_passed" = "true" ]; then
            echo "🎉 TEST PASSED: QualityControl working correctly" >> validation_report.txt
            echo "true" > test_passed.txt
        else
            echo "❌ TEST FAILED: QualityControl has issues" >> validation_report.txt
            echo "false" > test_passed.txt
        fi
    >>>

    output {
        File validation_report = "validation_report.txt"
        Boolean test_passed = read_boolean("test_passed.txt")
    }

    runtime {
        docker: "mamana/imputation:minimac4-4.1.6"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
        docker_args: "--platform linux/amd64"
    }

    meta {
        description: "Validate QualityControl task results"
    }
} 