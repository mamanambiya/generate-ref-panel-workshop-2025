version 1.0

import "../../tasks/extract_region.wdl" as ExtractRegionTask

## Unit test for ExtractRegion task
workflow TestExtractRegion {
    input {
        File test_vcf = "tests/test_chr22_region_bgzip.vcf.gz"
        String test_chromosome = "22"
        Int test_start = 16000000
        Int test_end = 16200000
        String test_prefix = "unit_test_extracted"
    }

    call ExtractRegionTask.ExtractRegion {
        input:
            input_vcf = test_vcf,
            input_vcf_index = test_vcf + ".csi",
            chromosome = test_chromosome,
            start_position = test_start,
            end_position = test_end,
            output_prefix = test_prefix
    }

    call ValidateExtraction {
        input:
            extracted_vcf = ExtractRegion.extracted_vcf,
            summary_file = ExtractRegion.summary_file,
            stats_file = ExtractRegion.stats_file,
            expected_chromosome = test_chromosome,
            expected_start = test_start,
            expected_end = test_end
    }

    output {
        File test_extracted_vcf = ExtractRegion.extracted_vcf
        File test_summary = ExtractRegion.summary_file
        File validation_report = ValidateExtraction.validation_report
        Boolean test_passed = ValidateExtraction.test_passed
    }

    meta {
        description: "Unit test for ExtractRegion task"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }
}

## Validate extraction results
task ValidateExtraction {
    input {
        File extracted_vcf
        File summary_file
        File stats_file
        String expected_chromosome
        Int expected_start
        Int expected_end
    }

    command <<<
        echo "ExtractRegion Unit Test Validation" > validation_report.txt
        echo "=================================" >> validation_report.txt
        echo "Test Date: $(date)" >> validation_report.txt
        echo "" >> validation_report.txt
        
        # Check if files exist
        test_passed=true
        
        echo "File Existence Tests:" >> validation_report.txt
        if [ -f "~{extracted_vcf}" ]; then
            echo "✅ Extracted VCF file exists" >> validation_report.txt
        else
            echo "❌ Extracted VCF file missing" >> validation_report.txt
            test_passed=false
        fi
        
        if [ -f "~{summary_file}" ]; then
            echo "✅ Summary file exists" >> validation_report.txt
        else
            echo "❌ Summary file missing" >> validation_report.txt
            test_passed=false
        fi
        
        if [ -f "~{stats_file}" ]; then
            echo "✅ Stats file exists" >> validation_report.txt
        else
            echo "❌ Stats file missing" >> validation_report.txt
            test_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "Content Validation Tests:" >> validation_report.txt
        
        # Validate VCF content
        variant_count=$(bcftools view -H "~{extracted_vcf}" | wc -l)
        echo "Variant count: $variant_count" >> validation_report.txt
        
        if [ "$variant_count" -gt 0 ]; then
            echo "✅ VCF contains variants" >> validation_report.txt
        else
            echo "❌ VCF is empty" >> validation_report.txt
            test_passed=false
        fi
        
        # Check chromosome consistency
        chromosomes=$(bcftools view -H "~{extracted_vcf}" | cut -f1 | sort -u)
        if [ "$chromosomes" = "~{expected_chromosome}" ]; then
            echo "✅ Chromosome matches expected (~{expected_chromosome})" >> validation_report.txt
        else
            echo "❌ Chromosome mismatch. Expected: ~{expected_chromosome}, Found: $chromosomes" >> validation_report.txt
            test_passed=false
        fi
        
        # Check position range
        min_pos=$(bcftools view -H "~{extracted_vcf}" | cut -f2 | sort -n | head -1)
        max_pos=$(bcftools view -H "~{extracted_vcf}" | cut -f2 | sort -n | tail -1)
        
        if [ "$min_pos" -ge "~{expected_start}" ] && [ "$max_pos" -le "~{expected_end}" ]; then
            echo "✅ Positions within expected range (~{expected_start}-~{expected_end})" >> validation_report.txt
            echo "   Actual range: $min_pos-$max_pos" >> validation_report.txt
        else
            echo "❌ Positions outside expected range" >> validation_report.txt
            echo "   Expected: ~{expected_start}-~{expected_end}" >> validation_report.txt
            echo "   Actual: $min_pos-$max_pos" >> validation_report.txt
            test_passed=false
        fi
        
        # Validate summary content
        if grep -q "Extracted.*variants" "~{summary_file}"; then
            echo "✅ Summary contains variant count" >> validation_report.txt
        else
            echo "❌ Summary missing variant count" >> validation_report.txt
            test_passed=false
        fi
        
        if grep -q "chr~{expected_chromosome}:~{expected_start}-~{expected_end}" "~{summary_file}"; then
            echo "✅ Summary contains correct region" >> validation_report.txt
        else
            echo "❌ Summary missing or incorrect region" >> validation_report.txt
            test_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "Overall Test Result:" >> validation_report.txt
        if [ "$test_passed" = "true" ]; then
            echo "🎉 TEST PASSED: ExtractRegion working correctly" >> validation_report.txt
            echo "true" > test_passed.txt
        else
            echo "❌ TEST FAILED: ExtractRegion has issues" >> validation_report.txt
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
        description: "Validate ExtractRegion task results"
    }
} 