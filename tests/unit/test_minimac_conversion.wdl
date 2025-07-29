version 1.0

import "../../tasks/minimac_conversion.wdl" as MinimacConversionTask

## Unit test for MinimacConversion task
workflow TestMinimacConversion {
    input {
        File test_vcf = "tests/test_chr22_region_bgzip.vcf.gz"
        String test_prefix = "unit_test_minimac"
        Int test_compression = 3
    }

    call MinimacConversionTask.MinimacConversion {
        input:
            input_vcf = test_vcf,
            input_vcf_index = test_vcf + ".csi",
            output_prefix = test_prefix,
            compression_level = test_compression
    }

    call ValidateMinimacConversion {
        input:
            input_vcf = test_vcf,
            msav_file = MinimacConversion.msav_file,
            conversion_summary = MinimacConversion.conversion_summary,
            all_output_files = MinimacConversion.all_output_files,
            expected_prefix = test_prefix
    }

    output {
        File test_msav_file = MinimacConversion.msav_file
        File test_conversion_summary = MinimacConversion.conversion_summary
        Array[File] test_all_files = MinimacConversion.all_output_files
        File validation_report = ValidateMinimacConversion.validation_report
        Boolean test_passed = ValidateMinimacConversion.test_passed
    }

    meta {
        description: "Unit test for MinimacConversion task"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }
}

## Validate Minimac conversion results
task ValidateMinimacConversion {
    input {
        File input_vcf
        File msav_file
        File conversion_summary
        Array[File] all_output_files
        String expected_prefix
    }

    command <<<
        echo "MinimacConversion Unit Test Validation" > validation_report.txt
        echo "=====================================" >> validation_report.txt
        echo "Test Date: $(date)" >> validation_report.txt
        echo "" >> validation_report.txt
        
        test_passed=true
        
        echo "File Existence Tests:" >> validation_report.txt
        if [ -f "~{msav_file}" ]; then
            echo "✅ MSAV file exists" >> validation_report.txt
        else
            echo "❌ MSAV file missing" >> validation_report.txt
            test_passed=false
        fi
        
        if [ -f "~{conversion_summary}" ]; then
            echo "✅ Conversion summary exists" >> validation_report.txt
        else
            echo "❌ Conversion summary missing" >> validation_report.txt
            test_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "MSAV File Validation Tests:" >> validation_report.txt
        
        # Check file format
        file_type=$(file "~{msav_file}")
        if echo "$file_type" | grep -q "Zstandard"; then
            echo "✅ MSAV file has correct Zstandard compression" >> validation_report.txt
        else
            echo "❌ MSAV file does not have Zstandard compression" >> validation_report.txt
            echo "   File type: $file_type" >> validation_report.txt
            test_passed=false
        fi
        
        # Check file size (should be reasonable)
        file_size=$(stat -f%z "~{msav_file}" 2>/dev/null || stat -c%s "~{msav_file}")
        if [ "$file_size" -gt 100 ]; then
            echo "✅ MSAV file has reasonable size: $file_size bytes" >> validation_report.txt
        else
            echo "❌ MSAV file too small: $file_size bytes" >> validation_report.txt
            test_passed=false
        fi
        
        # Check filename matches expected prefix
        filename=$(basename "~{msav_file}")
        expected_filename="~{expected_prefix}.msav"
        if [ "$filename" = "$expected_filename" ]; then
            echo "✅ MSAV filename matches expected: $expected_filename" >> validation_report.txt
        else
            echo "❌ MSAV filename mismatch. Expected: $expected_filename, Got: $filename" >> validation_report.txt
            test_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "Summary Content Validation:" >> validation_report.txt
        
        # Count input variants
        input_count=$(bcftools view -H "~{input_vcf}" | wc -l)
        
        # Validate summary content
        if grep -q "MSAV Conversion Summary" "~{conversion_summary}"; then
            echo "✅ Summary contains header" >> validation_report.txt
        else
            echo "❌ Summary missing header" >> validation_report.txt
            test_passed=false
        fi
        
        if grep -q "Input variants: $input_count" "~{conversion_summary}"; then
            echo "✅ Summary contains correct input variant count" >> validation_report.txt
        else
            echo "❌ Summary missing or incorrect input variant count" >> validation_report.txt
            test_passed=false
        fi
        
        if grep -q "Output MSAV file: ~{expected_prefix}.msav" "~{conversion_summary}"; then
            echo "✅ Summary contains correct output filename" >> validation_report.txt
        else
            echo "❌ Summary missing or incorrect output filename" >> validation_report.txt
            test_passed=false
        fi
        
        if grep -q "Conversion completed successfully" "~{conversion_summary}"; then
            echo "✅ Summary indicates successful conversion" >> validation_report.txt
        else
            echo "❌ Summary does not indicate successful conversion" >> validation_report.txt
            test_passed=false
        fi
        
        # Check file type is mentioned in summary
        if grep -q "File type:.*Zstandard" "~{conversion_summary}"; then
            echo "✅ Summary contains file type validation" >> validation_report.txt
        else
            echo "⚠️  Summary missing file type validation" >> validation_report.txt
        fi
        
        echo "" >> validation_report.txt
        echo "Output Files Validation:" >> validation_report.txt
        
        # Count output files
        output_file_count=$(echo '~{sep=" " all_output_files}' | wc -w)
        echo "Output file count: $output_file_count" >> validation_report.txt
        
        if [ "$output_file_count" -ge 1 ]; then
            echo "✅ At least one output file generated" >> validation_report.txt
        else
            echo "❌ No output files generated" >> validation_report.txt
            test_passed=false
        fi
        
        echo "" >> validation_report.txt
        echo "Overall Test Result:" >> validation_report.txt
        if [ "$test_passed" = "true" ]; then
            echo "🎉 TEST PASSED: MinimacConversion working correctly" >> validation_report.txt
            echo "true" > test_passed.txt
        else
            echo "❌ TEST FAILED: MinimacConversion has issues" >> validation_report.txt
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
        description: "Validate MinimacConversion task results"
    }
} 