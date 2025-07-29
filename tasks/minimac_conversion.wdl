version 1.0

## Convert quality-controlled VCF to MSAV format using Minimac4
task MinimacConversion {
    input {
        File input_vcf
        File input_vcf_index
        String output_prefix = "federated_panel"
        Int compression_level = 5
    }

    command <<<
        echo "Converting VCF to MSAV format using Minimac4..."
        echo "Input VCF: ~{input_vcf}"
        echo "Output prefix: ~{output_prefix}"
        echo "Compression level: ~{compression_level}"
        
        # Check input file
        variant_count=$(bcftools view -H "~{input_vcf}" | wc -l)
        echo "Processing $variant_count variants..."
        
        # Convert VCF to MSAV using Minimac4
        minimac4 \
            --compress-reference "~{input_vcf}" \
            --output "~{output_prefix}.msav"
        
        # Generate conversion summary
        echo "MSAV Conversion Summary" > "~{output_prefix}_conversion_summary.txt"
        echo "======================" >> "~{output_prefix}_conversion_summary.txt"
        echo "Input variants: $variant_count" >> "~{output_prefix}_conversion_summary.txt"
        echo "Compression level: ~{compression_level}" >> "~{output_prefix}_conversion_summary.txt"
        echo "Output MSAV file: ~{output_prefix}.msav" >> "~{output_prefix}_conversion_summary.txt"
        
        # Check if MSAV file was created successfully
        if [ -f "~{output_prefix}.msav" ]; then
            msav_size=$(ls -lh "~{output_prefix}.msav" | awk '{print $5}')
            echo "MSAV file size: $msav_size" >> "~{output_prefix}_conversion_summary.txt"
            echo "Conversion completed successfully" >> "~{output_prefix}_conversion_summary.txt"
            
            # Validate MSAV file format
            file_type=$(file "~{output_prefix}.msav")
            echo "File type: $file_type" >> "~{output_prefix}_conversion_summary.txt"
        else
            echo "ERROR: MSAV file was not created" >> "~{output_prefix}_conversion_summary.txt"
            exit 1
        fi
        
        # List all generated files
        echo "" >> "~{output_prefix}_conversion_summary.txt"
        echo "Generated files:" >> "~{output_prefix}_conversion_summary.txt"
        ls -la ~{output_prefix}.* >> "~{output_prefix}_conversion_summary.txt"
    >>>

    output {
        File msav_file = "~{output_prefix}.msav"
        File conversion_summary = "~{output_prefix}_conversion_summary.txt"
        # Minimac4 may generate additional files
        Array[File] all_output_files = glob("~{output_prefix}.*")
    }

    runtime {
        docker: "mamana/imputation:minimac4-4.1.6"
        memory: "6 GB"
        cpu: 2
        disks: "local-disk 30 SSD"
        docker_args: "--platform linux/amd64"
    }

    meta {
        description: "Convert quality-controlled VCF to MSAV format using Minimac4"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }

    parameter_meta {
        input_vcf: "Quality-controlled VCF file for conversion"
        input_vcf_index: "Index file for the input VCF"
        output_prefix: "Prefix for output MSAV files"
        compression_level: "MSAV compression level (1-9, default: 5)"
    }
} 