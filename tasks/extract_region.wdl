version 1.0

## Extract genomic regions from VCF files using bcftools
task ExtractRegion {
    input {
        File input_vcf
        File input_vcf_index
        String chromosome
        Int start_position
        Int end_position
        String output_prefix = "extracted_region"
    }

    command <<<
        echo "Extracting region chr~{chromosome}:~{start_position}-~{end_position}..."
        
        bcftools view \
            --regions "~{chromosome}:~{start_position}-~{end_position}" \
            --output-type z \
            --output "~{output_prefix}.vcf.gz" \
            "~{input_vcf}"
        
        # Index the output VCF
        bcftools index "~{output_prefix}.vcf.gz"
        
        # Generate statistics
        bcftools stats "~{output_prefix}.vcf.gz" > "~{output_prefix}_stats.txt"
        
        # Count extracted variants
        variant_count=$(bcftools view -H "~{output_prefix}.vcf.gz" | wc -l)
        echo "Extracted ${variant_count} variants" > "~{output_prefix}_summary.txt"
        echo "Region: chr~{chromosome}:~{start_position}-~{end_position}" >> "~{output_prefix}_summary.txt"
        echo "Extraction completed successfully" >> "~{output_prefix}_summary.txt"
    >>>

    output {
        File extracted_vcf = "~{output_prefix}.vcf.gz"
        File extracted_vcf_index = "~{output_prefix}.vcf.gz.csi"
        File stats_file = "~{output_prefix}_stats.txt"
        File summary_file = "~{output_prefix}_summary.txt"
    }

    runtime {
        docker: "mamana/imputation:minimac4-4.1.6"
        memory: "4 GB"
        cpu: 2
        disks: "local-disk 20 SSD"
        docker_args: "--platform linux/amd64"
    }

    meta {
        description: "Extract genomic regions from VCF files using bcftools"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }

    parameter_meta {
        input_vcf: "Input VCF file (bgzip compressed)"
        input_vcf_index: "Index file for the input VCF"
        chromosome: "Target chromosome (e.g., '22')"
        start_position: "Start genomic position"
        end_position: "End genomic position"
        output_prefix: "Prefix for output files"
    }
} 