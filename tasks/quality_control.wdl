version 1.0

## Apply quality control filters to VCF files
task QualityControl {
    input {
        File input_vcf
        File input_vcf_index
        Float maf_threshold = 0.01
        Float call_rate_threshold = 0.95
        Boolean remove_indels = true
        String output_prefix = "quality_controlled"
    }

    command <<<
        echo "Applying quality control filters..."
        echo "MAF threshold: ~{maf_threshold}"
        echo "Call rate threshold: ~{call_rate_threshold}"
        echo "Remove indels: ~{remove_indels}"
        
        # Start with the input VCF
        cp "~{input_vcf}" temp_input.vcf.gz
        cp "~{input_vcf_index}" temp_input.vcf.gz.csi
        
        # Apply filters step by step
        current_vcf="temp_input.vcf.gz"
        
        # Filter by MAF (Minor Allele Frequency)
        if [ "$(echo "~{maf_threshold} > 0" | bc)" -eq 1 ]; then
            echo "Filtering by MAF >= ~{maf_threshold}..."
            bcftools view \
                --min-af ~{maf_threshold} \
                --output-type z \
                --output temp_maf.vcf.gz \
                "$current_vcf"
            bcftools index temp_maf.vcf.gz
            current_vcf="temp_maf.vcf.gz"
        fi
        
        # Remove indels if requested (keep only SNPs)
        if [ "~{remove_indels}" = "true" ]; then
            echo "Removing indels (keeping only SNPs)..."
            bcftools view \
                --types snps \
                --output-type z \
                --output temp_snps.vcf.gz \
                "$current_vcf"
            bcftools index temp_snps.vcf.gz
            current_vcf="temp_snps.vcf.gz"
        fi
        
        # Apply call rate filter (remove variants with too many missing calls)
        if [ "$(echo "~{call_rate_threshold} < 1" | bc)" -eq 1 ]; then
            echo "Filtering by call rate >= ~{call_rate_threshold}..."
            # Calculate maximum missing genotypes allowed
            max_missing=$(echo "1 - ~{call_rate_threshold}" | bc -l)
            bcftools view \
                --exclude "F_MISSING > $max_missing" \
                --output-type z \
                --output temp_callrate.vcf.gz \
                "$current_vcf"
            bcftools index temp_callrate.vcf.gz
            current_vcf="temp_callrate.vcf.gz"
        fi
        
        # Final output
        cp "$current_vcf" "~{output_prefix}.vcf.gz"
        cp "${current_vcf}.csi" "~{output_prefix}.vcf.gz.csi"
        
        # Generate quality control statistics
        bcftools stats "~{output_prefix}.vcf.gz" > "~{output_prefix}_qc_stats.txt"
        
        # Generate summary report
        original_variants=$(bcftools view -H "~{input_vcf}" | wc -l)
        final_variants=$(bcftools view -H "~{output_prefix}.vcf.gz" | wc -l)
        
        echo "Quality Control Summary" > "~{output_prefix}_qc_summary.txt"
        echo "======================" >> "~{output_prefix}_qc_summary.txt"
        echo "Original variants: $original_variants" >> "~{output_prefix}_qc_summary.txt"
        echo "Final variants: $final_variants" >> "~{output_prefix}_qc_summary.txt"
        echo "Variants removed: $((original_variants - final_variants))" >> "~{output_prefix}_qc_summary.txt"
        echo "Retention rate: $(echo "scale=2; $final_variants * 100 / $original_variants" | bc)%" >> "~{output_prefix}_qc_summary.txt"
        echo "" >> "~{output_prefix}_qc_summary.txt"
        echo "Applied filters:" >> "~{output_prefix}_qc_summary.txt"
        echo "- MAF threshold: ~{maf_threshold}" >> "~{output_prefix}_qc_summary.txt"
        echo "- Call rate threshold: ~{call_rate_threshold}" >> "~{output_prefix}_qc_summary.txt"
        echo "- Remove indels: ~{remove_indels}" >> "~{output_prefix}_qc_summary.txt"
        
        # Cleanup temporary files
        rm -f temp_*.vcf.gz temp_*.vcf.gz.csi
    >>>

    output {
        File qc_vcf = "~{output_prefix}.vcf.gz"
        File qc_vcf_index = "~{output_prefix}.vcf.gz.csi"
        File qc_stats = "~{output_prefix}_qc_stats.txt"
        File qc_summary = "~{output_prefix}_qc_summary.txt"
    }

    runtime {
        docker: "mamana/imputation:minimac4-4.1.6"
        memory: "4 GB"
        cpu: 2
        disks: "local-disk 20 SSD"
        docker_args: "--platform linux/amd64"
    }

    meta {
        description: "Apply quality control filters to VCF files"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }

    parameter_meta {
        input_vcf: "Input VCF file to apply QC filters"
        input_vcf_index: "Index file for the input VCF"
        maf_threshold: "Minimum allele frequency threshold (0.0-1.0)"
        call_rate_threshold: "Minimum call rate threshold (0.0-1.0)"
        remove_indels: "Remove indels and keep only SNPs"
        output_prefix: "Prefix for output files"
    }
} 