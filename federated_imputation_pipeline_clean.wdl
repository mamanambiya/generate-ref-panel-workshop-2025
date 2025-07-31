version 1.0

## Federated Genotype Imputation Pipeline - Clean Version
## GA4GH Hackathon 2025 - African Genomics Team
##
## This workflow processes VCF files into MSAV format to enable federated genotype
## imputation networks while preserving data sovereignty.
##
## Clean version designed for engine-level Docker platform handling

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
        docker: "mamana/minimac4-all"
        memory: "4 GB"
        cpu: 2
        disks: "local-disk 20 SSD"
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

## Apply quality control filters to VCF files
task QualityControl {
    input {
        File input_vcf
        File input_vcf_index
        Float min_quality = 20.0
        Float min_af = 0.01
        Float max_missing = 0.1
        String output_prefix = "qc_filtered"
    }

    command <<<
        echo "Applying quality control filters..."
        echo "Min Quality: ~{min_quality}"
        echo "Min Allele Frequency: ~{min_af}"
        echo "Max Missing Rate: ~{max_missing}"
        
        # Apply quality filters
        bcftools view \
            --include "QUAL>=~{min_quality} && AF>=~{min_af} && F_MISSING<=~{max_missing}" \
            --output-type z \
            --output "~{output_prefix}.vcf.gz" \
            "~{input_vcf}"
        
        # Index the filtered VCF
        bcftools index "~{output_prefix}.vcf.gz"
        
        # Generate quality control statistics
        bcftools stats "~{output_prefix}.vcf.gz" > "~{output_prefix}_qc_stats.txt"
        
        # Count variants after QC
        variant_count=$(bcftools view -H "~{output_prefix}.vcf.gz" | wc -l)
        original_count=$(bcftools view -H "~{input_vcf}" | wc -l)
        
        echo "Original variants: ${original_count}" > "~{output_prefix}_qc_summary.txt"
        echo "Post-QC variants: ${variant_count}" >> "~{output_prefix}_qc_summary.txt"
        echo "Filtered out: $((original_count - variant_count))" >> "~{output_prefix}_qc_summary.txt"
        echo "Quality control completed successfully" >> "~{output_prefix}_qc_summary.txt"
    >>>

    output {
        File qc_vcf = "~{output_prefix}.vcf.gz"
        File qc_vcf_index = "~{output_prefix}.vcf.gz.csi"
        File qc_stats = "~{output_prefix}_qc_stats.txt"
        File qc_summary = "~{output_prefix}_qc_summary.txt"
    }

    runtime {
        docker: "mamana/minimac4-all"
        memory: "4 GB"
        cpu: 2
        disks: "local-disk 20 SSD"
    }

    meta {
        description: "Apply quality control filters to VCF files"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }

    parameter_meta {
        input_vcf: "Input VCF file to apply QC filters"
        input_vcf_index: "Index file for the input VCF"
        min_quality: "Minimum variant quality score (default: 20.0)"
        min_af: "Minimum allele frequency (default: 0.01)"
        max_missing: "Maximum missing rate (default: 0.1)"
        output_prefix: "Prefix for output files"
    }
}

## Convert quality-controlled VCF to MSAV format using Minimac4
task MinimacConversion {
    input {
        File input_vcf
        File input_vcf_index
        String output_prefix = "converted"
        Boolean compress_reference = true
    }

    command <<<
        echo "Converting VCF to MSAV format using Minimac4..."
        
        # Ensure input VCF is properly formatted
        echo "Input VCF: ~{input_vcf}"
        
        # Convert to MSAV format using Minimac4
        if [ "~{compress_reference}" = "true" ]; then
            minimac4 \
                --compress-reference \
                --refHaps "~{input_vcf}" \
                --output "~{output_prefix}"
        else
            minimac4 \
                --refHaps "~{input_vcf}" \
                --output "~{output_prefix}"
        fi
        
        # Generate conversion summary
        echo "VCF to MSAV conversion completed" > "~{output_prefix}_conversion_summary.txt"
        echo "Input: ~{input_vcf}" >> "~{output_prefix}_conversion_summary.txt"
        echo "Output prefix: ~{output_prefix}" >> "~{output_prefix}_conversion_summary.txt"
        echo "Compress reference: ~{compress_reference}" >> "~{output_prefix}_conversion_summary.txt"
        
        # List all generated files
        echo "Generated files:" >> "~{output_prefix}_conversion_summary.txt"
        ls -la ~{output_prefix}.* >> "~{output_prefix}_conversion_summary.txt" || true
    >>>

    output {
        File msav_file = "~{output_prefix}.msav"
        File conversion_summary = "~{output_prefix}_conversion_summary.txt"
        # Minimac4 may generate additional files
        Array[File] all_output_files = glob("~{output_prefix}.*")
    }

    runtime {
        docker: "mamana/minimac4-all"
        memory: "6 GB"
        cpu: 2
        disks: "local-disk 30 SSD"
    }

    meta {
        description: "Convert quality-controlled VCF to MSAV format using Minimac4"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
    }

    parameter_meta {
        input_vcf: "Quality-controlled VCF file for conversion"
        input_vcf_index: "Index file for the input VCF"
        output_prefix: "Prefix for output MSAV files"
        compress_reference: "Whether to compress the reference for Minimac4 (default: true)"
    }
}

## Complete federated imputation pipeline workflow
workflow FederatedImputationPipeline {
    input {
        File input_vcf
        File input_vcf_index
        String chromosome = "22"
        Int start_position = 16000000
        Int end_position = 16990000
        Float min_quality = 20.0
        Float min_af = 0.01
        Float max_missing = 0.1
        String output_prefix = "federated_imputation"
        Boolean compress_reference = true
    }

    # Step 1: Extract genomic region
    call ExtractRegion {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            chromosome = chromosome,
            start_position = start_position,
            end_position = end_position,
            output_prefix = output_prefix + "_extracted"
    }

    # Step 2: Apply quality control
    call QualityControl {
        input:
            input_vcf = ExtractRegion.extracted_vcf,
            input_vcf_index = ExtractRegion.extracted_vcf_index,
            min_quality = min_quality,
            min_af = min_af,
            max_missing = max_missing,
            output_prefix = output_prefix + "_qc"
    }

    # Step 3: Convert to MSAV format
    call MinimacConversion {
        input:
            input_vcf = QualityControl.qc_vcf,
            input_vcf_index = QualityControl.qc_vcf_index,
            output_prefix = output_prefix + "_msav",
            compress_reference = compress_reference
    }

    output {
        # Final MSAV output for federated imputation
        File msav_output = MinimacConversion.msav_file
        Array[File] all_msav_files = MinimacConversion.all_output_files
        
        # Intermediate files for debugging/validation
        File extracted_vcf = ExtractRegion.extracted_vcf
        File qc_vcf = QualityControl.qc_vcf
        
        # Summary and statistics files
        File extraction_summary = ExtractRegion.summary_file
        File extraction_stats = ExtractRegion.stats_file
        File qc_summary = QualityControl.qc_summary
        File qc_stats = QualityControl.qc_stats
        File conversion_summary = MinimacConversion.conversion_summary
    }

    meta {
        description: "Federated Genotype Imputation Pipeline - Clean Version (Engine-Level Platform Handling)"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
        email: "team@afrigenomics.org"
        version: "1.0.0"
        docker: "mamana/minimac4-all"
        execution_notes: "Use with cromwell.conf or run_cromwell_multiplatform.sh for platform handling"
    }

    parameter_meta {
        input_vcf: "Input VCF file containing genetic variants"
        input_vcf_index: "Index file for the input VCF"
        chromosome: "Chromosome to extract (default: 22)"
        start_position: "Start position for extraction (default: 16000000)"
        end_position: "End position for extraction (default: 16990000)"
        min_quality: "Minimum variant quality score (default: 20.0)"
        min_af: "Minimum allele frequency (default: 0.01)"
        max_missing: "Maximum missing rate (default: 0.1)"
        output_prefix: "Prefix for all output files (default: federated_imputation)"
        compress_reference: "Whether to compress the reference for Minimac4 (default: true)"
        msav_output: "Final MSAV file for federated imputation networks"
        all_msav_files: "All files generated by Minimac4 conversion"
    }

} 