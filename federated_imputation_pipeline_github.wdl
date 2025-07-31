version 1.0

## Federated Genotype Imputation Pipeline - GitHub Container Registry Version
## GA4GH Hackathon 2025 - African Genomics Team
##
## This workflow processes VCF files into MSAV format to enable federated genotype
## imputation networks while preserving data sovereignty.
## Uses GitHub Container Registry: ghcr.io/mamanambiya/federated-imputation:latest

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
        docker: "ghcr.io/mamanambiya/federated-imputation:latest"
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
        docker: "ghcr.io/mamanambiya/federated-imputation:latest"
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
        maf_threshold: "Minimum allele frequency threshold (0.0-1.0)"
        call_rate_threshold: "Minimum call rate threshold (0.0-1.0)"
        remove_indels: "Remove indels and keep only SNPs"
        output_prefix: "Prefix for output files"
    }
}

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
        docker: "ghcr.io/mamanambiya/federated-imputation:latest"
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
        compression_level: "MSAV compression level (1-9, default: 5)"
    }
}

## Main workflow that orchestrates all tasks
workflow FederatedImputationPipeline {
    input {
        # Required inputs
        File input_vcf
        String chromosome
        Int start_position  
        Int end_position
        
        # Optional parameters with defaults
        String output_prefix = "federated_panel"
        Float maf_threshold = 0.01
        Float call_rate_threshold = 0.95
        Boolean remove_indels = true
        Int compression_level = 5
    }

    # Step 1: Extract genomic region from input VCF
    call ExtractRegion {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf + ".csi",
            chromosome = chromosome,
            start_position = start_position,
            end_position = end_position,
            output_prefix = output_prefix + "_extracted"
    }

    # Step 2: Apply quality control filters
    call QualityControl {
        input:
            input_vcf = ExtractRegion.extracted_vcf,
            input_vcf_index = ExtractRegion.extracted_vcf_index,
            maf_threshold = maf_threshold,
            call_rate_threshold = call_rate_threshold,
            remove_indels = remove_indels,
            output_prefix = output_prefix + "_qc"
    }

    # Step 3: Convert to MSAV format using Minimac4
    call MinimacConversion {
        input:
            input_vcf = QualityControl.qc_vcf,
            input_vcf_index = QualityControl.qc_vcf_index,
            output_prefix = output_prefix,
            compression_level = compression_level
    }

    # Pipeline outputs
    output {
        # Final MSAV output for federated imputation
        File federated_panel_msav = MinimacConversion.msav_file
        Array[File] all_msav_files = MinimacConversion.all_output_files
        
        # Intermediate VCF files
        File extracted_vcf = ExtractRegion.extracted_vcf
        File quality_controlled_vcf = QualityControl.qc_vcf
        
        # Statistics and reports
        File extraction_stats = ExtractRegion.stats_file
        File extraction_summary = ExtractRegion.summary_file
        File qc_stats = QualityControl.qc_stats
        File qc_summary = QualityControl.qc_summary
        File conversion_summary = MinimacConversion.conversion_summary
    }

    meta {
        description: "Federated Genotype Imputation Pipeline - GitHub Container Registry Version"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
        email: "team@afrigenomics.org"
        version: "1.0.0"
        docker: "ghcr.io/mamanambiya/federated-imputation:latest"
    }

    parameter_meta {
        input_vcf: {
            description: "Input VCF file (bgzip compressed and indexed)",
            help: "Must be bgzip compressed (.vcf.gz) with corresponding index (.csi)",
            pattern: "*.vcf.gz"
        }
        chromosome: {
            description: "Target chromosome for extraction",
            help: "Chromosome identifier (e.g., '22', 'X', 'Y')",
            example: "22"
        }
        start_position: {
            description: "Start genomic position for extraction",
            help: "1-based genomic coordinate",
            example: 16000000
        }
        end_position: {
            description: "End genomic position for extraction", 
            help: "1-based genomic coordinate",
            example: 16200000
        }
        output_prefix: {
            description: "Prefix for all output files",
            help: "Used to name intermediate and final output files",
            example: "h3africa_chr22_panel"
        }
        maf_threshold: {
            description: "Minimum allele frequency threshold",
            help: "Variants with MAF below this threshold will be removed",
            example: 0.01
        }
        call_rate_threshold: {
            description: "Minimum call rate threshold",
            help: "Variants with call rate below this threshold will be removed",
            example: 0.95
        }
        remove_indels: {
            description: "Remove indels and keep only SNPs",
            help: "If true, only single nucleotide polymorphisms are retained",
            example: true
        }
        compression_level: {
            description: "MSAV compression level (1-9)",
            help: "Higher values produce smaller files but take longer to process",
            example: 5
        }
    }
} 