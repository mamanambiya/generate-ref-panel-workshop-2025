version 1.0

import "../tasks/extract_region.wdl" as ExtractRegionTask
import "../tasks/quality_control.wdl" as QualityControlTask  
import "../tasks/minimac_conversion.wdl" as MinimacConversionTask

## Federated Genotype Imputation Pipeline
## GA4GH Hackathon 2025 - African Genomics Team
##
## This workflow processes VCF files into MSAV format to enable federated genotype
## imputation networks while preserving data sovereignty.
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
    call ExtractRegionTask.ExtractRegion {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf + ".csi",
            chromosome = chromosome,
            start_position = start_position,
            end_position = end_position,
            output_prefix = output_prefix + "_extracted"
    }

    # Step 2: Apply quality control filters
    call QualityControlTask.QualityControl {
        input:
            input_vcf = ExtractRegion.extracted_vcf,
            input_vcf_index = ExtractRegion.extracted_vcf_index,
            maf_threshold = maf_threshold,
            call_rate_threshold = call_rate_threshold,
            remove_indels = remove_indels,
            output_prefix = output_prefix + "_qc"
    }

    # Step 3: Convert to MSAV format using Minimac4
    call MinimacConversionTask.MinimacConversion {
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
        description: "Federated Genotype Imputation Pipeline - Convert VCF to MSAV format"
        author: "GA4GH Hackathon 2025 - African Genomics Team"
        email: "team@afrigenomics.org"
        version: "1.0.0"
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