version 1.0

import "fastqc.wdl" as fastqc
import "bismark_index.wdl" as bismark_index
import "library_type.wdl" as library_type
import "trim_galore.wdl" as trimGalore
import "bismark.wdl" as bismark

workflow BisulfiteSeq_processing {
    input {
        File reads_1
        String lib_type = "unknown"
        File ref_genome
        File? ref_genome_index_tar
        Boolean wgbs = true
        File? reads_2
        Int? num_reads_align
        String? trim_galore_optional_args
        String? bismark_optional_args

        String? fastqc_docker_image
        String? trim_galore_docker_image
        String? bismark_docker_image
    }

    parameter_meta {
        reads_1:{
            help: "Fastq reads 1 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        lib_type:{
            help: "The library type (directional; non_directional; unknown)."
        }
        ref_genome:{
            help: "Reference genome.",
            patterns: ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz"]
        }
        ref_genome_index_tar:{
            help: "Reference genome index (tarball)."
        }
        wgbs:{
            help: "Sequencing strategy: whole genome shotgun BS-Seq (WGBS)."
        }
        reads_2:{
            help: "Fastq reads 2 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        num_reads_align:{
            help: "Number of reads aligned to the reference genome to detect library type."
        }
        trim_galore_optional_args:{
            help: "Additional TrimGalore! options (see https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md): If no additional options are specified TrimGalore! will use a set of default values."
        }
        bismark_optional_args:{
            help: "Additional Bismark options (see https://github.com/FelixKrueger/Bismark/tree/master/Docs Appendix section): If no additional options are specified Bismark will use a set of default values."
        }
        fastqc_docker_image:{
            help: "Docker image containing runtime environment for FastQC."
        }
        trim_galore_docker_image:{
            help: "Docker image containing runtime environment for TrimGalore!."
        }
        bismark_docker_image:{
            help: "Docker image containing runtime environment for Bismark."
        }
    }

    meta {
        title: "Bisulfite Sequencing Data Processing (QC, trimming, alignment, methylation calling)"
        description: "Trim and align bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single workflow."
    }

    String prefix_fq1 = basename(basename(basename(basename(reads_1, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")
    if(defined(reads_2))
    {
        String prefix_fq2 = basename(basename(basename(basename(reads_2, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")
    }

    # FastQC before trimming
    call fastqc.fastqc as fastqc {
    input:
        reads_1 = reads_1,
        filename_1 = prefix_fq1,
        reads_2 = reads_2,
        filename_2 = prefix_fq2,
        docker_im = fastqc_docker_image
    }

    # Reference genome index
    if(!defined(ref_genome_index_tar))
    {
        call bismark_index.bismark_index as bismark_index {
        input:
            ref_genome = ref_genome,
            docker_im = bismark_docker_image
        }
    }
    File rg_indx_tar = select_first([ref_genome_index_tar, bismark_index.ref_gen_index])

    # Library type detection
    if(lib_type == "unknown")
    {
        call library_type.library_type as library_type {
        input:
            reads_1 = reads_1,
            ref_genome = ref_genome,
            ref_genome_index_tar = rg_indx_tar,
            num_reads_align = num_reads_align,
            docker_im = bismark_docker_image
        } 
    }
    String lib_typ = select_first([library_type.library_type, lib_type])

   # TrimGalore (adapters and low quality reads trimming)
    call trimGalore.trimGalore as trim_galore {
    input:
        reads_1 = reads_1,
        filename_1 = prefix_fq1,
        reads_2 = reads_2,
        filename_2 = prefix_fq2,
        library_type = lib_typ,
        wgbs = wgbs,
        optional_args = trim_galore_optional_args,
        docker_im = trim_galore_docker_image,
    }
    if(defined(trim_galore.trimmed_reads_2))
    {
        File? reads_2_trimmed = trim_galore.trimmed_reads_2
    }

    # Bismark alignment
    call bismark.bismark as bismark {   
        input:
            reads_1 = trim_galore.trimmed_reads_1,
            reads_2 = reads_2_trimmed,
            ref_genome = ref_genome,
            ref_genome_index_tar = rg_indx_tar,
            wgbs = wgbs,
            bismark_opt_args = bismark_optional_args,
            library_type = lib_typ,
            docker_im = bismark_docker_image
    }

    output {
        File? fastQC_pre_trim_reads_1 = fastqc.fastqc_res_1
        File? fastQC_pre_trim_reads_2 = fastqc.fastqc_res_2

        File? reference_genome_index = bismark_index.ref_gen_index

        String library_type_detected = lib_typ
        File? lib_detec_align_report = library_type.bismark_alignment_report

        File trimmed_reads_1 = trim_galore.trimmed_reads_1
        File? trimmed_reads_2 = trim_galore.trimmed_reads_2
        Array[File] trimming_report = trim_galore.trimming_report
        File fastQC_post_trim_1 = trim_galore.fastqc_post_res_1
        File? fastQC_post_trim_2 = trim_galore.fastqc_post_res_2

        File bismark_align_report = bismark.bismark_alignment_report
        File bismark_aligned_reads = bismark.mapped_reads
        File? bismark_aligned_reads_deduplicated = bismark.mapped_reads_deduplicated
        File? bismark_deduplication_report = bismark.bismark_deduplication_report
        File bismark_methylation_extractor = bismark.methylation_extractor
        File bismark_summary_report = bismark.bismark_summary_report
    }
}