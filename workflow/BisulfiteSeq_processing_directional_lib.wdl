version 1.0

import "fastqc.wdl" as fastqc
import "trim_galore.wdl" as trimGalore
import "bwa_meth_index.wdl" as bwa_meth_index
import "bwa_meth.wdl" as bwa_meth
import "picard_mark_duplicates.wdl" as mark_duplicates
import "methyl_dackel.wdl" as methyl_dackel

workflow BisulfiteSeq_processing_directional_lib {
    input {
        File reads_1
        String lib_type
        File ref_genome
        File? ref_genome_index_tar
        Boolean wgbs = true
        File? reads_2
        Boolean QC_pre_trim = true
        Boolean trim_reads = true
        String? trim_galore_optional_args
        Boolean patterned_flowcell = true
        String? methyl_dackel_optional_args

        String? fastqc_docker_image
        String? trim_galore_docker_image
        String? bwameth_docker_image
    }

    parameter_meta {
        reads_1:{
            help: "Fastq reads 1 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        lib_type:{
            help: "The library type MUST be directional; if you are not sure choose bismark alignment tool with library type detection applet."
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
        QC_pre_trim:{
            help: "Perform quality control before trimming?"
        }
        trim_reads:{
            help: "Perform low-quality and adapter trimming of bisulfite reads?"
        }
        trim_galore_optional_args:{
            help: "Additional TrimGalore! options (see https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md): If no additional options are specified TrimGalore! will use a set of default values."
        }
        patterned_flowcell:{
            help: "A) Patterned flow cells have clusters with defined sizes, defined shapes, and ordered spacing. B) Nonpatterned flow cells have clusters with varied sizes, undefined shapes, and irregular spacing."
        }
        methyl_dackel_optional_args:{
            help: "Additional MethylDackel options (see https://github.com/dpryan79/methyldackel)"
        }
        fastqc_docker_image:{
            help: "Docker image containing runtime environment for FastQC."
        }
        trim_galore_docker_image:{
            help: "Docker image containing runtime environment for TrimGalore!."
        }
        bwameth_docker_image:{
            help: "Docker image containing runtime environment for BWA-meth."
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
    if(QC_pre_trim){
        call fastqc.fastqc as fastqc {
        input:
            reads_1 = reads_1,
            filename_1 = prefix_fq1,
            reads_2 = reads_2,
            filename_2 = prefix_fq2,
            docker_im = fastqc_docker_image
        }
    }
    
    # Reference genome index
    if(!defined(ref_genome_index_tar))
    {
        call bwa_meth_index.bwameth_index as bwa_meth_index {
        input:
            ref_genome = ref_genome,
            docker_im = bwameth_docker_image
        }
    }
    File rg_indx_tar = select_first([ref_genome_index_tar, bwa_meth_index.ref_gen_index])

   # TrimGalore (adapters and low quality reads trimming)
   if(trim_reads){
        call trimGalore.trimGalore as trim_galore {
        input:
            reads_1 = reads_1,
            filename_1 = prefix_fq1,
            reads_2 = reads_2,
            filename_2 = prefix_fq2,
            library_type = lib_type,
            optional_args = trim_galore_optional_args,
            docker_im = trim_galore_docker_image
        }
    }
    File reads_to_align_1 = select_first([trim_galore.trimmed_reads_1, reads_1])
    if(defined(reads_2))
    {
        File reads_to_align_2 = select_first([trim_galore.trimmed_reads_2, reads_2])
    }

    # BWA-meth alignment
    call bwa_meth.bwa_meth as bwa_meth {   
        input:
            reads_1 = reads_to_align_1,
            library_type = lib_type,
            reads_2 = reads_to_align_2,
            ref_genome_index_tar = rg_indx_tar,
            docker_im = bwameth_docker_image
    }

    if(wgbs)
    {
        # Picard MarkDuplicates
        call mark_duplicates.picard_mark_duplicates as mark_duplicates {   
            input:
                mapped_reads_sorted_bam = bwa_meth.mapped_reads_sorted_bam,
                patterned_flowcell_model = patterned_flowcell,
                docker_im = bwameth_docker_image
        }
    }
    File aligned_sorted_reads = select_first([mark_duplicates.dedup_aligned_reads_sorted_bam,bwa_meth.mapped_reads_sorted_bam])
    File aligned_sorted_reads_index = select_first([mark_duplicates.dedup_aligned_reads_sorted_bam_index,bwa_meth.mapped_reads_sorted_bam_index])
        
    # MethylDackel extract methylation metrics
    call methyl_dackel.MethylDackel as methyl_dackel {   
        input:
            mapped_reads_sorted_bam = aligned_sorted_reads,
            mapped_reads_sorted_bai = aligned_sorted_reads_index,
            ref_genome = ref_genome,
            methyl_dackel_opt_args = methyl_dackel_optional_args,
            docker_im = bwameth_docker_image
    }

    output {
        File? fastQC_pre_trim_reads_1 = fastqc.fastqc_res_1
        File? fastQC_pre_trim_reads_2 = fastqc.fastqc_res_2

        File? reference_genome_index = bwa_meth_index.ref_gen_index

        File? trimmed_reads_1 = trim_galore.trimmed_reads_1
        File? trimmed_reads_2 = trim_galore.trimmed_reads_2
        Array[File]? trimming_report = trim_galore.trimming_report
        File? fastQC_post_trim_1 = trim_galore.fastqc_post_res_1
        File? fastQC_post_trim_2 = trim_galore.fastqc_post_res_2

        File bwa_meth_alignment_report = bwa_meth.bwa_meth_alignment_report
        File bwa_meth_aligned_reads_sorted_bam = bwa_meth.mapped_reads_sorted_bam
        File bwa_meth_aligned_reads_sorted_bai = bwa_meth.mapped_reads_sorted_bam_index
        
        File? bwa_meth_dedup_aligned_reads_sorted_bam = mark_duplicates.dedup_aligned_reads_sorted_bam
        File? bwa_meth_dedup_aligned_reads_sorted_bai = mark_duplicates.dedup_aligned_reads_sorted_bam_index
        File? MarkDuplicates_deduplicate_metrics = mark_duplicates.deduplicate_metrics
        File? MarkDuplicates_deduplication_report = mark_duplicates.bwa_meth_alignment_report_dedup

        File methylation_bedgraph = methyl_dackel.out_bedgraph
    }
}