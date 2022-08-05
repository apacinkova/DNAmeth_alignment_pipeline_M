version 1.0

import "fastqc.wdl" as fastqc
import "trim_galore.wdl" as trimGalore

workflow bs_align {
    input {
        File reads_1
        File? reads_2
        String? docker_image
        Int? disk_space
        Int? cpu
        Int? mem
    }

    String prefix_fq1 = basename(basename(basename(basename(reads_1, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")
    if(defined(reads_2)){
        String prefix_fq2 = basename(basename(basename(basename(reads_2, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")
    }

    call fastqc.fastqc as fastqc {
        input:
            reads_1 = reads_1,
            reads_2 = reads_2,
            docker_im = docker_image,
            disk_sp = disk_space,
            actual_mem = mem 
    }

    call trimGalore.trimGalore as tg {
        input:
            reads_1 = reads_1,
            reads_2 = reads_2,
            sample_name = prefix_fq1,
            docker_im = docker_image,
            disk_sp = disk_space,
            cores = cpu,
            actual_mem = mem
    }
    
    output {
        File fastQC_pre_trim = fastqc.fastqc_res
        File? trimmed_reads = tg.trimmed_reads
        File? trimmed_reads_1 = tg.trimmed_reads_1
        File? trimmed_reads_2 = tg.trimmed_reads_2
        File trimming_report = tg.trimming_report
        File fastQC_post_trim = tg.fastqc_post_res
    }
}