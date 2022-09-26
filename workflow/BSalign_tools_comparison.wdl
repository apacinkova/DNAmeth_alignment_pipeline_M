version 1.0

import "fastqc.wdl" as fastqc
import "alignment_statistics.wdl" as alignment_statistics

workflow BSalign_tools_comparison {
    input {
        File bismark_bam
        File bwa_meth_bam
        File ref_genome

        String? docker_image
        Int? disk_sp
    }

    parameter_meta {
        bismark_bam:{
            help: "Aligned reads from Bismark."
        }
        bwa_meth_bam:{
            help: "Aligned reads from BWA-meth."
        }
        ref_genome:{
            help: "Reference genome.",
            patterns: ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz"]
        }
        docker_image: {
            help: "Docker image tarball containing execution environment. Please pass as a string in format dx://project-xxx:file-yyy."
        }
    }

    meta {
        title: "BS reads alignment evaluation"
        description: "Post alignment QC and comparison of Bismark and BWA-meth."
    }

    String prefix_1 = basename(bismark_bam, ".bam")
    String prefix_2 = basename(bwa_meth_bam, ".bam")

    # FastQC
    call fastqc.fastqc as fastqc_bismark {
    input:
        reads_1 = bismark_bam,
        filename_1 = prefix_1,
        docker_im = docker_image
    }

    # FastQC
    call fastqc.fastqc as fastqc_bwa_meth {
    input:
        reads_1 = bwa_meth_bam,
        filename_1 = prefix_2,
        docker_im = docker_image
    }
    
    call alignment_statistics.alignment_statistics as alignment_statistics {
    input: 
        bismark_bam = bismark_bam,
        bwa_meth_bam = bwa_meth_bam,
        ref_genome = ref_genome,
        docker_im = docker_image
    }

    output {
        File fastQC_bismark = fastqc_bismark.fastqc_res_1
        File fastQC_bwa_meth = fastqc_bwa_meth.fastqc_res_1

        File alignment_stat_res = alignment_statistics.alignment_stat_res
        File bismark_SummaryMetric = alignment_statistics.bismark_SummaryMetric
        File bwa_meth_SummaryMetric = alignment_statistics.bwa_meth_SummaryMetric
    }
}