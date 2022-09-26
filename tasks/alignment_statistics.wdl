version 1.0

task alignment_statistics {

    input {
        File bismark_bam
        File bwa_meth_bam
        File ref_genome

        String? docker_im
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
        docker_im: {
            help: "Docker image tarball containing execution environment. Please pass as a string in format dx://project-xxx:file-yyy."
        }
        disk_sp: {
            help: "Required disk space (in GB) to run the app."
        }
    }

    meta {
        title: "Alignment statistics"
        description: "Using BAM flags to evaluate alignment."
    }

    String bismark_sample_name = basename(bismark_bam, ".bam")
    String bwa_meth_sample_name = basename(bwa_meth_bam, ".bam")

    Int disk_space = select_first([disk_sp, 2*ceil(size(bismark_bam, "G")) + 2*ceil(size(bwa_meth_bam, "G")) + ceil(size(ref_genome, "G"))])
    
    command <<<
        set -exo pipefail
        
        # 1) input BAM files
        echo "Bismark: ~{bismark_bam}" > alignment_statistics.txt
        echo "BWA-meth: ~{bwa_meth_bam}" >> alignment_statistics.txt

        # 2) number of aligned reads in BAM file
        echo "number of aligned reads in bismark BAM:" >> alignment_statistics.txt
        samtools view -c -F 3840 ~{bismark_bam} >> alignment_statistics.txt

        echo "number of aligned reads in bwameth BAM:" >> alignment_statistics.txt
        samtools view -c -F 3840 ~{bwa_meth_bam} >> alignment_statistics.txt

        # 3) number of aligned reads with MAPQ>=10 in BAM file
        echo "number of aligned reads with MAPQ>=10 in bismark BAM:" >> alignment_statistics.txt
        samtools view -h -q 10 -c -F 3840 ~{bismark_bam} >> alignment_statistics.txt

        echo "number of aligned reads with MAPQ>=10 in bwa_meth BAM:" >> alignment_statistics.txt
        samtools view -h -q 10 -c -F 3840 ~{bwa_meth_bam} >> alignment_statistics.txt

        # 4) Sort output Bismark BAM by leftmost coordinates
        samtools sort -O BAM -@ $(nproc) ~{bismark_bam} > ~{bismark_sample_name}_bismark_sorted.bam

        # 5) Picard summary of alignment metrics from a BAM (for Illumina only)
        java -jar /usr/local/bin/picard.jar CollectAlignmentSummaryMetrics \
        I=~{bismark_sample_name}_bismark_sorted.bam \
        BS=true \
        O="~{bismark_sample_name}_bismark_SummaryMetrics.txt" \
        R=~{ref_genome}
                
        java -jar /usr/local/bin/picard.jar CollectAlignmentSummaryMetrics \
        I=~{bwa_meth_bam} \
        BS=true \
        O="~{bwa_meth_sample_name}_SummaryMetrics.txt" \
        R=~{ref_genome}

        # 6) MAPQ plot 
        samtools view ~{bismark_bam} | cut -f 5 | sort | uniq  -c | sort -n | awk '{printf("MAPQ:%s\t%d\n",$2,$1);}' | gnuplot -e " set terminal dumb ;set nokey; plot '-' using 2:xtic(1) with boxes" > ~{bismark_bam}_bismark_MAPQ_plot.txt
        samtools view ~{bwa_meth_bam} | cut -f 5 | sort | uniq  -c | sort -n | awk '{printf("MAPQ:%s\t%d\n",$2,$1);}' | gnuplot -e " set terminal dumb ;set nokey; plot '-' using 2:xtic(1) with boxes" > ~{bwa_meth_bam}_MAPQ_plot.txt

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GGJVZZj0pfpJbzBk1qpkgxF6"])
        cpu: 2
        memory: "32GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File alignment_stat_res = "alignment_statistics.txt"
        File bismark_SummaryMetric = "${bismark_sample_name}_bismark_SummaryMetrics.txt"
        File bwa_meth_SummaryMetric = "${bwa_meth_sample_name}_SummaryMetrics.txt"
        File bismark_MAPQ_plot = "${bismark_bam}_bismark_MAPQ_plot.txt"
        File bwa_meth_MAPQ_plot = "${bwa_meth_bam}_MAPQ_plot.txt"
    }
}