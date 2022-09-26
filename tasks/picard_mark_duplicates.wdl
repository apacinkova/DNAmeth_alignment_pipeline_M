version 1.0

task picard_mark_duplicates {
    input {
        File mapped_reads_sorted_bam
        Boolean patterned_flowcell_model = true

        String? docker_im
        Int? disk_sp
        Int? cores
    }

    meta {
        description: "Tool to mark duplicates in coordinate-sorted bam file"
    }

    parameter_meta {
        mapped_reads_sorted_bam: {
            help: "Input coordinate-sorted bam with aligned reads",
            patterns: ["*.bam"]
        }
        patterned_flowcell_model: {
            help: "A) Patterned flow cells have clusters with defined sizes, defined shapes, and ordered spacing. B) Nonpatterned flow cells have clusters with varied sizes, undefined shapes, and irregular spacing."
        }
        docker_im: {
            help: "Docker image tarball containing execution environment. Please pass as a string in format dx://project-xxx:file-yyy."
        }
        disk_sp: {
            help: "Required disk space (in GB) to run the app."
        }
        cores: {
            help: "Required number of CPU cores to be used for mapping."
        }
    }

    String sample_name = basename(mapped_reads_sorted_bam, "_sorted.bam")
    String output_bam = sample_name + '.deduplicated.bam'
    String output_sorted_bam = sample_name + '.sorted.deduplicated.bam'
    Int optical_dup_pixel_dist = if patterned_flowcell_model then 2500 else 100
    String marked_dup_metrics =  sample_name + '.duplicate_metrics'

    Int disk_space = select_first([disk_sp, 2 * ceil(size(mapped_reads_sorted_bam, "G"))])
    Int multi_core = select_first([cores, 2])
    Int mem = multi_core * 8

    command <<<

        set -uexo pipefail

        # 1) Remove duplicates from BWA-aligned BAM
        java -jar /usr/local/bin/picard.jar MarkDuplicates \
        INPUT=~{mapped_reads_sorted_bam} \
        OUTPUT=~{output_bam} \
        METRICS_FILE=~{marked_dup_metrics} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        CREATE_INDEX=false \
        READ_NAME_REGEX='[a-zA-Z0-9\-\_]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*' \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{optical_dup_pixel_dist} \
        VALIDATION_STRINGENCY=SILENT \
        PROGRAM_RECORD_ID=null

        # 2) Sort BAM output by leftmost coordinates
        samtools sort -O BAM -@ $(nproc) ~{output_bam} > ~{output_sorted_bam}
        
        # 3) Index coordinate-sorted BAM file for fast random access
        samtools index -@ $(nproc) ~{output_sorted_bam}
        
        # 4) Count the number of alignments for each FLAG type after deduplication
        samtools flagstat ~{output_sorted_bam} > ~{sample_name}_bwa_meth_sorted_alignment_dedup_report.txt

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GGFPf1j0pfpFKk7480z3BKPg"])
        cpu: multi_core
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File dedup_aligned_reads_sorted_bam = "${output_sorted_bam}"
        File dedup_aligned_reads_sorted_bam_index = "${output_sorted_bam}.bai"
        File deduplicate_metrics = "${marked_dup_metrics}"
        File bwa_meth_alignment_report_dedup = "${sample_name}_bwa_meth_sorted_alignment_dedup_report.txt"
    }
}
