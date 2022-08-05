version 1.0

task trimGalore {

    input {
        File reads_1
        String sample_name
        File? reads_2
        String? docker_im
        Int? disk_sp
        Int? cores
        Int? actual_mem
    }

    Int disk_space = select_first([
        disk_sp, 
        ceil(2 * ((size(reads_1, "G") + (size(reads_2, "G")))))
    ])

    parameter_meta {
        reads_1:{
            help: "Fastq reads 1 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        reads_2:{
            help: "Fastq reads 2 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        docker_im: {
            help: "Docker image containing runtime environment."
        }
        disk_sp: {
            help: "Required disk space (in GB) to run the app."
        }
        cores: {
            help: "Required number of CPU cores to be used for trimming. Default values: single-end --cores=9; paired-end --cores=15."
        }
        actual_mem: {
            help: "Required memory (in GB) to run the app."
        }
    }

    meta {
        title: "TrimGalore!"
        description: "Taking appropriate QC measures for RRBS-type or other -Seq applications with Trim Galore!"
    }

    command <<<
        set -e +x -o pipefail

        # create directory for outputs
        mkdir $HOME/fastQC_reports
        mkdir $HOME/trimmed_data

        # Run TrimGalore! using either paired-end or single-end mode
        if [[ -n "~{reads_2}" ]]; then
            trim_galore \
            --basename "~{sample_name}" \
            --phred33 \
            --fastqc \
            --fastqc_args "--extract --outdir $HOME/fastQC_reports/" \
            --paired \
            --cores $(nproc) \
            --output_dir $HOME/trimmed_data/ \
            --gzip \
            "~{reads_1}" \
            "~{reads_2}"

            mv $HOME/trimmed_data/"~{sample_name}"_val_1.fq.gz "~{sample_name}"_val_1.fq.gz
            mv $HOME/trimmed_data/"~{sample_name}"_val_2.fq.gz "~{sample_name}"_val_2.fq.gz

        else
            trim_galore \
            --basename "~{sample_name}" \
            --phred33 \
            --fastqc \
            --fastqc_args "--extract --outdir $HOME/fastQC_reports/" \
            --cores $(nproc) \
            --output_dir $HOME/trimmed_data/ \
            --gzip \
            "~{reads_1}"

            mv $HOME/trimmed_data/"~{sample_name}"_trimmed.fq.gz "~{sample_name}"_trimmed.fq.gz

        fi

        mv $HOME/fastQC_reports/*/fastqc_report.html fastqc_report_after_trim.html
        mv $HOME/trimmed_data/*_trimming_report.txt trimming_report.txt

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GFb34vj0pfp28G6y0Jzqp1FZ"])
        cpu: select_first([cores,if defined(reads_2) then 14 else 8])
        memory: select_first([actual_mem,if defined(reads_2) then 56 else 32]) + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File? trimmed_reads_1 = "${sample_name}_val_1.fq.gz"
        File? trimmed_reads_2 = "${sample_name}_val_2.fq.gz"
        File? trimmed_reads = "${sample_name}_trimmed.fq.gz"
        File trimming_report = "trimming_report.txt"
        File fastqc_post_res = "fastqc_report_after_trim.html"
    }
}