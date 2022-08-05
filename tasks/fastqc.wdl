version 1.0

task fastqc {
    input {
        File reads_1
        File? reads_2
        String? docker_im
        Int? disk_sp
        Int? actual_mem
    }

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
        actual_mem: {
            help: "Required memory (in GB) to run the app."
        }
    }

    Int disk_spc = select_first([
        disk_sp, 
        ceil(2 * ((size(reads_1, "G") + (size(reads_2, "G")))))
    ]) 

    command <<<
        set -e +x -o pipefail
        mkdir $HOME/fastQC_reports

        # Run FastQC using either paired-end or single-end mode
        if [[ -n "~{reads_2}" ]]; then
            fastqc \
            "~{reads_1}" \
            "~{reads_2}" \
            --extract \
            --outdir $HOME/fastQC_reports
        else
            fastqc \
            "~{reads_1}" \
            --extract \
            --outdir $HOME/fastQC_reports
        fi

        mv $HOME/fastQC_reports/*/fastqc_report.html fastqc_report.html
    >>>

   runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GFb34vj0pfp28G6y0Jzqp1FZ"])
        cpu: 1
        memory: select_first([actual_mem, 8]) + "GB"
        disks: "local-disk ${disk_spc} SSD"
    }

   output {
        File fastqc_res = "fastqc_report.html"
    }
}