version 1.0

task fastp {

    input {
        File reads_1
        File? reads_2
        String prefix_1
        String? prefix_2
        String? docker_im
        Int? disk_sp
        Int? actual_cpu
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
        prefix_1: {
            help: "File name fastq reads 1 without the suffix."
        }
        prefix_2: {
            help: "File name fastq reads 2 without the suffix."
        }
        docker_im: {
            help: "Docker image containing runtime environment."
        }
        disk_sp: {
            help: "Required disk space (in GB) to run the app."
        }
        actual_cpu: {
            help: "Required cpu to run the app."
        }
        actual_mem: {
            help: "Required memory (in GB) to run the app."
        }
    }

    meta {
        title: "fastp"
        description: "A tool for fast quality control of FastQ files (trimming and filtering is disabled)."
    }


    Int actual_disk_sp = select_first([
    disk_sp, 
    ceil(2 * ((size(reads_1, "G") + (size(reads_2, "G")))))
    ])

    command <<<
        set -e +x -o pipefail

        # check if reads_2 are defined to run single or paired-end mode
        if [[ -n "~{reads_2}" ]]; then
            fastp -i ~{reads_1} -I ~{reads_2} \
            --overrepresentation_analysis \
            -o "~{prefix_1}".1.fastq.gz -O "~{prefix_2}".2.fastq.gz \
            --thread $(nproc)\
            --dont_overwrite \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --json "~{prefix_1}".json \
             --html "~{prefix_1}".fastp.html
        else
            fastp -i ~{reads_1} \
            --overrepresentation_analysis \
            -o "~{prefix_1}".1.fastq.gz \
            --thread $(nproc)\
            --dont_overwrite \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --json "~{prefix_1}".json \
            --html "~{prefix_1}".fastp.html
        fi
    >>>

    runtime{
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GFVxG0j0pfpGxKfFB3BBj0X9"])
        cpu: select_first([actual_cpu, 1])
        memory: select_first([actual_mem, 4]) + "GB"
        disks: "local-disk ${actual_disk_sp} SSD"
    }
    
    output {
        File fastp_html = "${prefix_1}.fastp.html"
        File fastp_json = "${prefix_1}.json"
    } 
}