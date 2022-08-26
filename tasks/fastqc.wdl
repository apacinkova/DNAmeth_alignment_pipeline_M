version 1.0

task fastqc {
    input {
        File reads_1
        String? filename_1
        File? reads_2
        String? filename_2
        String? docker_im
        Int? disk_sp
        Int? actual_mem
    }

    parameter_meta {
        reads_1:{
            help: "Fastq reads 1 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        filename_1:{
            help: "Reads 1 file name - needs to be the whole file name without suffix."
        }
        reads_2:{
            help: "Fastq reads 2 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        filename_2:{
            help: "Reads 2 file name - needs to be the whole file name without suffix."
        }
        docker_im: {
            help: "Docker image tarball containing execution environment. Please pass as a string in format dx://project-xxx:file-yyy."
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

    String filename = select_first([filename_1,basename(basename(basename(basename(reads_1, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")])
    
    command <<<
        set -e +x -o pipefail
        mkdir $HOME/fastQC_reports

        # Run FastQC using either paired-end or single-end mode
        if [[ -n "~{reads_2}" ]]; then
            fastqc \
            ~{reads_1} \
            ~{reads_2} \
            --outdir $HOME/fastQC_reports

            mv $HOME/fastQC_reports/~{filename_2}_fastqc.html ~{filename_2}_fastqc_report.html
            
        else
            fastqc \
            ~{reads_1} \
            --outdir $HOME/fastQC_reports

        fi
        mv $HOME/fastQC_reports/~{filename}_fastqc.html ~{filename}_fastqc_report.html

    >>>

   runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GG1j3x80pfp2QxJy4q7q2PQ5"])
        cpu: 1
        memory: select_first([actual_mem, 8]) + "GB"
        disks: "local-disk ${disk_spc} SSD"
    }

   output {
        File fastqc_res_1 = "${filename}_fastqc_report.html"
        File? fastqc_res_2 = "${filename_2}_fastqc_report.html"
    }
}