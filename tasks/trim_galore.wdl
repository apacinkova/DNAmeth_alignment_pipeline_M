version 1.0

task trimGalore {

    input {
        File reads_1
        String? filename_1
        String library_type
        Boolean wgbs = true
        File? reads_2
        String? filename_2
        String? optional_args

        String? docker_im
        Int? disk_sp
        Int? cores
        Int? actual_mem
    }

    Int disk_space = select_first([
        disk_sp, 
        ceil(2 * ((size(reads_1, "G") + (size(reads_2, "G")))))
    ])

    String filename = select_first([filename_1,basename(basename(basename(basename(reads_1, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")])

    parameter_meta {
        reads_1:{
            help: "Fastq reads 1 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        filename_1:{
            help: "Reads 1 file name - needs to be the whole file name without suffix."
        }
        library_type:{
            help: "The library type (directional; non_directional)."
        }
        wgbs:{
            help: "Sequencing strategy: whole genome shotgun BS-Seq (WGBS). Note for RRBS using MseI: DNA material digested with MseI (recognition motif: TTAA), use wgbs=true."
        }
        reads_2:{
            help: "Fastq reads 2 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        filename_2:{
            help: "Reads 2 file name - needs to be the whole file name without suffix."
        }
        optional_args:{
            help: "Additional bismark options (see https://github.com/FelixKrueger/Bismark/tree/master/Docs Appendix section): If no additional options are specified Bismark will use a set of default values."
        }
        docker_im: {
            help: "Docker image tarball containing execution environment. Please pass as a string in format dx://project-xxx:file-yyy."
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
        title: "TrimGalore"
        description: "Taking appropriate QC measures for RRBS-type or other -Seq applications with Trim Galore"
    }

    command <<<
        set -e +x -o pipefail

        # generate key TrimGalore! arguments
        if [ ~{wgbs} != true ];then
            if [ "~{library_type}" = "non_directional" ];then 
                ARGS='--non_directional --rrbs'
            else 
                ARGS='--rrbs'
            fi
        fi

        if [[ -n "~{optional_args}" ]];then
            ARGS+=" ~{optional_args}"
        fi

        # create directory for outputs
        mkdir $HOME/fastQC_reports
        mkdir $HOME/trimmed_data

        # Run TrimGalore! using either paired-end or single-end mode
        if [[ -n "~{reads_2}" ]]; then
            trim_galore \
            --phred33 \
            --fastqc \
            --fastqc_args "--outdir $HOME/fastQC_reports/" \
            --paired \
            --cores $(nproc) \
            --output_dir $HOME/trimmed_data/ \
            --gzip \
            $ARGS \
            ~{reads_1} \
            ~{reads_2}

            mv $HOME/trimmed_data/~{filename}_val_1.fq.gz ~{filename}_trimmed.fq.gz
            mv $HOME/trimmed_data/~{filename_2}_val_2.fq.gz ~{filename_2}_trimmed.fq.gz
            mv $HOME/fastQC_reports/~{filename}_val_1_fastqc.html ~{filename}_trimmed_fastqc_report.html
            mv $HOME/fastQC_reports/~{filename_2}_val_2_fastqc.html ~{filename_2}_trimmed_fastqc_report.html

        else
            trim_galore \
            --phred33 \
            --fastqc \
            --fastqc_args "--outdir $HOME/fastQC_reports/" \
            --cores $(nproc) \
            --output_dir $HOME/trimmed_data/ \
            --gzip \
            $ARGS \
            ~{reads_1}

            mv $HOME/trimmed_data/~{filename}_trimmed.fq.gz ~{filename}_trimmed.fq.gz
            mv $HOME/fastQC_reports/~{filename}_trimmed_fastqc.html ~{filename}_trimmed_fastqc_report.html
            
        fi

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GG499g00pfpB3Ybb80X2G5Yg"])
        cpu: select_first([cores,if defined(reads_2) then 14 else 8])
        memory: select_first([actual_mem,if defined(reads_2) then 56 else 32]) + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File trimmed_reads_1 = "${filename}_trimmed.fq.gz"
        File? trimmed_reads_2 = "${filename_2}_trimmed.fq.gz"
        Array[File] trimming_report = glob("*_trimming_report.txt")
        File fastqc_post_res_1 = "${filename}_trimmed_fastqc_report.html"
        File? fastqc_post_res_2 = "${filename_2}_trimmed_fastqc_report.html"
    }
}