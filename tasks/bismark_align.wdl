version 1.0

task bismark_align {

    input {
        File reads_1
        File ref_genome
        File ref_genome_index_tar
        String library_type
        File? reads_2
        String? bismark_opt_args

        String? docker_im
        Int? disk_sp
        Int? cores
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
        ref_genome:{
            help: "Reference genome.",
            patterns: ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz"]
        }
        ref_genome_index_tar:{
            help: "Reference genome index (tarball)."
        }
        bismark_opt_args:{
            help: "Additional bismark options (see https://github.com/FelixKrueger/Bismark/tree/master/Docs Appendix section): If no additional options are specified Bismark will use a set of default values."
        }
        library_type:{
            help: "The library type (directional; non_directional)."
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

    meta {
        title: "Bismark Bisulfite Mapper"
        description: "Map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step."
    }

    String sample_name = basename(basename(basename(basename(reads_1, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")

    Int disk_space = select_first([disk_sp, 6 * ceil(size(reads_1, "G") + size(reads_2, "G")) + 3 * ceil(size(ref_genome_index_tar, "G"))])
    Int multi_core = select_first([cores, 16])
    Int mem = if (multi_core * 4) < 108 then 108 else (multi_core * 4)
    
    command <<<
        set -exo pipefail

        # CPU cores for alignment need at least 18 GB RAM per core
        # https://github.com/FelixKrueger/Bismark/issues/111
        nproc_alignment=$(expr $(free -g | grep Mem: | awk '{print $2}') / 18)

        if [ "${nproc_alignment}" -lt 3 ];then
            echo "WARNING: CPU for alignment was 0, setting to 3..."
            nproc_alignment=3 # minimum for alignment is 2 (set to use at least 3)
        fi

        # generate key bismark arguments
        if [ "~{library_type}" = "non_directional" ];then 
            ARGS='--non_directional'
        fi

        if [[ -n $ARGS ]]; then 
            if [[ -n "~{bismark_opt_args}" ]];then
                ARGS+=" ~{bismark_opt_args}"
            fi
        else
            if [[ -n "~{bismark_opt_args}" ]];then
                ARGS="~{bismark_opt_args}"
            fi
        fi
        
        # create directory for reference genome
        mkdir -p $HOME/ref_genome
        mv "~{ref_genome}" $HOME/ref_genome/
        
        # 1) Untar index to $HOME/ref_genome
        tar -xf ~{ref_genome_index_tar} -C $HOME/ref_genome/
        
        ## Specifying --basename in conjuction with --multicore is currently not supported (but we are aiming to fix this soon).
        # 2) Bismark read alignment and methylation calling: PE versus SE reads
        if [[ -n "~{reads_2}" ]]; then
            bismark --genome $HOME/ref_genome/ \
            --gzip \
            $ARGS \
            --parallel $nproc_alignment \
            -1 ~{reads_1} \
            -2 ~{reads_2}

            BAM_ALIGNED="~{sample_name}"_bismark_bt2_pe.bam
            ALIGNMENT_REPORT="~{sample_name}"_bismark_bt2_PE_report.txt
        else    
            bismark --genome $HOME/ref_genome/ \
            --gzip \
            $ARGS \
            --parallel $nproc_alignment \
            ~{reads_1}

            BAM_ALIGNED="~{sample_name}"_bismark_bt2.bam
            ALIGNMENT_REPORT="~{sample_name}"_bismark_bt2_SE_report.txt
        fi

        mv $BAM_ALIGNED "~{sample_name}".bam
        mv $ALIGNMENT_REPORT "~{sample_name}"_bismark_alignment_report.txt

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GG6B1280pfp6ZxfkJVKKzgF9"])
        cpu: multi_core
        dx_timeout: "168H"
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File bismark_alignment_report = "${sample_name}_bismark_alignment_report.txt"
        File mapped_reads = "${sample_name}.bam"
    }
}