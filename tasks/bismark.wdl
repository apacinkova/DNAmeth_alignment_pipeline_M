version 1.0

task bismark {

    input {
        File reads_1
        File ref_genome
        File ref_genome_index_tar
        Boolean wgbs = true
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
        wgbs:{
            help: "Sequencing strategy: whole genome shotgun BS-Seq (WGBS)."
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

    Int disk_space = select_first([disk_sp, ceil(2 * (size(reads_1, "G") + size(reads_2, "G")) + 4 * size(ref_genome_index_tar, "G"))])
    Int multi_core = select_first([cores, 16])
    Int mem = if (multi_core * 4) < 54 then 54 else (multi_core * 4)
    
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

            mv "~{sample_name}"_bismark_bt2_pe.bam "~{sample_name}".bam
            mv "~{sample_name}"_bismark_bt2_PE_report.txt "~{sample_name}"_bismark_alignment_report.txt
            EXTRACTOR_ARGS=--paired-end
        else    
            bismark --genome $HOME/ref_genome/ \
            --gzip \
            $ARGS \
            --parallel $nproc_alignment \
            ~{reads_1}

            mv "~{sample_name}"_bismark_bt2.bam "~{sample_name}".bam
            mv "~{sample_name}"_bismark_bt2_SE_report.txt "~{sample_name}"_bismark_alignment_report.txt
            EXTRACTOR_ARGS=--single-end
        fi

        INPUT_BAM="~{sample_name}".bam

        # 3) Remove alignments to the same position in the genome which can arise by e.g. PCR amplification (WGBS or PBAT)
        if [ ~{wgbs} = true ] || [[ $ARGS = *"pbat"* ]];then
            deduplicate_bismark --bam $INPUT_BAM
            INPUT_BAM="~{sample_name}".deduplicated.bam
        fi

        # 4) Bismark methylation extractor
        mkdir meth_extractor
        bismark_methylation_extractor $EXTRACTOR_ARGS --output meth_extractor --gzip --bedGraph $INPUT_BAM
        tar -czvf meth_extractor_res.tar.gz meth_extractor        

        # 5) Generate a graphical HTML overall summary
        # choose the .bam file that is NOT deduplicated
        INPUT_BAM_NOTDUP=$(find . -type f -name "*.bam" | grep -vw "deduplicated" | awk '{print substr( $0, 3 )}')
        cp meth_extractor/* ./
        bismark2summary $INPUT_BAM_NOTDUP

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GG4GPQQ0pfpBf22y8gzP9Z9k"])
        cpu: multi_core
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File bismark_alignment_report = "${sample_name}_bismark_alignment_report.txt"
        File mapped_reads = "${sample_name}.bam"
        File? mapped_reads_deduplicated = "${sample_name}.deduplicated.bam"
        File? bismark_deduplication_report = "${sample_name}.deduplication_report.txt"
        File methylation_extractor = "meth_extractor_res.tar.gz"
        File bismark_summary_report = "bismark_summary_report.html"
    }
}