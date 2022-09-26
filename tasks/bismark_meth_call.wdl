version 1.0

task bismark_meth_call {

    input {
        File mapped_reads_bam
        Boolean wgbs = true
        File bismark_alignment_report
        String? bismark_opt_args

        String? docker_im
        Int? disk_sp
        Int? cores
    }

    parameter_meta {
        mapped_reads_bam:{
            help: "Mapped reads using Bismark aligner.",
            patterns: ["*.bam"]
        }
        wgbs:{
            help: "Sequencing strategy: whole genome shotgun BS-Seq (WGBS)."
        }
        bismark_alignment_report:{
            help: "Bismark alignment report."
        }
        bismark_opt_args:{
            help: "Additional bismark options (see https://github.com/FelixKrueger/Bismark/tree/master/Docs Appendix section): If no additional options are specified Bismark will use a set of default values."
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

    String sample_name = basename(mapped_reads_bam, ".bam,")

    Int disk_space = select_first([disk_sp, 3 * ceil(size(mapped_reads_bam, "G"))])
    Int multi_core = select_first([cores, 16])
    Int mem = if (multi_core * 4) < 108 then 108 else (multi_core * 4)
    
    command <<<
        set -exo pipefail
        WORKDIR=$PWD

        # CPU cores for methylation calling need at least 18 GB RAM per core
        # https://github.com/FelixKrueger/Bismark/issues/111
        nproc_meth_call=$(expr $(free -g | grep Mem: | awk '{print $2}') / 18)

        if [ "${nproc_meth_call}" -lt 3 ];then
            echo "WARNING: CPU for methylation calling was 0, setting to 3..."
            nproc_meth_call=3 # minimum for methylation calling is 2 (set to use at least 3)
        fi

        mv ~{bismark_alignment_report} $WORKDIR/~{sample_name}_SE_report.txt
        BAM_ALIGNED=~{mapped_reads_bam}
        INPUT_BAM_ME=~{mapped_reads_bam}
        ARGS="~{bismark_opt_args}"

        # 1) Remove alignments to the same position in the genome which can arise by e.g. PCR amplification (WGBS or PBAT)
        if [ "~{wgbs}" = "true" ] || [[ $ARGS = *"pbat"* ]];then
            deduplicate_bismark --bam $BAM_ALIGNED
            mv $(ls *deduplicated.bam) "~{sample_name}".deduplicated.bam
            mv $(ls *deduplication_report.txt) "~{sample_name}.deduplication_report.txt"
            INPUT_BAM_ME="~{sample_name}".deduplicated.bam
        fi

        # 2) Bismark methylation extractor
        mkdir meth_extractor
        bismark_methylation_extractor --output meth_extractor --gzip --bedGraph $INPUT_BAM_ME
        tar -czvf meth_extractor_res.tar.gz meth_extractor        

        # 3) Generate a graphical HTML overall summary using original BAM file
        cp meth_extractor/* ./
        #bismark2summary ~{mapped_reads_bam}

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GG6B1280pfp6ZxfkJVKKzgF9"])
        cpu: multi_core
        dx_timeout: "168H"
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File? mapped_reads_deduplicated = "${sample_name}.deduplicated.bam"
        File? bismark_deduplication_report = "${sample_name}.deduplication_report.txt"
        File methylation_extractor = "meth_extractor_res.tar.gz"
        #File bismark_summary_report = "bismark_summary_report.html"
    }
}