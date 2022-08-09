version 1.0

task bismark {

    input {
        File reads_1
        String sample_name
        File? reads_2
        File? ref_genome
        File? ref_genome_index_tar
        Boolean wgbs = true
        String? bismark_optional_args
        String? docker_im
        Int? disk_sp
        Int? cores
    }

    parameter_meta {
        reads_1:{
            help: "Fastq reads 1 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        sample_name:{
            help: "Sample name - needs to be the whole file name without suffix."
        }
        reads_2:{
            help: "Fastq reads 2 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        ref_genome:{
            help: "Reference genome."
        }
        ref_genome_index_tar:{
            help: "Reference genome index (tarball)."
        }
        wgbs:{
            help: "Sequencing strategies: whole genome shotgun BS-Seq (WGBS)."
        }
        bismark_optional_args:{
            help: "Additional bismark options (see https://github.com/FelixKrueger/Bismark/tree/master/Docs Appendix section): If no additional options are specified Bismark will use a set of default values."
        }
        docker_im: {
            help: "Docker image containing runtime environment."
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

    Int disk_space = if defined(ref_genome_index_tar) then select_first([disk_sp, ceil(2 * (size(reads_1, "G") + size(reads_2, "G")) + size(ref_genome_index_tar, "G"))]) else select_first([disk_sp, ceil(2 * ((size(reads_1, "G") + (size(reads_2, "G"))))) + 16])
    Int multi_core = select_first([cores, 10])
    Int mem = multi_core * 24
    String ref_genome_name = if defined(ref_genome) then basename(basename(basename(ref_genome, ".gz"), ".fasta"), ".fa") else "reference_genome"
    
    command <<<
        set -e +x -o pipefail
        # create directory for reference genome
        mkdir -p $HOME/ref_genome
        
        # 1) Run Bismark genome preparation
        if [[ -n "~{ref_genome}" ]]; then
            mkdir -p $HOME/ref_genome_tar
            mv ~{ref_genome} $HOME/ref_genome
            bismark_genome_preparation $HOME/ref_genome/
            tar -czvf $HOME/ref_genome_tar/"~{ref_genome_name}".tar.gz $HOME/ref_genome/
            mv $HOME/ref_genome_tar/"~{ref_genome_name}".tar.gz .
        else
            # Untar index to $HOME/ref_genome
            tar -xf ~{ref_genome_index_tar} -C $HOME/ref_genome/ --strip-components 2
        fi

        # 2) Bismark read alignment and methylation calling
        ## Specifying --basename in conjuction with --multicore is currently not supported (but we are aiming to fix this soon).
        if [[ -n "~{reads_2}" ]]; then
            bismark --genome $HOME/ref_genome/ \
            #--basename "~{sample_name}" \
            --gzip \
            ~{bismark_optional_args} \
            --parallel $(nproc) \
            -1 ~{reads_1} \
            -2 ~{reads_2}  > "~{sample_name}".bismark.out
            mv "~{sample_name}"_PE_report.txt "~{sample_name}"_bismark_alignment_report.txt
        else
            bismark --genome $HOME/ref_genome/ \
            #--basename "~{sample_name}" \
            --gzip \
            ~{bismark_optional_args} \
            --parallel $(nproc) \
            ~{reads_1}  > "~{sample_name}".bismark.out
             mv "~{sample_name}"_bismark_SE_report.txt "~{sample_name}"_bismark_alignment_report.txt
        fi

        if [ ~{wgbs} = true ];then
            deduplicate_bismark --bam "~{sample_name}".bam

            # 3) Bismark methylation extractor
            bismark_methylation_extractor --gzip --bedGraph "~{sample_name}".deduplicated.bam

            # 4) Generate a graphical HTML overall summary
            bismark2summary --output "~{sample_name}".bismark_summary_report.html \
            --alignment_report "~{sample_name}"_bismark_alignment_report.txt 
        else
            # 3) Bismark methylation extractor
            bismark_methylation_extractor --gzip --bedGraph "~{sample_name}".bam

            # 4) Generate a graphical HTML overall summary
            bismark2summary --output "~{sample_name}".bismark_summary_report.html \
            --alignment_report "~{sample_name}"_bismark_alignment_report.txt 
        fi
    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GFbvjp00pfp8JZv52Zv5xxvP"])
        cpu: multi_core
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File? ref_gen_index = "${ref_genome_name}.tar.gz"
        File bismark_map_out = "${sample_name}.bismark.out"
        File bismark_alignment_report = "${sample_name}_bismark_alignment_report.txt"
        File mapped_reads = "${sample_name}.bam"
        File? mapped_reads_deduplicated = "${sample_name}.deduplicated.bam"
        File bismark_summary_report = "${sample_name}.bismark_summary_report.html"
        File bedgraph = "bedgraph.html"
        File meth_CpG_context = "CpG_context_${sample_name}_bismark_bt2.txt.gz"
        File meth_CHG_context = "CHG_context_${sample_name}_bismark_bt2.txt.gz"
        File meth_CHH_context = "CHH_context_${sample_name}_bismark_bt2.txt.gz"
    }
}