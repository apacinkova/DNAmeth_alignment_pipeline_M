version 1.0

task bwa_meth {

    input {
        File reads_1
        String? sample_name
        File ref_genome_index_tar
        String library_type
        File? reads_2

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
            help: "Sample name."
        }
        ref_genome_index_tar:{
            help: "Reference genome index (tarball).",
            patterns: ["*.tar.gz"]
        }
        library_type:{
            help: "The library type (directional; non_directional)."
        }
        reads_2:{
            help: "Fastq reads 2 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
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
        title: "BWA-METH Bisulfite Mapper"
        description: "Map bisulfite treated sequencing reads (directional library) to a genome of interest."
    }

    String sample_nam = select_first([sample_name,basename(basename(basename(basename(reads_1, ".fastq.gz"), ".fastq"), ".fq"), ".fq.gz")])
    String ref_index_prefix = basename(ref_genome_index_tar, "_index_bwa.tar.gz")

    Int disk_space = select_first([disk_sp, ceil(4 * (size(reads_1, "G") + size(reads_2, "G")) + 4 * size(ref_genome_index_tar, "G"))])
    Int multi_core = select_first([cores, 32])
    Int mem = multi_core * 3
    
    command <<<
        set -exo pipefail

        # library type check
        if [[ "~{library_type}" != "directional" ]]; then
            echo "The library must be directional! For non-directional library use Bismark."
            exit
        fi
        
        # CPU cores for samtools sort need at least 5 GB RAM per core (use the same for alignment)
        # https://www.biostars.org/p/18933/
        nproc_alignment=$(expr $(free -g | grep Mem: | awk '{print $2}') / 5)

        if [ "${nproc_alignment}" -lt 3 ];then
            echo "WARNING: CPU for alignment was 0, setting to 3..."
            nproc_alignment=3 # set minimum as 3
        fi

        # 1) Untar reference genome index
        tar -xf ~{ref_genome_index_tar}

        # 2) BWA-METH read alignment, sort alignments by leftmost coordinates and output as BAM
        bwameth.py --threads $nproc_alignment \
        --reference index_bwa/~{ref_index_prefix} \
        ~{reads_1} ~{reads_2} > ~{sample_nam}_bwa_meth.sam

        # 3) Sort output SAM by leftmost coordinates and save it as BAM
        samtools sort -O BAM -@ $nproc_alignment ~{sample_nam}_bwa_meth.sam > ~{sample_nam}_bwa_meth_sorted.bam
        
        # 4) Index coordinate-sorted BAM file for fast random access
        samtools index -@ $nproc_alignment ~{sample_nam}_bwa_meth_sorted.bam 

        # 5) Count the number of alignments for each FLAG type
        samtools flagstat ~{sample_nam}_bwa_meth_sorted.bam > ~{sample_nam}_bwa_meth_sorted_alignment_report.txt

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GGFPf1j0pfpFKk7480z3BKPg"])
        cpu: multi_core
        dx_timeout: "96H"
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File mapped_reads_sorted_bam = "${sample_nam}_bwa_meth_sorted.bam"
        File mapped_reads_sorted_bam_index = "${sample_nam}_bwa_meth_sorted.bam.bai"
        File bwa_meth_alignment_report = "${sample_nam}_bwa_meth_sorted_alignment_report.txt"
    }
}