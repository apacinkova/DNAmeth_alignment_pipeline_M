version 1.0

task MethylDackel {
    input {
       File mapped_reads_sorted_bam
       File mapped_reads_sorted_bai
       File ref_genome
       String? methyl_dackel_opt_args

       String? docker_im
       Int? disk_sp
       Int? cores
    }

    meta {
        description: "Tool to extract per-base methylation metrics in coordinate-sorted bam file"
    }

    parameter_meta {
        mapped_reads_sorted_bam: {
            help: "Input coordinate-sorted bam with aligned reads",
            patterns: ["*.bam"]
        }
        mapped_reads_sorted_bai: {
            help: "Index of the input coordinate-sorted bam with aligned reads"
        }
        ref_genome:{
            help: "Reference genome.",
            patterns: ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz"]
        }
        methyl_dackel_opt_args:{
            help: "Additional MethylDackel options (https://github.com/dpryan79/methyldackel)"
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

    String output_filename = basename(mapped_reads_sorted_bam, ".bam") 
    String ref_genome_unzip = basename(ref_genome, ".gz")

    Int disk_space = select_first([disk_sp, 2 * ceil(size(mapped_reads_sorted_bam, "G"))])
    Int multi_core = select_first([cores, 1])
    Int mem = multi_core * 4

    command <<<
        set -exo pipefail      
        
        WORKDIR=$HOME

        # Untar reference genome
        gunzip -c ~{ref_genome} > $WORKDIR/ref_genome.fa

        # Index reference sequence in the FASTA format
        samtools faidx $WORKDIR/ref_genome.fa

        # Rename xxx.bam.bai to xxx.bai to make it visible for MethylDackel
        file=~{mapped_reads_sorted_bai}
        mv "$file" "${file//bam.bai/bai}"

        # Extract per-base methylation metrics
        MethylDackel extract \
            $WORKDIR/ref_genome.fa \
            ~{mapped_reads_sorted_bam} \
            -o ~{output_filename} \
            ~{methyl_dackel_opt_args}
            
    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GGFPf1j0pfpFKk7480z3BKPg"])
        cpu: multi_core
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }
    
    output {
        File out_bedgraph = "~{output_filename}"+'_CpG.bedGraph'
    }

}