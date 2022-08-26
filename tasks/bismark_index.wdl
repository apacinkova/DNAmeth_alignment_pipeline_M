version 1.0

task bismark_index {

    input {
        File ref_genome
        String? docker_im
        Int? disk_sp
        Int? actual_mem
    }

    parameter_meta {
        ref_genome:{
            help: "Reference genome.",
            patterns: ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz"]
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

    meta {
        title: "Build bismark reference genome index."
    }

    Int disk_space = select_first([disk_sp, 30])

    command <<<
        set -e +x -o pipefail
        # create directory for reference genome
        WORKDIR=$PWD
        mkdir -p $WORKDIR/ref_genome
        mv "~{ref_genome}" $WORKDIR/ref_genome/
        
        # run Bismark genome preparation
        bismark_genome_preparation $WORKDIR/ref_genome/

        cd $WORKDIR/ref_genome/
        tar -czvf ref_genome_index_bt2.tar.gz Bisulfite_Genome
        mv ref_genome_index_bt2.tar.gz $WORKDIR
    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GG4GPQQ0pfpBf22y8gzP9Z9k"])
        cpu: 8
        memory: select_first([actual_mem, 16]) + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File ref_gen_index = "ref_genome_index_bt2.tar.gz"
    }
}