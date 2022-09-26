version 1.0

task bwameth_index {

    input {
        File ref_genome
        String? docker_im
        Int? disk_sp
        Int? actual_mem
        Int? cores
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
        cores: {
            help: "Required number of CPU cores to be used for mapping."
        }
    }

    meta {
        title: "Build bwa reference genome index for bwa_meth."
    }

    Int disk_space = select_first([disk_sp, 8*ceil(size(ref_genome, "GB"))])

    command <<<
        set -e +x -o pipefail

        # create directory for reference genome index
        WORKDIR=$HOME
        mkdir -p $WORKDIR/index_bwa
        
        # run bwa genome index
        bwameth.py index ~{ref_genome}

        mv ~{ref_genome}.bwa* $WORKDIR/index_bwa/
        tar -C $WORKDIR -czvf "~{ref_genome}"_index_bwa.tar.gz index_bwa
    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GGFPf1j0pfpFKk7480z3BKPg"])
        cpu: select_first([cores, 8])
        memory: select_first([actual_mem, 8]) + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        File ref_gen_index = "${ref_genome}_index_bwa.tar.gz"
    }
}