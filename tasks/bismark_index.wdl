version 1.0

task bismark_index {

    input {
        File ref_genome
        String ref_genome_name = "reference_genome"
        String? docker_im
        Int? disk_sp
        Int? actual_mem
    }

    parameter_meta {
        ref_genome:{
            help: "Reference genome."
        }
        docker_im: {
            help: "Docker image containing runtime environment."
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

    command <<<
        set -e +x -o pipefail
        # create directory for reference genome
        mkdir -p $HOME/ref_genome
        mv "~{ref_genome}" $HOME/ref_genome/
        bismark_genome_preparation $HOME/ref_genome/
        tar -czvf "~{ref_genome_name}".tar.gz $PWD

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GFbvjp00pfp8JZv52Zv5xxvP"])
        cpu: 1
        memory: select_first([actual_mem, 16]) + "GB"
        disks: "local-disk ${disk_sp} SSD"
    }

    output {
        File ref_gen_index = "${ref_genome_name}.tar.gz"
    }
}