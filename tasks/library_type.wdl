version 1.0

task library_type {

    input {
        File reads_1
        File ref_genome
        File ref_genome_index_tar
        Int num_reads_align = 100000

        String? docker_im
        Int? disk_sp
        Int? cores
    }

    parameter_meta {
        reads_1:{
            help: "Fastq reads 1 in fastq file.",
            patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        }
        ref_genome:{
            help: "Reference genome.",
            patterns: ["*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz"]
        }
        ref_genome_index_tar:{
            help: "Reference genome index (tarball)."
        }
        num_reads_align:{
            help: "Number fo reads from the original fastq file used for library type detection."
        }
        docker_im: {
            help: "Docker image tarball containing execution environment. Please pass as a string in format dx://project-xxx:file-yyy."
        }
        disk_sp: {
            help: "Required disk space (in GB) to run the app."
        }
        cores: {
            help: "Required number of CPU cores to be used for mapping (minimum is 2)."
        }
    }

    meta {
        title: "Library type detection"
        description: "Map bisulfite treated sequencing reads to detect the library type."
    }

    Int disk_space = select_first([disk_sp, ceil(4 * size(reads_1, "G") + 5 * size(ref_genome_index_tar, "G") + size(ref_genome, "G"))])
    Int multi_core = if select_first([cores, 2]) < 2 then 2 else select_first([cores, 2])
    Int mem = multi_core * 24

    command <<<
        set -e +x -o pipefail
        
        # randomly choose num_redas_align reads
        seqtk sample -s100 ~{reads_1} ~{num_reads_align} > random_reads.fastq

        # create directory for reference genome
        mkdir -p $HOME/ref_genome
        mv ~{ref_genome} $HOME/ref_genome/

        # 1) Untar index to $HOME/ref_genome
        tar -xf ~{ref_genome_index_tar} -C $HOME/ref_genome/

        # 2) Library type detection
        bismark --genome $HOME/ref_genome/ \
        --gzip \
        --non_directional \
        --parallel ~{multi_core} \
        random_reads.fastq > bismark_alignment_report.txt

        # number of reads aligned to the original strands
        OT="$(grep '((converted) top strand)' bismark_alignment_report.txt | grep -o -E '[0-9]+' | awk -F "," 'NR!=1{Total=Total+$1} END{print Total}')"
        OB="$(grep '((converted) bottom strand)' bismark_alignment_report.txt | grep -o -E '[0-9]+' | awk -F "," 'NR!=1{Total=Total+$1} END{print Total}')"
        SUM_ORIGINAL=$(($OT + $OB))

        # number of reads aligned to the complementary strands
        CTOT="$(grep '(complementary to (converted) top strand)' bismark_alignment_report.txt | grep -o -E '[0-9]+' | awk -F "," 'NR!=1{Total=Total+$1} END{print Total}')"
        CTOB="$(grep '(complementary to (converted) bottom strand)' bismark_alignment_report.txt | grep -o -E '[0-9]+' | awk -F "," 'NR!=1{Total=Total+$1} END{print Total}')"
        SUM_COMPLEMENTARY=$(($CTOT + $CTOB))

       # if the proportion of reads aligned to the complementary strands > than 25% -> non_directional library
        if python -c "print( round($SUM_COMPLEMENTARY / float($SUM_COMPLEMENTARY + $SUM_ORIGINAL),3) < 0.25)"
            then
                echo 'directional' > library_type
            else
                echo 'non_directional' > library_type
        fi

        mv bismark_alignment_report.txt bismark_alignment_report_library_type_detect.txt

    >>>

    runtime {
        docker: select_first([docker_im, "dx://project-GFBQvF80pfpKzXz1FyzF8Zyj:file-GG1fXFj0pfp5BxPb4bJQK6qy"])
        cpu: multi_core
        memory: mem + "GB"
        disks: "local-disk ${disk_space} SSD"
    }

    output {
        String library_type = read_string("library_type")
        File bismark_alignment_report = "bismark_alignment_report_library_type_detect.txt"
    }
}