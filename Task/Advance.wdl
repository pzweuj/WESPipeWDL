version 1.0
# panzhaowen
# 20220523
# 其他分析

# ExpansionHunter
# https://github.com/Illumina/ExpansionHunter
task ExpansionHunter {
    input {
        String sample
        File bam
        File bai
    }

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File refFai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

    command <<<
        ExpansionHunterPipe.py -i ~{bam} \
            -o ~{sample} \
            -r ~{reference} \
            -p false
    >>>

    output {
        File ehResults = "~{sample}.EH.txt"
    }

    runtime {
        docker: "pzweuj/expansionhunter:5.0.0"
    }

}



