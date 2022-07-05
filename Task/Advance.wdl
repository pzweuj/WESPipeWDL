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

# UPD analysis
# 仅适用与父母子三人样本
# vcf文件需包含人群频率AF信息，该信息储存于INFO的CSQ tag下，参考VEP注释结果
# https://github.com/bjhall/upd
task UPD {
    input {
        String proband
        String father
        String mother
        File vcf
    }

    command <<<
        upd --vcf ~{vcf} \
            --proband ~{proband} \
            --father ~{father} \
            --mother ~{mother} \
            --af-tag AF \
            regions > ~{proband}.trio.upd.txt
    >>>

    output {
        File updResults = "~{proband}.trio.upd.txt"
    }

}

