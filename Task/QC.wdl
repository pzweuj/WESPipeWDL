version 1.0
# 20220412
# panzhaowen
# WES QC

# fastp
# https://github.com/OpenGene/fastp
task Fastp {
    input {
        String sample
        File rawRead1
        File rawRead2
        Int threads
    }

    command <<<
        fastp \
            -i ~{rawRead1} \
            -I ~{rawRead2} \
            -o ~{sample}.clean_R1.fastq.gz \
            -O ~{sample}.clean_R2.fastq.gz \
            -w ~{threads} \
            -j ~{sample}.json \
            -h ~{sample}.html -t 1
    >>>

    output {
        File cleanRead1 = "~{sample}.clean_R1.fastq.gz"
        File cleanRead2 = "~{sample}.clean_R2.fastq.gz"
        File jsonReport = "~{sample}.json"
        File htmlReport = "~{sample}.html"
    }

    runtime {
        docker: "pzweuj/gencore:v2"
        cpus: threads
    }

}

# 胚系样本QC整理
task YKQCGermline {
    input {
        String sample
        File fastpJson
        File bamdstDepth
        File bamdstCoverage
        File gender
    }

    File qcstat = "/home/novelbio/pipeline/WESpipeWDL/Script/qcstat_germline.py"

    command <<<
        python3 ~{qcstat} ~{sample} ~{fastpJson} ~{bamdstCoverage} ~{bamdstDepth} ~{gender} ~{sample}.QC.txt
    >>>

    output {
        File YKQCReport = "~{sample}.QC.txt"
    }
}

# 性别验证
task Gender {
    input {
        String sample
        File bam
        File bai
    }

    command <<<
        samtools view ~{bam} Y:2654896-2655723 | wc -l > ~{sample}.SRY.counts.txt
    >>>

    output {
        File genderPredict = "~{sample}.SRY.counts.txt"
    }

}


