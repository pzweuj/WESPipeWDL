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
        File fold80
        File insertsize
    }

    String qcstat = "/home/novelbio/pipeline/WESpipeWDL/Script/qcstat_germline.py"

    command <<<
        python3 ~{qcstat} \
            -i ~{sample} \
            -o ~{sample}.QC.txt \
            -j ~{fastpJson} \
            -b ~{bamdstCoverage} \
            -d ~{bamdstDepth} \
            -g ~{gender} \
            -f ~{fold80} \
            -s ~{insertsize}
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

# 亲缘校检
# 基于plink2
task KingShip {
    input {
        String sample
        File vcf
    }

    String script = "/home/novelbio/pipeline/WESpipeWDL/Script/run_plink2.sh"

    command <<<
        mkdir ~{sample}_kingship
        ~{script} -i ~{vcf} -o ~{sample}_kingship -p ~{sample}
    >>>

    output {
        File kingResults = "~{sample}_kingship/~{sample}.kinship"
    }

}

# 插入片段大小
task InsertSize {
    input {
        String sample
        File bam
        File bai
    }

    command <<<
        gatk CollectInsertSizeMetrics \
            -I ~{bam} \
            -O ~{sample}.insertsize.txt \
            -H ~{sample}.insertsize.pdf
        cat ~{sample}.insertsize.txt | grep -A 1 MEDIAN_INSERT_SIZE | cut -f 1,6 > ~{sample}.insertsize.select.txt
    >>>

    output {
        File insertSizeResults = "~{sample}.insertsize.txt"
        File insertSizeSelect = "~{sample}.insertsize.select.txt"
        File insertSizePDF = "~{sample}.insertsize.pdf"
    }

    runtime {
        docker: "aperdriau/gatk4.2.0:latest"
    }    

}

# HsMetrics
task HsMetrics {
    input {
        String sample
        File bam
        File bai
        File probe
        File bed
    }

    File refDict = "/slurm/databases/b37/human_g1k_v37_decoy.dict"

    command <<<
        gatk BedToIntervalList \
            -I ~{bed} \
            -O b.interval_list \
            -SD ~{refDict}
        gatk BedToIntervalList \
            -I ~{probe} \
            -O p.interval_list \
            -SD ~{refDict}
        gatk CollectHsMetrics \
            -BI p.interval_list \
            -TI b.interval_list \
            -I ~{bam} \
            -O ~{sample}.hs.txt
        cat ~{sample}.hs.txt | grep -A 1 FOLD_80_BASE_PENALTY | cut -f 45 > ~{sample}.fold80.txt
    >>>

    output {
        File hsResults = "~{sample}.hs.txt"
        File fold80Results = "~{sample}.fold80.txt"
    }

    runtime {
        docker: "aperdriau/gatk4.2.0:latest"
    }  

}


