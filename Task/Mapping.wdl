version 1.0
# panzhaowen
# 20220413
# WES Mapping

# bwa
# https://github.com/lh3/bwa
# 同时使用samblaster进行去重
task Bwa {
    input {
        String sample
        File cleanRead1
        File cleanRead2
        Int threads
    }

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_amb = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.amb"
    File ref_ann = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.ann"
    File ref_bwt = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.bwt"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"
    File ref_pac = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.pac"
    File ref_sa = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.sa"
    Int memory = ceil(threads * 4)

    command <<<
        bwa mem -t ~{threads} \
            -R "@RG\tPL:illumina\tSM:~{sample}\tPU:YK\tID:~{sample}" \
            -v 1 -M ~{reference} ~{cleanRead1} ~{cleanRead2} \
            | sambamba view /dev/stdin -S -h -f bam -t ~{threads} -o ~{sample}.bam
        mkdir ~{sample}_tmp
        sambamba sort \
            -t ~{threads} \
            ~{sample}.bam \
            -o ~{sample}.sort.bam \
            --tmpdir ~{sample}_tmp -m ~{memory}G
        rm ~{sample}.bam*
    >>>

    output {
        File sortBam = "~{sample}.sort.bam"
        File sortBamBai = "~{sample}.sort.bam.bai"
    }

    runtime {
        docker: "pzweuj/mapping:2022May"
        cpus: threads
    }
}

# BWA and samblaster pipe
task BwaMarkDup {
    input {
        String sample
        File cleanRead1
        File cleanRead2
        Int threads
    }

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_amb = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.amb"
    File ref_ann = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.ann"
    File ref_bwt = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.bwt"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"
    File ref_pac = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.pac"
    File ref_sa = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.sa"
    Int memory = ceil(threads * 4)

    command <<<
        bwa mem -t ~{threads} \
            -R "@RG\tPL:illumina\tSM:~{sample}\tPU:YK\tID:~{sample}" \
            -v 1 -M ~{reference} ~{cleanRead1} ~{cleanRead2} \
            | samblaster -M | sambamba view /dev/stdin -S -h -f bam -t ~{threads} -o ~{sample}.bam
        mkdir ~{sample}_tmp
        sambamba sort \
            -t ~{threads} \
            ~{sample}.bam \
            -o ~{sample}.sort.bam \
            --tmpdir ~{sample}_tmp -m ~{memory}G
        rm ~{sample}.bam*
    >>>

    output {
        File sortBam = "~{sample}.sort.bam"
        File sortBamBai = "~{sample}.sort.bam.bai"
    }

    runtime {
        docker: "pzweuj/mapping:2022May"
        cpus: threads
    }
}

# sambamba markdup
task sambambaMarkDup {
    input {
        String sample
        File sortBam
        File sortBamBai
        Int threads
    }

    command <<<
        mkdir ~{sample}_markdups_tmp
        sambamba markdup \
            -t ~{threads} \
            --tmpdir ~{sample}_markdups_tmp \
            ~{sortBam} ~{sample}.marked.bam
    >>>

    output {
        File markBam = "~{sample}.marked.bam"
        File markBamBai = "~{sample}.marked.bam.bai"        
    }

    runtime {
        docker: "pzweuj/mapping:2022May"
        cpus: threads
    }
}

# Markduplicates
# gatk
task MarkDuplicates {
    input {
        String sample
        File sortBam
        File sortBamBai
    }

    command <<<
        mkdir ~{sample}_markdups_tmp
        gatk MarkDuplicates \
            -I ~{sortBam} \
            -O ~{sample}.marked.bam \
            -M ~{sample}.dups.txt \
            --CREATE_INDEX true \
            --TMP_DIR ~{sample}_markdups_tmp
        mv ~{sample}.marked.bai ~{sample}.marked.bam.bai
        rm -rf ~{sample}_markdups_tmp
    >>>

    output {
        File markBam = "~{sample}.marked.bam"
        File markBamBai = "~{sample}.marked.bam.bai"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }
}

# BQSR
task BQSR {
    input {
        String sample
        File bam
        File bai
    }

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File millsIndel = "/home/novelbio/databases/b37/gatk_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
    File millsIndelIdx = "/home/novelbio/databases/b37/gatk_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi"
    File genome1k = "/home/novelbio/databases/b37/gatk_bundle/1000G_phase1.indels.b37.vcf.gz"
    File genome1kIdx = "/home/novelbio/databases/b37/gatk_bundle/1000G_phase1.indels.b37.vcf.gz.tbi"
    File dbsnp = "/home/novelbio/databases/b37/gatk_bundle/dbsnp_138.b37.vcf.gz"
    File dbsnpIdx = "/home/novelbio/databases/b37/gatk_bundle/dbsnp_138.b37.vcf.gz.tbi"

    command <<<
        gatk BaseRecalibrator \
            --known-sites ~{millsIndel} \
            --known-sites ~{genome1k} \
            --known-sites ~{dbsnp} \
            -R ~{reference} \
            -I ~{bam} \
            -O ~{sample}.recal.table
        gatk ApplyBQSR \
            -R ~{reference} \
            --bqsr-recal-file ~{sample}.recal.table \
            -I ~{bam} \
            -O ~{sample}.Realign.bam
        mv ~{sample}.Realign.bai ~{sample}.Realign.bam.bai
    >>>

    output {
        File realignBam = "~{sample}.Realign.bam"
        File realignBamBai = "~{sample}.Realign.bam.bai"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }
}

# bamdst
# https://github.com/shiquan/bamdst
# 用于统计bam比对效果
task Bamdst {
    input {
        String sample
        File bam
        File bai
        File bed
    }

    command <<<
        mkdir ~{sample}_tmp
        bamdst \
            -p ~{bed} \
            -o ~{sample}_tmp \
            ~{bam}
        
        python3 <<CODE
        bamdstReportFile = open("~{sample}_tmp/coverage.report", "r")
        bamdstReport = open("~{sample}.bamdst.txt", "w")
        bamdstReport.write("~{sample} Bamdst QC Report\n")
        for line in bamdstReportFile:
            if line.startswith("#"):
                continue
            else:
                lines = line.lstrip()
                bamdstReport.write(lines)
        bamdstReport.close()
        bamdstReportFile.close()
        CODE
    >>>

    output {
        File bamdstReport = "~{sample}.bamdst.txt"
        File coverageReport = "~{sample}_tmp/coverage.report"
        File depthReport = "~{sample}_tmp/depth.tsv.gz"
        File regionReport = "~{sample}_tmp/region.tsv.gz"
    }

    runtime {
        docker: "pzweuj/bamdst:v1.0.9"
    }
}
