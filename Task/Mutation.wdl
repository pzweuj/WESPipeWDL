version 1.0
# panzhaowen
# 20220413
# WES SNV Indel



################ 线粒体部分 ######################
# Mutect2 MT
# 用于线粒体
# https://github.com/broadinstitute/gatk
task Mutect2 {
    input {
        File bam
        File bai
        String sample
        Int threads
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk Mutect2 \
            -R ~{reference} \
            -I ~{bam} \
            -O ~{sample}.chrM.vcf \
            --native-pair-hmm-threads ~{threads} \
            -L chrM \
            -A Coverage -A GenotypeSummaries \
            -mbq 15 --force-active true \
            --max-reads-per-alignment-start 0 \
            --mitochondria-mode
    >>>

    output {
        File vcf = "~{sample}.chrM.vcf"
        File vcfIdx = "~{sample}.chrM.vcf.idx"
        File vcfStats = "~{sample}.chrM.vcf.stats"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
        cpus: threads
    }
}

# FilterMutect2
task FilterMutect2 {
    input {
        String sample
        File vcf
        File vcfIdx
        File vcfStats
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk FilterMutectCalls \
            -O ~{sample}.mutect2.multiAelle.filter.vcf \
            -R ~{reference} \
            -V ~{vcf} \
            --mitochondria-mode true
    >>>

    output {
        File filterVcf = "~{sample}.mutect2.filter.vcf"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }

}


##################### 通用模块 ########################
# 左对齐及拆分突变
task LeftAlignMutect2 {
    input {
        String sample
        File vcf
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk LeftAlignAndTrimVariants \
            -R ~{reference} \
            -V ~{vcf} \
            -O ~{sample}.leftAlign.vcf \
            --split-multi-allelics --keep-original-ac -no-trim
    >>>

    output {
        File leftVcf = "~{sample}.leftAlign.vcf"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }

}

# 左对齐及拆分突变 bcftools
task LeftAlignBcftools {
    input {
        String sample
        File vcf
    }

    command <<<
        bcftools norm -m -any ~{vcf} \
            > ~{sample}.leftAlign.vcf    
    >>>

    output {
        File leftVcf = "~{sample}.leftAlign.vcf"
    }

    runtime {
        docker: "llcondocker/freebayes:1.3.4"
    }

}


######################## freebayes pipe ###########################
# freebayes
# https://github.com/freebayes/freebayes
task Freebayes {
    input {
        String sample
        File bam
        File bai
        File bed
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        freebayes -f ~{reference} \
            ~{bam} \
            -t ~{bed} \
            --min-coverage 10 \
            --genotyping-max-banddepth 2 \
            > ~{sample}.freebayes.vcf
    >>>

    output {
        File vcf = "~{sample}.freebayes.vcf"
    }

    runtime {
        docker: "llcondocker/freebayes:1.3.4"
    }
}

# freebayes 二人
task FreebayesTwo {
    input {
        String proband
        File probandBam
        File probandBamBai
        File p2Bam
        File p2BamBai
        File bed
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        freebayes -f ~{reference} \
            ~{probandBam} ~{p2Bam} \
            -t ~{bed} \
            --min-coverage 10 \
            --genotyping-max-banddepth 2 \
            > ~{proband}.trio.freebayes.vcf
    >>>

    output {
        File vcf = "~{proband}.freebayes.vcf"
    }

    runtime {
        docker: "llcondocker/freebayes:1.3.4"
    }
}

# freebayes 三人
task FreebayesThree {
    input {
        String proband
        File probandBam
        File probandBamBai
        File p2Bam
        File p2BamBai
        File p3Bam
        File p3BamBai
        File bed
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        freebayes -f ~{reference} \
            ~{probandBam} ~{p2Bam} ~{p3Bam} \
            -t ~{bed} \
            --min-coverage 10 \
            --genotyping-max-banddepth 2 \
            > ~{proband}.trio.freebayes.vcf
    >>>

    output {
        File vcf = "~{proband}.freebayes.vcf"
    }

    runtime {
        docker: "llcondocker/freebayes:1.3.4"
    }
}

# freebayes 四人
task FreebayesFour {
    input {
        String proband
        File probandBam
        File probandBamBai
        File p2Bam
        File p2BamBai
        File p3Bam
        File p3BamBai
        File p4Bam
        File p4BamBai
        File bed
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        freebayes -f ~{reference} \
            ~{probandBam} ~{p2Bam} ~{p3Bam} ~{p4Bam} \
            -t ~{bed} \
            --min-coverage 10 \
            --genotyping-max-banddepth 2 \
            > ~{proband}.trio.freebayes.vcf
    >>>

    output {
        File vcf = "~{proband}.freebayes.vcf"
    }

    runtime {
        docker: "llcondocker/freebayes:1.3.4"
    }
}

# 过滤
task FilterFreebayes {
    input {
        String sample
        File vcf
        Int minDP
        Float minAF
    }

    command <<<
        bcftools view -e "FORMAT/DP<~{minDP} || FORMAT/AF<~{minAF}" \
            ~{vcf} > ~{sample}.filtered.vcf
    >>>

    output {
        File filterVcf = "~{sample}.filtered.vcf"
    }

}


###################### HaplotypeCaller Pipe #######################
# gatk胚系流程
task HaplotypeCaller {
    input {
        String sample
        File bam
        File bai
        Int threads
        File? bed
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk HaplotypeCaller \
            -R ~{reference} \
            -I ~{bam} \
            -O ~{sample}.g.vcf.gz \
            --native-pair-hmm-threads ~{threads} \
            -L ~{bed} \
            -ERC GVCF
    >>>

    output {
        File combineGVCF = "~{sample}.g.vcf.gz"
        File combineGVCFTbi = "~{sample}.g.vcf.gz.tbi"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }
}

# 家系合并
# 二人
task CombineGVCFsTwo {
    input {
        String proband
        File probandGVCF
        File probandGVCFTbi
        File p2GVCF
        File p2GVCFTbi
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk CombineGVCFs \
            -R ~{reference} \
            -V ~{probandGVCF} \
            -V ~{p2GVCF} \
            -O ~{proband}.trio.g.vcf.gz
    >>>

    output {
        File combineGVCF = "~{proband}.trio.g.vcf.gz"
        File combineGVCFTbi = "~{proband}.trio.g.vcf.gz.tbi"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }

}

# 三人
task CombineGVCFsThree {
    input {
        String proband
        File probandGVCF
        File probandGVCFTbi
        File p2GVCF
        File p2GVCFTbi
        File p3GVCF
        File p3GVCFTbi
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk CombineGVCFs \
            -R ~{reference} \
            -V ~{probandGVCF} \
            -V ~{p2GVCF} \
            -V ~{p3GVCF} \
            -O ~{proband}.trio.g.vcf.gz
    >>>

    output {
        File combineGVCF = "~{proband}.trio.g.vcf.gz"
        File combineGVCFTbi = "~{proband}.trio.g.vcf.gz.tbi"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }

}

# 四人
task CombineGVCFsFour {
    input {
        String proband
        File probandGVCF
        File probandGVCFTbi
        File p2GVCF
        File p2GVCFTbi
        File p3GVCF
        File p3GVCFTbi
        File p4GVCF
        File p4GVCFTbi
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk CombineGVCFs \
            -R ~{reference} \
            -V ~{probandGVCF} \
            -V ~{p2GVCF} \
            -V ~{p3GVCF} \
            -V ~{p4GVCF} \
            -O ~{proband}.trio.g.vcf.gz
    >>>

    output {
        File combineGVCF = "~{proband}.trio.g.vcf.gz"
        File combineGVCFTbi = "~{proband}.trio.g.vcf.gz.tbi"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }

}

# 从gvcf到vcf
task GenotypeGVCFs {
    input {
        String sample
        File gvcf
        File gvcfTbi
    }

    File reference = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa"
    File ref_dict = "/home/novelbio/databases/hs37d5/hs37d5.chr.dict"
    File ref_fai = "/home/novelbio/databases/hs37d5/hs37d5.chr.fa.fai"

    command <<<
        gatk GenotypeGVCFs \
            -R ~{reference} \
            -V ~{gvcf} \
            -O ~{sample}.HaplotypeCaller.vcf
    >>>

    output {
        File vcf = "~{sample}.HaplotypeCaller.vcf"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }
}

# 过滤
task FilterHaplotypeCaller {
    input {
        String sample
        File vcf
        File bed
        Int minDP
    }

    command <<<
        gatk SelectVariants \
            -V ~{vcf} \
            -O ~{sample}.indel.vcf \
            -select-type INDEL
        gatk SelectVariants \
            -V ~{vcf} \
            -O ~{sample}.snp.vcf \
            -select-type SNP
        gatk VariantFiltration \
            -V ~{sample}.snp.vcf \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -filter "DP < ~{minDP}" --filter-name "DP~{minDP}" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -O ~{sample}.snp.filter.vcf
        gatk VariantFiltration \
            -V ~{sample}.indel.vcf \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "SOR > 10.0" --filter-name "SOR10" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -filter "DP < ~{minDP}" --filter-name "DP~{minDP}" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -O ~{sample}.indel.filter.vcf
        gatk MergeVcfs \
            -I ~{sample}.snp.filter.vcf \
            -I ~{sample}.indel.filter.vcf \
            -O ~{sample}.merge.vcf
        gatk SelectVariants \
            -V ~{sample}.merge.vcf \
            -O ~{sample}.filtered.vcf \
            -L ~{bed} \
            --exclude-filtered true
    >>>

    output {
        File filterVcf = "~{sample}.filtered.vcf"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }
}

