version 1.0
# panzhaowen
# 20220413
# WES SNV Indel



################ 线粒体部分 ######################
# Mutect2 MT
# 用于线粒体
# https://github.com/broadinstitute/gatk
task Mutect2MT {
    input {
        String sample
        File bam
        File bai
        Int threads
    }

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

    command <<<
        gatk Mutect2 \
            -R ~{reference} \
            -I ~{bam} \
            -O ~{sample}.MT.vcf \
            --native-pair-hmm-threads ~{threads} \
            -L MT \
            -mbq 15 --force-active true \
            --max-reads-per-alignment-start 0 \
            --mitochondria-mode true \
            --callable-depth 20
    >>>

    output {
        File vcf = "~{sample}.MT.vcf"
        File vcfIdx = "~{sample}.MT.vcf.idx"
        File vcfStats = "~{sample}.MT.vcf.stats"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
        cpus: threads
    }
}

# FilterMutect2
task FilterMutect2MT {
    input {
        String sample
        File vcf
        File vcfIdx
        File vcfStats
    }

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

    command <<<
        gatk FilterMutectCalls \
            -O ~{sample}.MT.filter.vcf \
            -R ~{reference} \
            -V ~{vcf} \
            --mitochondria-mode true
        gatk VariantFiltration \
            -filter "DP < 20" --filter-name "DP20" \
            -V ~{sample}.MT.filter.vcf \
            -O ~{sample}.MT.filterDP.vcf
        gatk SelectVariants \
            -V ~{sample}.MT.filterDP.vcf \
            -O ~{sample}.MT.filtered.vcf \
            --exclude-filtered true
    >>>

    output {
        File filterVcf = "~{sample}.MT.filtered.vcf"
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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

    command <<<
        gatk LeftAlignAndTrimVariants \
            -R ~{reference} \
            -V ~{vcf} \
            -O ~{sample}.HC.leftAlign.vcf \
            --split-multi-allelics --keep-original-ac -no-trim
    >>>

    output {
        File leftVcf = "~{sample}.HC.leftAlign.vcf"
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

# 非GATK的gvcf到vcf
# task GVcf2Vcf {
#     input {
#         String sample
#         File gvcf
#     }

#     command <<<
#         gzip -dc ~{gvcf} | extract_variants \
#             | bgzip -c > ~{sample}.vcf.gz
#     >>>

#     output {
#         File vcf = "~{sample}.vcf.gz"
#     }

# }



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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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
        bcftools view -e "FORMAT/DP<~{minDP} || INFO/AF<~{minAF}" \
            ~{vcf} > ~{sample}.freebayes.filtered.vcf
    >>>

    output {
        File filterVcf = "~{sample}.freebayes.filtered.vcf"
    }

}

# Freebayes作为HC的补充，分析HC中未发现的突变
task FreeSubtractHC {
    input {
        String sample
        File vcfFree
        File vcfHC
    }

    command <<<
        bedtools subtract -a ~{vcfFree} -b ~{vcfHC} > ~{sample}.subtract.vcf
    >>>

    output {
        File subtractVcf = "~{sample}.subtract.vcf"
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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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
        File HCgVcf = "~{sample}.g.vcf.gz"
        File HCgVcfTbi = "~{sample}.g.vcf.gz.tbi"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
        cpus: threads
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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

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

    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File ref_dict = "/home/novelbio/databases/b37/human_g1k_v37_decoy.dict"
    File ref_fai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

    command <<<
        gatk GenotypeGVCFs \
            -R ~{reference} \
            -V ~{gvcf} \
            -O ~{sample}.HC.vcf
    >>>

    output {
        File vcf = "~{sample}.HC.vcf"
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
            -O ~{sample}.indel.vcf.gz \
            -select-type INDEL
        gatk SelectVariants \
            -V ~{vcf} \
            -O ~{sample}.snp.vcf.gz \
            -select-type SNP
        gatk VariantFiltration \
            -V ~{sample}.snp.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -filter "DP < ~{minDP}" --filter-name "DP~{minDP}" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -O ~{sample}.snp.filter.vcf.gz
        gatk VariantFiltration \
            -V ~{sample}.indel.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "SOR > 10.0" --filter-name "SOR10" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -filter "DP < ~{minDP}" --filter-name "DP~{minDP}" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -O ~{sample}.indel.filter.vcf.gz
        gatk MergeVcfs \
            -I ~{sample}.snp.filter.vcf.gz \
            -I ~{sample}.indel.filter.vcf.gz \
            -O ~{sample}.merge.vcf.gz
        gatk SelectVariants \
            -V ~{sample}.merge.vcf.gz \
            -O ~{sample}.HC.filtered.vcf \
            -L ~{bed} \
            --exclude-filtered true
    >>>

    output {
        File filterVcf = "~{sample}.HC.filtered.vcf"
    }

    runtime {
        docker: "broadinstitute/gatk:4.2.6.1"
    }
}

