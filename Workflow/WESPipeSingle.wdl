version 1.0
# pzw
# 20220426
# 全外显子分析流程 单人样本

# Test
# import "../Task/QC.wdl" as qc
# import "../Task/Mapping.wdl" as mapping
# import "../Task/Mutation.wdl" as mutation
# import "../Task/Annotation.wdl" as annotation
# import "../Task/Advance.wdl" as advance
import "QC.wdl" as qc
import "Mapping.wdl" as mapping
import "Mutation.wdl" as mutation
import "Annotation.wdl" as annotation
import "Advance.wdl" as advance

workflow WESPipeSingle {
    input {
        String sample
        File rawRead1
        File rawRead2
        Int threads
        File probe
        File bed
    }

    # 质控
    call qc.Fastp as QC {input: sample=sample, rawRead1=rawRead1, rawRead2=rawRead2, threads=threads}

    # 比对
    call mapping.BwaMarkDup as Mapping {input: sample=sample, cleanRead1=QC.cleanRead1, cleanRead2=QC.cleanRead2, threads=threads}
    # call mapping.MarkDuplicates as MarkDup {input: sample=sample, sortBam=Mapping.sortBam, sortBamBai=Mapping.sortBamBai}
    call mapping.BQSR as BQSR {input: sample=sample, bam=Mapping.sortBam, bai=Mapping.sortBamBai}
    
    # 质控报告
    call mapping.Bamdst as BamStat {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai, bed=probe}
    call qc.Gender as Gender {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai}
    call qc.YKQCGermline as YKQC {input: sample=sample, fastpJson=QC.jsonReport, bamdstDepth=BamStat.depthReport, bamdstCoverage=BamStat.coverageReport, gender=Gender.genderPredict}

    # 线粒体检测
    call mutation.Mutect2MT as MT {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai, threads=threads}
    call mutation.FilterMutect2MT as MTFilter {input: sample=sample, vcf=MT.vcf, vcfIdx=MT.vcfIdx, vcfStats=MT.vcfStats}
    call mutation.LeftAlignMutect2 as MTLeft {input: sample=sample, vcf=MTFilter.filterVcf}

    # 变异检测
    ## HC 流程
    call mutation.HaplotypeCaller as HaplotypeCaller {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai, threads=threads, bed=bed}
    call mutation.GenotypeGVCFs as GenotypeGVCFs {input: sample=sample, gvcf=HaplotypeCaller.HCgVcf, gvcfTbi=HaplotypeCaller.HCgVcfTbi}
    call mutation.FilterHaplotypeCaller as FilterHaplotypeCaller {input: sample=sample, vcf=GenotypeGVCFs.vcf, bed=bed, minDP=10}
    call mutation.LeftAlignMutect2 as HCLeft {input: sample=sample, vcf=FilterHaplotypeCaller.filterVcf}
    ## Freebayes 流程
    call mutation.Freebayes as Freebayes {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai, bed=bed}
    call mutation.FilterFreebayes as FilterFreebayes {input: sample=sample, vcf=Freebayes.vcf, minDP=10, minAF=0.2}
    call mutation.LeftAlignBcftools as FreebayesLeft {input: sample=sample, vcf=FilterFreebayes.filterVcf}
    ## Freebayes HC差异提取
    call mutation.FreeSubtractHC as FreeSubtractHC {input: sample=sample, vcfFree=FreebayesLeft.leftVcf, vcfHC=HCLeft.leftVcf}

    # 注释
    call annotation.GeneCoverage as GeneCoverage {input: sample=sample, regionGz=BamStat.regionReport, bed=bed}
    ## 线粒体注释
    call annotation.MTAnnotation as MTAnnotation {input: sample=sample, vcf=MTLeft.leftVcf, threads=threads}
    ## HC结果注释
    call annotation.Snpeff as HCSnpeff {input: sample=sample, vcf=HCLeft.leftVcf}
    call annotation.Annovar as HCAnnotation {input: sample=sample, vcf=HCSnpeff.annoVcf, threads=threads}
    call annotation.AnnotationFix as HCAnnoFix {input: sample=sample, annoFile=HCAnnotation.annovarResults, geneCoverFile=GeneCoverage.geneCover}
    ## 差异结果注释
    call annotation.Snpeff as SubSnpeff {input: sample=sample, vcf=FreeSubtractHC.subtractVcf}
    call annotation.Annovar as SubAnnotation {input: sample=sample, vcf=SubSnpeff.annoVcf, threads=threads}
    call annotation.AnnotationFix as SubAnnoFix {input: sample=sample, annoFile=SubAnnotation.annovarResults, geneCoverFile=GeneCoverage.geneCover}
    
    # 动态突变
    call advance.ExpansionHunter as ExpansionHunter {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai}
    

}