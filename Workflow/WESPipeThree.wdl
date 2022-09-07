version 1.0
# pzw
# 20220714
# 全外显子分析流程 三人家系

# Test
import "../Task/QC.wdl" as qc
import "../Task/Mapping.wdl" as mapping
import "../Task/Mutation.wdl" as mutation
import "../Task/Annotation.wdl" as annotation
import "../Task/Advance.wdl" as advance
# import "QC.wdl" as qc
# import "Mapping.wdl" as mapping
# import "Mutation.wdl" as mutation
# import "Annotation.wdl" as annotation
# import "Advance.wdl" as advance

workflow WESPipeThree {
    input {
        String proband
        String father
        String mother
        File probandRead1
        File probandRead2
        File fatherRead1
        File fatherRead2
        File motherRead1
        File motherRead2
        File probe
        File bed
        Int threads
    }

    # 质控
    call qc.Fastp as QCProband {input: sample=proband, rawRead1=probandRead1, rawRead2=probandRead2, threads=threads}
    call qc.Fastp as QCFather {input: sample=father, rawRead1=fatherRead1, rawRead2=fatherRead2, threads=threads}
    call qc.Fastp as QCMother {input: sample=mother, rawRead1=motherRead1, rawRead2=motherRead2, threads=threads}

    # 比对
    call mapping.Bwa as MappingProband {input: sample=proband, cleanRead1=QCProband.cleanRead1, cleanRead2=QCProband.cleanRead2, threads=threads}
    call mapping.Bwa as MappingFather {input: sample=father, cleanRead1=QCFather.cleanRead1, cleanRead2=QCFather.cleanRead2, threads=threads}
    call mapping.Bwa as MappingMother {input: sample=mother, cleanRead1=QCMother.cleanRead1, cleanRead2=QCMother.cleanRead2, threads=threads}
    call mapping.MarkDuplicates as MarkDupProband {input: sample=proband, sortBam=MappingProband.sortBam, sortBamBai=MappingProband.sortBamBai}
    call mapping.MarkDuplicates as MarkDupFather {input: sample=father, sortBam=MappingFather.sortBam, sortBamBai=MappingFather.sortBamBai}
    call mapping.MarkDuplicates as MarkDupMother {input: sample=mother, sortBam=MappingMother.sortBam, sortBamBai=MappingMother.sortBamBai}
    call mapping.BQSR as BQSRProband {input: sample=proband, bam=MarkDupProband.markBam, bai=MarkDupProband.markBamBai}
    call mapping.BQSR as BQSRFather {input: sample=father, bam=MarkDupFather.markBam, bai=MarkDupFather.markBamBai}
    call mapping.BQSR as BQSRMother {input: sample=mother, bam=MarkDupMother.markBam, bai=MarkDupMother.markBamBai}

    # 质控报告
    call mapping.Bamdst as BamStatProband {input: sample=proband, bam=MarkDupProband.markBam, bai=MarkDupProband.markBamBai, bed=probe}
    call mapping.Bamdst as BamStatFather {input: sample=father, bam=MarkDupFather.markBam, bai=MarkDupFather.markBamBai, bed=probe}
    call mapping.Bamdst as BamStatMother {input: sample=mother, bam=MarkDupMother.markBam, bai=MarkDupMother.markBamBai, bed=probe}
    call qc.Gender as GenderProband {input: sample=proband, bam=MarkDupProband.markBam, bai=MarkDupProband.markBamBai}
    call qc.Gender as GenderFather {input: sample=father, bam=MarkDupFather.markBam, bai=MarkDupFather.markBamBai}
    call qc.Gender as GenderMother {input: sample=mother, bam=MarkDupMother.markBam, bai=MarkDupMother.markBamBai}
    call qc.InsertSize as InsertSizeProband {input: sample=proband, bam=MarkDupProband.markBam, bai=MarkDupProband.markBamBai}
    call qc.InsertSize as InsertSizeFather {input: sample=father, bam=MarkDupFather.markBam, bai=MarkDupFather.markBamBai}
    call qc.InsertSize as InsertSizeMother {input: sample=mother, bam=MarkDupMother.markBam, bai=MarkDupMother.markBamBai}
    call qc.HsMetrics as HsMetricsProband {input: sample=proband, bam=MarkDupProband.markBam, bai=MarkDupProband.markBamBai, probe=probe, bed=bed}
    call qc.HsMetrics as HsMetricsFather {input: sample=father, bam=MarkDupFather.markBam, bai=MarkDupFather.markBamBai, probe=probe, bed=bed}
    call qc.HsMetrics as HsMetricsMother {input: sample=mother, bam=MarkDupMother.markBam, bai=MarkDupMother.markBamBai, probe=probe, bed=bed}
    call qc.YKQCGermline as YKQCProband {
        input:
            sample=proband,
            fastpJson=QCProband.jsonReport,
            bamdstDepth=BamStatProband.depthReport,
            bamdstCoverage=BamStatProband.coverageReport,
            gender=GenderProband.genderPredict,
            fold80=HsMetricsProband.fold80Results,
            insertsize=InsertSizeProband.insertSizeSelect
    }
    call qc.YKQCGermline as YKQCFather {
        input:
            sample=father,
            fastpJson=QCFather.jsonReport,
            bamdstDepth=BamStatFather.depthReport,
            bamdstCoverage=BamStatFather.coverageReport,
            gender=GenderFather.genderPredict,
            fold80=HsMetricsFather.fold80Results,
            insertsize=InsertSizeFather.insertSizeSelect
    }
    call qc.YKQCGermline as YKQCMother {
        input:
            sample=mother,
            fastpJson=QCMother.jsonReport,
            bamdstDepth=BamStatMother.depthReport,
            bamdstCoverage=BamStatMother.coverageReport,
            gender=GenderMother.genderPredict,
            fold80=HsMetricsMother.fold80Results,
            insertsize=InsertSizeMother.insertSizeSelect
    }

    # 线粒体检测（仅先证者）
    call mutation.Mutect2MT as MT {input: sample=proband, bam=BQSRProband.realignBam, bai=BQSRProband.realignBamBai, threads=threads}
    call mutation.FilterMutect2MT as MTFilter {input: sample=proband, vcf=MT.vcf, vcfIdx=MT.vcfIdx, vcfStats=MT.vcfStats}
    call mutation.LeftAlignMutect2 as MTLeft {input: sample=proband, vcf=MTFilter.filterVcf}

    # 变异检测
    ## HC流程
    call mutation.HaplotypeCaller as HCProband {input: sample=proband, bam=BQSRProband.realignBam, bai=BQSRProband.realignBamBai, threads=threads, bed=bed}
    call mutation.HaplotypeCaller as HCFather {input: sample=father, bam=BQSRFather.realignBam, bai=BQSRFather.realignBamBai, threads=threads, bed=bed}
    call mutation.HaplotypeCaller as HCMother {input: sample=mother, bam=BQSRMother.realignBam, bai=BQSRMother.realignBamBai, threads=threads, bed=bed}
    call mutation.CombineGVCFsThree as CombineGVCFs {
        input:
            proband=proband,
            probandGVCF=HCProband.HCgVcf,
            probandGVCFTbi=HCProband.HCgVcfTbi,
            p2GVCF=HCFather.HCgVcf,
            p2GVCFTbi=HCFather.HCgVcfTbi,
            p3GVCF=HCMother.HCgVcf,
            p3GVCF=HCMother.HCgVcfTbi
    }
    call mutation.GenotypeGVCFs as GenotypeGVCFs {input: sample=proband, gvcf=CombineGVCFs.combineGVCF, gvcfTbi=CombineGVCFs.combineGVCFTbi}
    call mutation.FilterHaplotypeCaller as FilterHaplotypeCaller {input: sample=proband, vcf=GenotypeGVCFs.vcf, bed=bed, minDP=10}
    call mutation.LeftAlignMutect2 as HCLeft {input: sample=proband, vcf=FilterHaplotypeCaller.filterVcf}
    ## Freebayes流程
    call mutation.FreebayesThree as Freebayes {
        input:
            proband=proband,
            probandBam=BQSRProband.realignBam,
            probandBamBai=BQSRProband.realignBamBai,
            p2Bam=BQSRFather.realignBam,
            p2BamBai=BQSRFather.realignBamBai,
            p3Bam=BQSRMother.realignBam,
            p2BamBai=BQSRMother.realignBamBai,
            bed=bed
    }
    call mutation.FilterFreebayes as FilterFreebayes {input: sample=proband, vcf=Freebayes.vcf, minDP=10, minAF=0.2}
    call mutation.LeftAlignBcftools as FreebayesLeft {input: sample=proband, vcf=FilterFreebayes.filterVcf}
    ## Freebayes/HC差异提取
    call mutation.FreeSubtractHC as FreeSubtractHC {input: sample=proband, vcfFree=FreebayesLeft.leftVcf, vcfHC=HCLeft.leftVcf}

    # 注释
    call annotation.GeneCoverage as GeneCoverage {input: sample=proband, regionGz=BamStatProband.regionReport, bed=bed}
    ## 线粒体注释
    call annotation.MTAnnotation as MTAnnotation {input: sample=proband, vcf=MTLeft.leftVcf}
    ## HC结果注释
    call annotation.Snpeff as HCSnpeff {input: sample=proband, vcf=HCLeft.leftVcf}
    call annotation.Annovar as HCAnnotation {input: sample=proband, vcf=HCSnpeff.annoVcf}
    call annotation.AnnotationFix as HCAnnoFix {input: sample=proband, annoFile=HCAnnotation.annovarResults, geneCoverFile=GeneCoverage.geneCover}
    ## 差异结果注释
    String subProband = proband + "_subFreebayes"
    call annotation.Snpeff as SubSnpeff {input: sample=subProband, vcf=FreeSubtractHC.subtractVcf}
    call annotation.Annovar as SubAnnotation {input: sample=subProband, vcf=SubSnpeff.annoVcf}
    call annotation.AnnotationFix as SubAnnoFix {input: sample=subProband, annoFile=SubAnnotation.annovarResults, geneCoverFile=GeneCoverage.geneCover}
    
    # 动态突变
    call advance.ExpansionHunter as ExpansionHunter {input: sample=proband, bam=BQSRProband.realignBam, bai=BQSRProband.realignBamBai}

}

