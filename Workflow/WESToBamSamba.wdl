version 1.0
# pzw
# 20220426
# 全外显子分析流程 仅比对

# Test
# import "../Task/QC.wdl" as qc
# import "../Task/Mapping.wdl" as mapping
# import "../Task/Mutation.wdl" as mutation
# import "../Task/Annotation.wdl" as annotation
import "QC.wdl" as qc
import "Mapping.wdl" as mapping
import "Mutation.wdl" as mutation
import "Annotation.wdl" as annotation

workflow WESPipeToBam {
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
    call mapping.Bwa as Mapping {input: sample=sample, cleanRead1=QC.cleanRead1, cleanRead2=QC.cleanRead2, threads=threads}
    call mapping.sambambaMarkDup as MarkDup {input: sample=sample, sortBam=Mapping.sortBam, sortBamBai=Mapping.sortBamBai, threads=threads}
    call mapping.BQSR as BQSR {input: sample=sample, bam=MarkDup.markBam, bai=MarkDup.markBamBai}
    
    # 质控报告
    call mapping.Bamdst as BamStat {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai, bed=probe}
    call qc.Gender as Gender {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai}
    call qc.YKQCGermline as YKQC {input: sample=sample, fastpJson=QC.jsonReport, bamdstDepth=BamStat.depthReport, bamdstCoverage=BamStat.coverageReport, gender=Gender.genderPredict}
}