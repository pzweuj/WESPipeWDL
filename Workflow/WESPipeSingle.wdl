version 1.0
# pzw
# 20220426
# 全外显子分析流程 单人样本

import "QC.wdl" as qc
import "Mapping.wdl" as mapping

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
    call mapping.Bwa as Mapping {input: sample=sample, cleanRead1=QC.cleanRead1, cleanRead2=QC.cleanRead2, threads=threads}
    call mapping.MarkDuplicates as MarkDup {input: sample=sample, sortBam=Mapping.sortBam, sortBamBai=Mapping.sortBamBai}
    call mapping.BQSR as BQSR {input: sample=sample, bam=MarkDup.markBam, bai=MarkDup.markBamBai}
    
    # 质控报告
    call mapping.Bamdst as BamStat {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai, bed=probe}
    call QC.Gender as Gender {input: sample=sample, bam=BQSR.realignBam, bai=BQSR.realignBamBai}
    call QC.YKQCGermline as YKQC {input: sample=sample, fastpJson=QC.jsonReport, bamdstDepth=BamStat.depthReport, bamdstCoverage=BamStat.coverageReport, gender=Gender.genderPredict}

}