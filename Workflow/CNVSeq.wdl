version 1.0
# pzw
# 20220927
# CNV-seq 低深度全基因组
# 华大平台 PE100

# Test
# import "../Task/QC.wdl" as qc
# import "../Task/Mapping.wdl" as mapping
# import "../Task/Annotation.wdl" as annotation
# import "../Task/CNV.wdl" as cnv
# import "../Task/Advance.wdl" as advance

import "QC.wdl" as qc
import "Mapping.wdl" as mapping
import "Annotation.wdl" as annotation
import "CNV.wdl" as cnv
# import "Advance.wdl" as advance

workflow CNVSeqPipe {
    input {
        String sample
        String sex
        File rawRead1
        File rawRead2
        Int threads
    }

    # 质控
    call qc.Fastp as QC {input: sample=sample, rawRead1=rawRead1, rawRead2=rawRead2, threads=threads}

    # 比对
    call mapping.Bwa as Mapping {input: sample=sample, cleanRead1=QC.cleanRead1, cleanRead2=QC.cleanRead2, threads=threads}
    call mapping.MarkDuplicatesRM as RemoveDup {input: sample=sample, sortBam=Mapping.sortBam, sortBamBai=Mapping.sortBamBai}

    # 质控报告
    call qc.InsertSize as InsertSize {input: sample=sample, bam=RemoveDup.rmDupBam, bai=RemoveDup.rmDupBamBai}
    call qc.Genomecov as Genomecov {input: sample=sample, bam=RemoveDup.rmDupBam, bai=RemoveDup.rmDupBamBai}
    call qc.Gender as Gender {input: sample=sample, bam=RemoveDup.rmDupBam, bai=RemoveDup.rmDupBamBai}
    call qc.YKQCWGS as YKQC {input: sample=sample, fastpJson=QC.jsonReport, covFile=Genomecov.covFile, dupFile=RemoveDup.dupResult, gender=Gender.genderPredict, insertsize=InsertSize.insertSizeSelect}

    # CNV分析
    call cnv.Freec as Freec {input: sample=sample, bam=RemoveDup.rmDupBam, bai=RemoveDup.rmDupBamBai, threads=threads, sex=sex}
    call cnv.FreecPlot as FreecPlot {input: sample=sample, ratioFile=Freec.ratioResult, cnvResult=Freec.cnvResult}

    # CNV注释
    call annotation.AnnotSV as AnnotSV {input: sample=sample, cnvResult=Freec.cnvResult}

}


