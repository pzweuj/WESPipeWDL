# coding=utf-8
# pzw
# 20220629

"""
SampleID|rawReads|rawBases|cleanReads|cleanBases|duplicatesRate|cleanQ20|cleanQ30|cleanGC
|mappedReads|mappedBases|mappedRate|targetReads|targetBases|targetRate|depth|depthRmdups|Uniformity|Fold80
|1XCov|10XCov|20XCov|50XCov|1XCov_rmdup|10XCov_rmdup|20XCov_rmdup|50XCov_rmdup|medianInsertSize|meanInsertSize|PredictGender
"""

import sys
import json
import gzip
import argparse

# json解析
# fastpJsonReport = "202111.json"
def jsonAnalysis(fastpJsonReport):
    jsonFile = json.load(open(fastpJsonReport, "r"))
    rawReads = jsonFile["summary"]["before_filtering"]["total_reads"]
    rawBases = jsonFile["summary"]["before_filtering"]["total_bases"]
    cleanReads = jsonFile["summary"]["after_filtering"]["total_reads"]
    cleanBases = jsonFile["summary"]["after_filtering"]["total_bases"]
    # duplicationRate = "%.2f" % (jsonFile["duplication"]["rate"] * 100) + "%"
    cleanQ20Bases = jsonFile["summary"]["after_filtering"]["q20_bases"]
    cleanQ20 = "%.2f" % (cleanQ20Bases / cleanBases * 100) + "%"
    cleanQ30Bases = jsonFile["summary"]["after_filtering"]["q30_bases"]
    cleanQ30 = "%.2f" % (cleanQ30Bases / cleanBases * 100) + "%"
    cleanGC = "%.2f" % (jsonFile["summary"]["after_filtering"]["gc_content"] * 100) + "%"
    readLength1 = jsonFile["summary"]["after_filtering"]["read1_mean_length"]
    readLength2 = jsonFile["summary"]["after_filtering"]["read2_mean_length"]
    readLength = (readLength1 + readLength2) / 2
    return [rawReads, rawBases, cleanReads, cleanBases, cleanQ20, cleanQ30, cleanGC, readLength]

# bamdst解析
# bamdstReport = "coverage.report"
def coverageAnalysis(bamdstReport, cleanBases, readLength):
    bamdstCoverageReport = open(bamdstReport, "r")
    mappedReads = 0
    targetReads = 0
    averageDepth = 0
    averageRmDepth = 0
    for line in bamdstCoverageReport:
        if "[Total] Mapped Reads" in line:
            mappedReads = int(line.replace("\n", "").split("\t")[1])
        elif "[Target] Target Reads" in line:
            targetReads = int(line.replace("\n", "").split("\t")[1])
        elif "[Target] Average depth(rmdup)" in line:
            averageRmDepth = float(line.replace("\n", "").split("\t")[1])
        elif "[Target] Average depth" in line:
            averageDepth = float(line.replace("\n", "").split("\t")[1])
        else:
            continue
    bamdstCoverageReport.close()
    mappedBases = int(mappedReads * readLength)
    mappedRate = "%.2f" % (mappedBases / cleanBases * 100) + "%"
    targetBases = int(targetReads * readLength)
    targetRate = "%.2f" % (targetBases / mappedBases * 100) + "%"

    # 重复率
    if averageDepth != 0:
        duplicationRate = "%.2f" % (((averageDepth - averageRmDepth) / averageDepth) * 100) + "%"
    else:
        duplicationRate = "-"
    return [duplicationRate, mappedReads, mappedBases, mappedRate, targetReads, targetBases, targetRate, averageDepth, averageRmDepth]

# 以reads计算比对率及捕获率
# mappedRate = "%.2f" % (mappedReads / cleanReads * 100) + "%"
# targetRate = "%.2f" % (targetReads / mappedReads * 100) + "%"

# 深度
# depthFile = "depth.tsv.gz"
def uniformityAnalysis(depthFile, averageDepth):
    depthReport = gzip.open(depthFile, "rt")
    targetLength = 0
    uniform = 0
    C1 = C10 = C20 = C50 = 0
    C1_rm = C10_rm = C20_rm = C50_rm = 0
    targetDepthBases = 0
    targetRmDepthBases = 0
    for line in depthReport:
        if not line.startswith("#"):
            targetLength += 1
            lines = line.split("\t")
            # chrom = lines[0]
            # pos = lines[1]
            rawDepth = int(lines[2])
            rmDepth = int(lines[3])
            targetDepthBases += rawDepth
            targetRmDepthBases += rmDepth

            if rawDepth >= (averageDepth * 0.2):
                uniform += 1
            if rawDepth >= 50:
                C50 += 1
                C20 += 1
                C10 += 1
                C1 += 1
            elif rawDepth >= 20:
                C20 += 1
                C10 += 1
                C1 += 1
            elif rawDepth >= 10:
                C10 += 1
                C1 += 1
            elif rawDepth >= 1:
                C1 += 1
            else:
                pass

            # 去重
            if rmDepth >= 50:
                C50_rm += 1
                C20_rm += 1
                C10_rm += 1
                C1_rm += 1
            elif rmDepth >= 20:
                C20_rm += 1
                C10_rm += 1
                C1_rm += 1
            elif rmDepth >= 10:
                C10_rm += 1
                C1_rm += 1
            elif rmDepth >= 1:
                C1_rm += 1
            else:
                pass


    uniformity = "%.2f" % (uniform / targetLength * 100) + "%"
    cover50 = "%.2f" % (C50 / targetLength * 100) + "%"
    cover20 = "%.2f" % (C20 / targetLength * 100) + "%"
    cover10 = "%.2f" % (C10 / targetLength * 100) + "%"
    cover1 = "%.2f" % (C1 / targetLength * 100) + "%"
    cover50_rm = "%.2f" % (C50_rm / targetLength * 100) + "%"
    cover20_rm = "%.2f" % (C20_rm / targetLength * 100) + "%"
    cover10_rm = "%.2f" % (C10_rm / targetLength * 100) + "%"
    cover1_rm = "%.2f" % (C1_rm / targetLength * 100) + "%"
    return [uniformity, cover1, cover10, cover20, cover50, cover1_rm, cover10_rm, cover20_rm, cover50_rm]

# 性别校检
def genderAnalysis(genderFile):
    genderFileO = open(genderFile, "r", encoding="utf-8")
    SRY = 0
    for line in genderFileO:
        if line.replace("\n", "") != "":
            try:
                SRY = int(line)
            except:
                pass
    genderFileO.close()
    # 阈值
    if SRY <= 30:
        gender = "Female"
    elif SRY >= 200:
        gender = "Male"
    else:
        gender = "UNKNOWN"
    return gender

# fold 80
def fold80Analysis(fold80File):
    fold80FileO = open(fold80File, "r")
    fold80 = "-"
    check = False
    n = 0
    for line in fold80FileO:
        n += 1
        if n == 1:
            if "FOLD_80_BASE_PENALTY" in line:
                check = True
        if n == 2:
            if check:
                try:
                    fold80 = "%.2f" % float(line.replace("\n", ""))
                except:
                    pass
    fold80FileO.close()
    return fold80

# insert size
def insertSizeAnalysis(insertSizeFile):
    insertSizeFileO = open(insertSizeFile, "r")
    medianInsertSize = meanInsertSize = "-"
    check = False
    n = 0
    for line in insertSizeFileO:
        n += 1
        if n == 1:
            if "MEDIAN_INSERT_SIZE" in line:
                check = True
        if n == 2:
            if check:
                lines = line.replace("\n", "").split("\t")
                medianInsertSize = lines[0]
                meanInsertSize = "%.2f" % float(lines[1])
    insertSizeFileO.close()
    return [medianInsertSize, meanInsertSize]

# 输出结果
def main(sample, fastpJsonReport, bamdstReport, depthFile, genderFile, fold80File, insertSizeFile, outputResults):
    rawReads, rawBases, cleanReads, cleanBases, cleanQ20, cleanQ30, cleanGC, readLength = jsonAnalysis(fastpJsonReport)
    duplicationRate, mappedReads, mappedBases, mappedRate, targetReads, targetBases, targetRate, averageDepth, averageRmDepth = coverageAnalysis(bamdstReport, cleanBases, readLength)
    uniformity, cover1, cover10, cover20, cover50, cover1_rm, cover10_rm, cover20_rm, cover50_rm = uniformityAnalysis(depthFile, averageDepth)
    gender = genderAnalysis(genderFile)
    fold80 = fold80Analysis(fold80File)
    medianInsertSize, meanInsertSize = insertSizeAnalysis(insertSizeFile)
    outputFile = open(outputResults, "w")
    output = [
        sample, str(rawReads), str(rawBases),
        str(cleanReads), str(cleanBases), duplicationRate,
        cleanQ20, cleanQ30, cleanGC,
        str(mappedReads), str(mappedBases), mappedRate,
        str(targetReads), str(targetBases), targetRate,
        str(averageDepth), str(averageRmDepth), uniformity, fold80,
        cover1, cover10, cover20, cover50, cover1_rm, cover10_rm, cover20_rm, cover50_rm, medianInsertSize, meanInsertSize, gender
    ]
    outputFile.write("sampleID\trawReads\trawBases\tcleanReads\tcleanBases\tduplicatesRate\tcleanQ20\tcleanQ30\tcleanGC"\
        "\tmappedReads\tmappedBases\tmappedRate\ttargetReads\ttargetBases\ttargetRate\tdepth\tdepthRmdups\tUniformity\tFold80"\
        "\t1XCov\t10Xcov\t20XCov\t50XCov\t1XCov_rm\t10XCov_rm\t20XCov_rm\t50XCov_rm\tmedianInsertSize\tmeanInsertSize\tPredictGender\n")
    outputFile.write("\t".join(output) + "\n")
    outputFile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            qcstat germline pipe
        """,
        prog="qcstat_germline.py",
        usage="python3 qcstat_germline.py [-h] -i <sampleID> -o <outputFile> -j <fastp json report> -b <bamdst report> -d <bamdst depth file> -g <gender report> -f <fold 80> -s <insert size>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version",
        version="Version 0.2 20220629")
    parser.add_argument("-i", "--input", type=str,
        help="样本编号")
    parser.add_argument("-o", "--output", type=str,
        help="结果文件")
    parser.add_argument("-j", "--json", type=str,
        help="fastp json报告")
    parser.add_argument("-b", "--bamstat", type=str,
        help="bamdst报告")
    parser.add_argument("-d", "--depth", type=str,
        help="bamdst depth.tsv.gz")
    parser.add_argument("-g", "--gender", type=str,
        help="性别预测报告")
    parser.add_argument("-f", "--fold", type=str,
        help="fold80 报告")
    parser.add_argument("-s", "--size", type=str,
        help="insert size 报告")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(sample=args.input, fastpJsonReport=args.json, bamdstReport=args.bamstat, depthFile=args.depth, genderFile=args.gender, fold80File=args.fold, insertSizeFile=args.size, outputResults=args.output)

# end
