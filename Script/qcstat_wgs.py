# coding=utf-8
# pzw
# 20220928

"""
SampleID|rawReads|rawBases|cleanReads|cleanBases|duplicatesRate|cleanQ20|cleanQ30|cleanGC
|depth|depthRmdups|1XCov|10XCov|20XCov|30XCov|medianInsertSize|meanInsertSize|PredictGender
"""

import sys
import json
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

# dup解析
def dupAnalysis(dupFile):
    dupRate = "-"
    with open(dupFile, "r") as dups:
        line = dups.readlines()
        for l in range(len(line)):
            if line[l].startswith("LIBRARY"):
                dupRate = float(line[l+1].split("\t")[8])
    return dupRate

# genome信息
def genomecovAnalysis(genomecovFile):
    lines = []
    n = 1
    while n <= 30:
        lines.append("-")
        n += 1
    with open(genomecovFile, "r") as g:
        for line in g:
            if not line.startswith("1XCov"):
                if not line == "\n":
                    lines = line.replace("\n", "").split("\t")
    return lines  

# 性别校检
def genderAnalysis(genderFile, cutoff):
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
    if SRY < cutoff:
        gender = "Female"
    else:
        gender = "Male"
    return gender


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
def main(sample, fastpJsonReport, genderFile, dupFile, covFile, insertSizeFile, outputResults):
    rawReads, rawBases, cleanReads, cleanBases, cleanQ20, cleanQ30, cleanGC, readLength = jsonAnalysis(fastpJsonReport)
    dupRate = dupAnalysis(dupFile)
    duplicateRate = "%.2f" % (dupRate * 100) + "%"
    genomecovList = genomecovAnalysis(covFile)
    cov1X, cov10X, cov20X, cov30X, rmDepth, chrM, chrX, chrY = genomecovList[0:8]
    try:
        depth = "%.2f" % (float(rmDepth) / (1 - dupRate))
    except:
        depth = "-"
    gender = genderAnalysis(genderFile, 10)
    medianInsertSize, meanInsertSize = insertSizeAnalysis(insertSizeFile)

    # 输出
    outputFile = open(outputResults, "w")
    output = [
        sample, str(rawReads), str(rawBases),
        str(cleanReads), str(cleanBases), duplicateRate,
        cleanQ20, cleanQ30, cleanGC,
        depth, rmDepth, cov1X, cov10X, cov20X, cov30X, medianInsertSize, meanInsertSize, gender,
        chrM, chrX, chrY
    ]
    chromHeader = []
    for i in range(1, 23):
        output.append(genomecovList[i+7])
        chromHeader.append("chr" + str(i))
    outputString = "sampleID\trawReads\trawBases\tcleanReads\tcleanBases\tduplicatesRate\tcleanQ20\tcleanQ30\tcleanGC"\
        "\tdepth\tdepthRmdups\t1XCov\t10Xcov\t20XCov\t30XCov\tmedianInsertSize\tmeanInsertSize\tPredictGender\tchrM\tchrX\tchrY"
    outputString = outputString + "\t" + "\t".join(chromHeader)
    outputFile.write(outputString + "\n")
    outputFile.write("\t".join(output) + "\n")
    outputFile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            qcstat wgs pipe
        """,
        prog="qcstat_wgs.py",
        usage="python3 qcstat_wgs.py [-h] -i <sampleID> -o <outputFile> -j <fastp json report> -c <cov report> -d <dup report> -g <gender report> -s <insert size>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version",
        version="Version 0.1 20220928")
    parser.add_argument("-i", "--input", type=str,
        help="样本编号")
    parser.add_argument("-o", "--output", type=str,
        help="结果文件")
    parser.add_argument("-j", "--json", type=str,
        help="fastp json报告")
    parser.add_argument("-c", "--cov", type=str,
        help="cov report")
    parser.add_argument("-d", "--dup", type=str,
        help="picard/gatk dups metrix")
    parser.add_argument("-g", "--gender", type=str,
        help="性别预测报告")
    parser.add_argument("-s", "--size", type=str,
        help="insert size 报告")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(sample=args.input, fastpJsonReport=args.json, dupFile=args.dup, covFile=args.cov, genderFile=args.gender, insertSizeFile=args.size, outputResults=args.output)

# end
