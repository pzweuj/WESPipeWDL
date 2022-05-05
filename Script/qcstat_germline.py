# coding=utf-8
# pzw
# 20220426

"""
SampleID|rawReads|rawBases|cleanReads|cleanBases|duplicatesRate|cleanQ20|cleanQ30|cleanGC
|mappedReads|mappedBases|mappedRate|targetReads|targetBases|targetRate|depth|depthRmdups|Uniformity
|1XCov|10XCov|20XCov|50XCov
"""

import sys
import json
import gzip

sample = sys.argv[1]
fastpJsonReport = sys.argv[2]
bamdstReport = sys.argv[3]
depthFile = sys.argv[4]
genderFile = sys.argv[5]
outputResults = sys.argv[6]
outputFile = open(outputResults, "w")


# json解析
# fastpJsonReport = "202111.json"
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

# bamdst解析
# bamdstReport = "coverage.report"
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

# 以reads计算比对率及捕获率
# mappedRate = "%.2f" % (mappedReads / cleanReads * 100) + "%"
# targetRate = "%.2f" % (targetReads / mappedReads * 100) + "%"

# 深度
# depthFile = "depth.tsv.gz"
depthReport = gzip.open(depthFile, "rt")
targetLength = 0
uniform = 0
C1 = 0
C10 = 0
C20 = 0
C50 = 0
targetDepthBases = 0
targetRmDepthBases = 0
for line in depthReport:
    if not line.startswith("#"):
        targetLength += 1
        lines = line.split("\t")
        chrom = lines[0]
        pos = lines[1]
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
uniformity = "%.2f" % (uniform / targetLength * 100) + "%"
cover50 = "%.2f" % (C50 / targetLength * 100) + "%"
cover20 = "%.2f" % (C20 / targetLength * 100) + "%"
cover10 = "%.2f" % (C10 / targetLength * 100) + "%"
cover1 = "%.2f" % (C1 / targetLength * 100) + "%"

# 性别校检
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

# 输出结果
output = [
    sample, str(rawReads), str(rawBases),
    str(cleanReads), str(cleanBases), duplicationRate,
    cleanQ20, cleanQ30, cleanGC,
    str(mappedReads), str(mappedBases), mappedRate,
    str(targetReads), str(targetBases), targetRate,
    str(averageDepth), str(averageRmDepth), uniformity,
    cover1, cover10, cover20, cover50, gender
]

outputFile.write("sampleID\trawReads\trawBases\tcleanReads\tcleanBases\tduplicatesRate\tcleanQ20\tcleanQ30\tcleanGC"\
    "\tmappedReads\tmappedBases\tmappedRate\ttargetReads\ttargetBases\ttargetRate\tdepth\tdepthRmdups\tUniformity"\
    "\t1XCov\t10Xcov\t20XCov\t50XCov\tPredictGender\n")
outputFile.write("\t".join(output) + "\n")
outputFile.close()

