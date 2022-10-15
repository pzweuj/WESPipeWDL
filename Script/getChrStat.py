# coding=utf-8
# pzw
# 2020928
# 统计bedtools genomecov生产的文件



"""
bedtools genomecov -ibam test.bam -bga | gzip - > cov.txt.gz
"""
import sys
import gzip


sample = gzip.open(sys.argv[1], "rt")
output = open(sys.argv[2], "w")

# header
header = ["1XCov", "10XCov", "20XCov", "30XCov", "rmDup_depth", "chrM_depth", "chrX_depth", "chrY_depth"]
for i in range(1, 23):
    header.append("chr" + str(i) + "_depth")
output.write("\t".join(header) + "\n")

# 构建关注染色体列表
chromList = ["X", "Y", "MT"]
for i in range(1, 23):
    chromList.append(str(i))
chromDict = {}
for c in chromList:
    chromDict[c] = [0, 0]

# 计算每个染色体的长度及覆盖碱基
cov1X = cov10X = cov20X = cov30X = 0

for line in sample:
    lines = line.replace("\n", "").split("\t")
    chrom = lines[0]
    start = lines[1]
    end = lines[2]
    covBase = int(lines[3])
    regionLen = int(end) - int(start)

    # 30X Cov
    if covBase >= 30:
        cov1X += regionLen
        cov10X += regionLen
        cov20X += regionLen
        cov30X += regionLen
    # 20X Cov
    elif covBase >= 20:
        cov1X += regionLen
        cov10X += regionLen
        cov20X += regionLen
    # 10X Cov
    elif covBase >= 10:
        cov1X += regionLen
        cov10X += regionLen
    elif covBase >= 1:
        cov1X += regionLen
    else:
        pass

    # 各个染色体
    if chrom in chromList:
        chromDict[chrom][0] += regionLen
        chromDict[chrom][1] += (regionLen * covBase)

# 计算基因组覆盖度
genomeLen = 0
sumBases = 0
for k in chromDict:
    genomeLen += chromDict[k][0]
    sumBases += chromDict[k][1]
cov1XPercent = "%.2f" % (float(cov1X) / genomeLen * 100) + "%"
cov10XPercent = "%.2f" % (float(cov10X) / genomeLen * 100) + "%"
cov20XPercent = "%.2f" % (float(cov20X) / genomeLen * 100) + "%"
cov30XPercent = "%.2f" % (float(cov30X) / genomeLen * 100) + "%"
rmDupDepth = "%.2f" % (float(sumBases) / genomeLen)
outputStringList = [cov1XPercent, cov10XPercent, cov20XPercent, cov30XPercent, rmDupDepth]

# 计算各个染色体平均深度
depthDict = {}
for k in chromDict:
    depthDict[k] = "%.2f" % ((float(chromDict[k][1]) / chromDict[k][0]))

outputStringList.append(depthDict["MT"])
outputStringList.append(depthDict["X"])
outputStringList.append(depthDict["Y"])
for i in range(1, 23):
    outputStringList.append(depthDict[str(i)])

output.write("\t".join(outputStringList) + "\n")
output.close()

print(cov1XPercent, cov10XPercent, cov20XPercent, cov30XPercent)
print(depthDict)

