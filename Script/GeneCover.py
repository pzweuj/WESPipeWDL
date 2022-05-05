# coding=utf-8
# pzw
# 20220505
# 统计每个基因的覆盖度及平均深度

"""
环境中需有zcat、bedtools

适用于
bamdst -p gene.bed input.bam -o output_dir
生成的region.tsv.gz文件

其中bed文件需求格式为
chrom\tstart\tend\tgeneSymbol_x
参考/slurm/databases/b37/bed/IDTxGenExome2_exonsCDSmerged_b37.bed
"""

import os
import sys
import shutil

# 获得结果
def ExonCover(regionGz, bed, geneCoverFile):
    cmd = """
        zcat {regionGz} \\
            | bedtools intersect -a {bed} -b - -wb \\
            | awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$10"\\t"$13}}' \\
            > {geneCover}
    """.format(regionGz=regionGz, bed=bed, geneCover=geneCoverFile)
    os.system(cmd)

# 处理结果文件
def GeneCover(geneCoverFile):
    geneCover = open(geneCoverFile, "r", encoding="utf-8")
    geneDict = {}
    gene = "-"
    geneLength = 0
    sumDepth = 0
    covLength = 0
    for line in geneCover:
        lines = line.replace("\n", "").split("\t")
        geneSymbol = lines[3].split("_")[0]
        start = int(lines[1])
        end = int(lines[2])
        exonLength = end - start + 1
        avgDepth = float(lines[4])
        coverage = float(lines[5]) / 100
        if geneSymbol != gene:
            try:
                avgDepthGene = "%.2f" % (sumDepth / geneLength)
                covGene = "%.2f" % (covLength / geneLength * 100) + "%"
                geneDict[gene] = [avgDepthGene, covGene]
            except:
                pass
            geneLength = 0
            sumDepth = 0
            covLength = 0
            gene = geneSymbol
        geneLength += exonLength
        sumDepth += (exonLength * avgDepth)
        covLength += (exonLength * coverage)
    return geneDict

def main(regionGz, bed, output):
    if not os.path.exists("tmpDir"):
        os.makedirs("tmpDir")
    ExonCover(regionGz, bed, "tmpDir/tmp.txt")
    geneDict = GeneCover("tmpDir/tmp.txt")
    op = open(output, "w", encoding="utf-8")
    op.write("#Gene\tAvgDepth\tCoverage%\n")
    for k in geneDict.keys():
        if k != "-":
            op.write("\t".join([k, geneDict[k][0], geneDict[k][1]]) + "\n")
    op.close()
    shutil.rmtree("tmpDir")

####################################################
try:
    region = sys.argv[1]
    bedFile = sys.argv[2]
    outputFile = sys.argv[3]
    main(region, bedFile, outputFile)
except:
    print("Usage: python3 GeneCover.py <region.tsv.gz> <bed> <output>")

# end
