# coding=utf-8
# pzw
# 20220408

import os
import sys

def getAbsPath():
    now = os.path.abspath(os.path.dirname(sys.argv[0]))
    return now

try:
    clinvarZip = sys.argv[1]
    output = sys.argv[2]
except:
    print("python3 GetPathogenicHotRegion.py <clinvar.vcf.gz> <output.bed>")
    exit()

# 获得Pathogenic或Likely_pathogenic
tmpDir = getAbsPath() + "/tmp"
if not os.path.exists(tmpDir):
    os.makedirs(tmpDir)
cmd = """
        zcat {clinvarZip} | grep 'athogenic' | \\
            grep -v Conflicting_interpretations_of_pathogenicity | \\
            grep -v 'CLNSIG=Benign' | \\
            grep -v 'CLNSIG=Likely_benign' | \\
            grep missense_variant > {tmpDir}/tmp.vcf
        zcat {clinvarZip} | grep '^#' > {tmpDir}/header
        cat {tmpDir}/header {tmpDir}/tmp.vcf > {tmpDir}/clinvar.vcf
    """.format(clinvarZip=clinvarZip, tmpDir=tmpDir)
os.system(cmd)

# 获得bed文件
clinvarO = open(tmpDir + "/clinvar.vcf", "r", encoding="utf-8")
clinvarHotRegion = open(tmpDir + "/clinvar.hotregion.bed", "w", encoding="utf-8")
for line in clinvarO:
    if not line.startswith("#"):
        lines = line.split("\t")
        chrom = lines[0]
        start = lines[1]
        ref = lines[3]
        end_r = str(int(start) + len(ref) - 1 + 15)
        start_r = str(int(start) - 15)
        clinvarHotRegion.write("\t".join([chrom, start_r, end_r]) + "\n")
clinvarHotRegion.close()
clinvarO.close()
cmd = """
    bedtools intersect -a {tmpDir}/clinvar.hotregion.bed \\
        -b {tmpDir}/clinvar.vcf -c | awk '{{if($4>=3) print $0}}' \\
        | bedtools merge -i - > {output}
    cat {output} | awk '{{m+=$3-$2}}END{{print m}}'
    rm -rf {tmpDir}
    mv {output} {output}.tmp
""".format(tmpDir=tmpDir, output=output)
os.system(cmd)

# 染色体编号调整
beforeFix = open(output + ".tmp", "r", encoding="utf-8")
afterFix = open(output, "w", encoding="utf-8")
for line in beforeFix:
    afterFix.write("chr" + line.replace("MT", "M").replace("\n", "\tY\n"))
afterFix.close()
beforeFix.close()
os.system("rm {}.tmp".format(output))

