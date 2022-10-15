# coding=utf-8
# pzw
# 20221011
# 归档CNVseq分析结果

import os
import sys
import argparse

# 转存样本
def restore(sample, id, direction):
    cromwellDir = "/home/novelbio/project/cromwell/cromwell-executions"
    project = "CNVSeqPipe"
    idPath = cromwellDir + "/" + project + "/" + id
    outDir = direction + "/" + sample
    
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("成功创建文件夹", outDir)
    else:
        print("文件夹已存在")

    if project == "CNVSeqPipe":
        cmd = """
            cp {idPath}/call-YKQC/execution/{sample}.QC.txt {outDir}
            cp {idPath}/call-RemoveDup/execution/{sample}.rmdup.bam* {outDir}
            cp {idPath}/call-Freec/execution/{sample}.ratio.bed {outDir}
            cp {idPath}/call-Freec/execution/{sample}.ratio.circos {outDir}
            cp {idPath}/call-Freec/execution/{sample}.ratio.txt {outDir}
            cp {idPath}/call-Freec/execution/{sample}.CNVs.txt {outDir}
            cp {idPath}/call-FreecPlot/execution/{sample}.cnv.p.value.txt {outDir}
            cp {idPath}/call-FreecPlot/execution/{sample}.freec.log2.png {outDir}
            cp {idPath}/call-FreecPlot/execution/{sample}.freec.png {outDir}
            cp {idPath}/call-AnnotSV/execution/{sample}.annotsv.txt {outDir}
            cp -r {idPath}/call-FreecPlot/execution/{sample}_chrom_pics {outDir}
        """.format(idPath=idPath, sample=sample, outDir=outDir)
        os.system(cmd)
        print(sample, "完成数据转移")
    
    else:
        print("未找到项目名", project)
        exit()


# 主流程
def main(group, outputDir):

    tmpdir = "/home/novelbio/pipeline/WESpipeWDL/Config"

    cmd = """
        /home/novelbio/.local/bin/oliver st -d -g {group} --grid-style xx | tr -s ' ' '\t' > {tmpdir}/tmp.txt
    """.format(group=group, tmpdir=tmpdir)
    os.system(cmd)

    tmpFile = open(tmpdir + "/tmp.txt", "r", encoding="utf-8")
    for line in tmpFile:
        if line.startswith("\tJob\tName"):
            continue
        elif line.startswith("Job\tName"):
            continue
        elif line.startswith("-------"):
            continue
        else:
            lines = line.split("\t")
            status = lines[4]
            jobName = lines[0]
            id = lines[2]

            if jobName == "":
                status = lines[5]
                jobName = lines[1]
                id = lines[3]

            if status != "Succeeded":
                print(jobName, id, status, "\tNot Succeeded, skip")
            else:
                restore(jobName, id, outputDir)
    tmpFile.close()


# 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="RestoreJobGroup",
        prog="RestoreJobGroup.py",
        usage="python3 RestoreJobGroup.py [-h] -g <group name> -o <outputDir>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version",
        version="Version 0.2 20220808")
    parser.add_argument("-g", "--group", type=str,
        help="样本批次")
    parser.add_argument("-o", "--output", type=str,
        help="输出路径")

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(group=args.group, outputDir=args.output)
