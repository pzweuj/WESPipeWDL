# coding=utf-8
# pzw
# 20220502

"""
适用于snpeff+annovar注释后下游的注释整理流程
仅用于线粒体

java -jar snpeff.jar \\
    -c snpeff.config \\
    MT raw.vcf > snpeff.vcf

convert2annovar.pl -format vcf4 \\
    -allsample -withfreq \\
    snpeff.vcf --includeinfo > avinput

table_annovar.pl avinput \\
    /home/novelbio/databases/humandb \\
    -buildver hg19 \\
    -out anno -remove \\
    -protocol refGene,avsnp150,clinvar_20220320,mitomapB37 \\
    -operation g,f,f,f \\
    -nastring - -thread 8 -otherinfo
"""

import sys
import argparse


# vcf基因型可读
def checkGTPercent(percent):
    if percent >= 0.97:
        gt = "Hom"
    elif percent >= 0.03:
        gt = "Het"
    elif percent >= 0:
        gt = "Wt"
    else:
        gt = "NA"
    return gt

# 拆分不同的基因
def SplitLine(line, annoRank):
    lines = line.split("\t")
    annotation = lines[annoRank]
    annotations = annotation.split(";")
    snpeff = "-"
    lineList = [line]
    ANNRank = 0
    for anno in range(len(annotations)):
        if annotations[anno].startswith("ANN="):
            snpeff = annotations[anno]
            ANNRank = anno

    if snpeff != "-":
        geneDict = {}
        snpeffs = snpeff.replace("ANN=", "").split(",")
        for s in snpeffs:
            sSplit = s.split("|")
            gene = sSplit[3]
            try:
                geneDict[gene].append(s)
            except:
                geneDict[gene] = []
                geneDict[gene].append(s)
        if len(geneDict.keys()) > 1:
            lineList = []
            for k in geneDict.keys():
                annotations[ANNRank] = "ANN=" + ",".join(geneDict[k])
                annotation = ";".join(annotations)
                lines[annoRank] = annotation
                lineList.append("\t".join(lines))
    
    return lineList

# 坐标处理
def FixEndLocation(end, ref, alt):
    endFix = end
    if ref == "-":
        t = "Insertion"
        endFix = "-"
    elif alt == "-":
        t = "Deletion"
    elif len(ref) == 1 and len(alt) == 1:
        t = "SNV"
    else:
        t = "Complex"
    return [endFix, t]

# 氨基酸简写
def AATranslate(seq):
    transDict = {
        "Ala": "A", "Arg": "R", "Asn": "N",
        "Asp": "D", "Cys": "C", "Gln": "Q",
        "Glu": "E", "Gly": "G", "His": "H",
        "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W",
        "Tyr": "Y", "Val": "V"
    }
    for i in transDict.keys():
        if i in seq:
            seq = seq.replace(i, transDict[i])
    
    return seq

# 主流程
def main(annovarResultsFile, resultsFile, trioName, silent):
    annovarResults = open(annovarResultsFile, "r")
    results = open(resultsFile, "w")

    lineDict = {}
    lineDictList = []
    for line in annovarResults:
        if line.startswith("Chr\tStart"):
            lines = line.replace("\n", "").split("\t")
            for k in lines:
                lineDict[k] = "-"
                lineDictList.append(k)
    annovarResults.close()

    annovarResults = open(annovarResultsFile, "r")
    headCount = 0
    for line in annovarResults:
        if not line.startswith("Chr\tStart"):
            lineSplitList = SplitLine(line, 30)
            for l in lineSplitList:
                lines = l.replace("\n", "").split("\t")
                for i in range(len(lines)):
                    lineDict[lineDictList[i]] = lines[i]

                # 开始处理
                fe = FixEndLocation(lineDict["End"], lineDict["Ref"], lineDict["Alt"])
                lineDict["Type"] = fe[1]
                lineDict["End"] = fe[0]

                # 处理基因
                otif11 = lineDict["Otherinfo11"].split(";")
                snpeff_anno = "ANN=?|N|N|N|N|transcript|N|N|N|N|N|N|N|N|N|"
                for o11 in otif11:
                    if o11.startswith("ANN="):
                        snpeff_anno = o11
                snpeff_annos = snpeff_anno.replace("ANN=", "").split(",")
                
                # 线粒体默认取第一个
                sa = snpeff_annos[0]
                s = sa.split("|")
                for ss in range(len(s)):
                    if s[ss] == "":
                        s[ss] = "-"
                
                lineDict["Gene"] = s[3].replace("-CHR_END", "")
                lineDict["Feature"] = s[4].replace("-CHR_END", "")
                lineDict["cHGVS"] = s[9].replace("n.", "m.")
                lineDict["pHGVS"] = s[10]

                # 处理VAF
                trioNames = trioName.split(",")
                otif12 = lineDict["Otherinfo12"].split(":")
                tnNum = 1

                # 处理VAF
                trioNames = trioName.split(",")
                otif12 = lineDict["Otherinfo12"].split(":")
                tnNum = 1

                try:
                    DPList = []
                    for tn in trioNames:
                        otifCheck = lineDict["Otherinfo" + str(12 + tnNum)].split(":")
                        checkBand = lineDict["Otherinfo" + str(12 + tnNum)]
                        checkGTs = "NA"
                        tnNum += 1
                        format_zip = list(zip(otif12, otifCheck))
                        format_dict = {}
                        for f in format_zip:
                            format_dict[f[0]] = f[1]
                        DPList.append(int(format_dict["DP"]))
                        lineDict[tn + "_DP"] = format_dict["DP"]
                        lineDict[tn + "_AltAD"] = format_dict["AD"]
                        if not "," in lineDict[tn + "_AltAD"]:
                            lineDict[tn + "_AltAD"] = lineDict[tn + "_AltAD"]
                        else:
                            lineDict[tn + "_AltAD"] = lineDict[tn + "_AltAD"].split(",")[1]
                        lineDict[tn + "_RefAD"] = str(int(lineDict[tn + "_DP"]) - int(lineDict[tn + "_AltAD"]))

                        try:
                            AF = "%.2f" % ((float(lineDict[tn + "_AltAD"]) / float(lineDict[tn + "_DP"])) * 100) + "%"
                            checkGTs = checkGTPercent(float(lineDict[tn + "_AltAD"]) / float(lineDict[tn + "_DP"]))
                        except Exception:
                            AF = "-"
                        lineDict[tn + "_VAF"] = AF

                        lineDict[tn + "_INFO"] = checkBand
                        lineDict[tn + "_GT"] = checkGTs

                except:
                    continue
                
                # 列名
                if headCount == 0:
                    header = ["Chr", "Start", "End", "Ref", "Alt",
                        "Gene", "Feature", "cHGVS", "pHGVS", "avsnp150",
                        "CLNDN", "CLNSIG", "MitoMap_Disease",
                        "MitoMap_Quartile", "MitoMap_PMID", "MitoMap_AF",
                        "MAPQ", "VcfInfo"
                    ]
                    for tn in trioNames:
                        header.append(tn + "_DP")
                        header.append(tn + "_RefAD")
                        header.append(tn + "_AltAD")
                        header.append(tn + "_VAF")

                    for tn in trioNames:
                        header.append(tn + "_INFO")

                    for tn in trioNames:
                        header.append(tn + "_GT")

                    results.write("\t".join(header) + "\n")
                    headCount = 1

                output = [lineDict["Chr"], lineDict["Start"], lineDict["End"], lineDict["Ref"],
                    lineDict["Alt"], lineDict["Gene"], lineDict["Feature"], lineDict["cHGVS"],
                    AATranslate(lineDict["pHGVS"]), lineDict["avsnp150"], 
                    lineDict["CLNDN"], lineDict["CLNSIG"],
                    lineDict["MitoMap_Disease"], lineDict["MitoMap_Quartile"], lineDict["MitoMap_PMID"], lineDict["MitoMap_AF"], lineDict["Otherinfo3"], lineDict["Otherinfo11"]
                ]

                for tn in trioNames:
                    output.append(lineDict[tn + "_DP"])
                    output.append(lineDict[tn + "_RefAD"])
                    output.append(lineDict[tn + "_AltAD"])
                    output.append(lineDict[tn + "_VAF"])

                output.append(lineDict["Otherinfo12"])

                for tn in trioNames:
                    output.append(lineDict[tn + "_INFO"])

                for tn in trioNames:
                    output.append(lineDict[tn + "_GT"])

                if silent != "True":
                    print("\t".join(output))
                results.write("\t".join(output) + "\n")

    results.close()
    annovarResults.close()                     

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            AnnoMT

            适用于snpeff+annovar注释后结果整理。
            其中annovar需要注释数据库为
                    refGene,avsnp150,clinvar_20220320,mitomapB37
        
        """,
        prog="AnnoMT.py",
        usage="python3 AnnoMT.py [-h] -i <annoFile> -o <annoFixFile> -t <trioNameList>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version",
        version="Version 0.1 20220502")
    parser.add_argument("-i", "--anno", type=str,
        help="注释结果文件")
    parser.add_argument("-o", "--output", type=str,
        help="注释结果整理")
    parser.add_argument("-t", "--trio", type=str,
        help="受检者名称，按照vcf结果中名称顺序填写，以‘,’分隔")
    parser.add_argument("-s", "--silent", type=str,
        help="是否静默运行，默认是", default="True")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(annovarResultsFile=args.anno, resultsFile=args.output, trioName=args.trio, silent=args.silent)

# end
        
