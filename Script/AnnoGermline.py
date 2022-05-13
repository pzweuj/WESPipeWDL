# coding=utf-8
# pzw
# 20220502

"""
适用于snpeff+annovar注释后下游的注释整理流程
用于HC与Freebayes结果

java -jar snpeff.jar \\
    -c snpeff.config \\
    hg19 raw.vcf > snpeff.vcf

convert2annovar.pl -format vcf4 \\
    -allsample -withfreq \\
    snpeff.vcf --includeinfo > avinput

table_annovar.pl avinput \\
    /home/novelbio/databases/humandb \\
    -buildver hg19 \\
    -out anno -remove \\
    -protocol refGene,cytoBand,avsnp150,gnomad211_genome,gnomad211_exome,1000g2015aug_all,
        1000g2015aug_eas,exac03,esp6500siv2_all,clinvar_20220320,dbnsfp42a,dbscsnv11,
        intervar_20180118,SpliceAI \\
    -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f \\
    -nastring - -thread 8 -otherinfo
"""

import os
import sys
import argparse


# vcf基因型可读
def checkGTPercent(percent):
    if percent >= 0.95:
        gt = "Hom"
    elif percent >= 0.05:
        gt = "Het"
    elif percent >= 0:
        gt = "Wt"
    else:
        gt = "NA"
    return gt

# 获得参考转录本字典
def refTranscript(ts):
    refT = open(ts, "r")
    refDict = {}
    for r in refT:
        if not r.startswith("#"):
            gene = r.split("\t")[0]
            ts = r.split("\t")[1]
            refDict[gene] = ts
    refT.close()
    return refDict

# 补充OMIM注释
def appendOMIM(input, output, omim):
    inputO = open(input, "r", encoding="utf-8")
    outputO = open(output, "w", encoding="utf-8")
    omim = open(omim, "r", encoding="utf-8")
    geneDict = {}
    for line in omim:
        if not line.startswith("#"):
            lines = line.replace("\n", "").split("\t")
            gene = lines[0]
            disease = lines[1]
            phenotype = lines[2]
            geneDict[gene] = [disease, phenotype]
    omim.close()
    for line in inputO:
        lines = line.split("\t")
        gene = lines[5]
        if line.startswith("Chr\t"):
            lines.insert(11, "Phenotype")
            lines.insert(11, "Disease")
        else:
            try:
                lines.insert(11, geneDict[gene][1])
                lines.insert(11, geneDict[gene][0])
            except:
                lines.insert(11, "-")
                lines.insert(11, "-")
        outputO.write("\t".join(lines))
    inputO.close()
    outputO.close()

# 判断人群频率是否均低于某值
def checkAF(afList, freq):
    n = 0
    for a in afList:
        if a == "." or a == "-":
            af = 0
        else:
            af = float(a)
        
        if af > freq:
            n += 1
    if n > 0:
        return "F"
    else:
        return "T"

# InterVar合并
def MergeInterVar(analysisDict):
    analysisDict["InterVar_sig"] = "PVS1=" + analysisDict["PVS1"] + ";" + \
        "PS1=" + analysisDict["PS1"] + ";" + \
        "PS2=" + analysisDict["PS2"] + ";" + \
        "PS3=" + analysisDict["PS3"] + ";" + \
        "PS4=" + analysisDict["PS4"] + ";" + \
        "PM1=" + analysisDict["PM1"] + ";" + \
        "PM2=" + analysisDict["PM2"] + ";" + \
        "PM3=" + analysisDict["PM3"] + ";" + \
        "PM4=" + analysisDict["PM4"] + ";" + \
        "PM5=" + analysisDict["PM5"] + ";" + \
        "PM6=" + analysisDict["PM6"] + ";" + \
        "PP1=" + analysisDict["PP1"] + ";" + \
        "PP2=" + analysisDict["PP2"] + ";" + \
        "PP3=" + analysisDict["PP3"] + ";" + \
        "PP4=" + analysisDict["PP4"] + ";" + \
        "PP5=" + analysisDict["PP5"] + ";" + \
        "BA1=" + analysisDict["BA1"] + ";" + \
        "BS1=" + analysisDict["BS1"] + ";" + \
        "BS2=" + analysisDict["BS2"] + ";" + \
        "BS3=" + analysisDict["BS3"] + ";" + \
        "BS4=" + analysisDict["BS4"] + ";" + \
        "BP1=" + analysisDict["BP1"] + ";" + \
        "BP2=" + analysisDict["BP2"] + ";" + \
        "BP3=" + analysisDict["BP3"] + ";" + \
        "BP4=" + analysisDict["BP4"] + ";" + \
        "BP5=" + analysisDict["BP5"] + ";" + \
        "BP6=" + analysisDict["BP6"] + ";" + \
        "BP7=" + analysisDict["BP7"]
    return analysisDict

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

# 根据参考转录本选择需解析的snpeff结果
def SnpeffRefTranscriptSelect(snpeffAnnotation, refTranscriptFile):
    snpeffAnnos = snpeffAnnotation.replace("ANN=", "").split(",")
    snpeffSelect = snpeffAnnos[0]
    if len(snpeffAnnos) > 1:
        gene = snpeffAnnos[0].split("|")[3]
        refTrans = refTranscript(refTranscriptFile)

        try:
            ts = refTrans[gene]
            for sa in snpeffAnnos:
                if ts.split(".")[0] in sa:
                    if ts.split("|")[6] != "":
                        snpeffSelect = sa
        except:
            pass
    return snpeffSelect

# 坐标处理并获得Type
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

# 同义突变过滤
# 过滤人群频率＞5% 或
# clinvar为B/LB 或
# 测序深度小于10 的同义突变
def SynonymousFilter(DP, AF, EAS, clinvar):
    if AF == "-":
        AF = "0"
    if EAS == "-":
        EAS = "0"
    
    if "enign" in clinvar:
        f = True
    elif DP < 10:
        f = True
    elif float(AF) > 0.05:
        f = True
    elif float(EAS) > 0.05:
        f = True
    else:
        f = False
    return f

# PP2
# 补充一个先检的PP2基因列表
def appendPP2Gene(input, output, PP2Data, omim, geneCov):
    PP2List = []
    rank = 44
    if omim != "True" or omim != "true":
        rank = 42
        if geneCov != "False":
            rank = 44
    else:
        if geneCov != "False":
            rank = 46

    with open(PP2Data, "r", encoding="utf-8") as p:
        for line in p:
            PP2List.append(line.replace("\n", ""))
    inputO = open(input, "r", encoding="utf-8")
    outputO = open(output, "w", encoding="utf-8")
    for line in inputO:
        lines = line.split("\t")
        gene = lines[5]
        if line.startswith("Chr\t"):
            lines.insert(rank, "PP2Gene")
        else:
            if gene in PP2List:
                lines.insert(rank, "Y")
            else:
                lines.insert(rank, "-")
        outputO.write("\t".join(lines))
    outputO.close()
    inputO.close()

# Consequence
def Consequence(analysisDict, snpeffAnnotation):
    # Complex的情况，仍需要判断氨基酸突变
    # NCCL2022 将Complex_mutation删除了
    if "missense" in snpeffAnnotation:
        analysisDict["Consequence"] = "Missense_substitution"
    elif "splice" in snpeffAnnotation:
        analysisDict["Consequence"] = "Splice_Site_mutation"
    elif "synonymous" in snpeffAnnotation:
        analysisDict["Consequence"] = "Synonymous_substitution"
    elif "inframe_deletion" in snpeffAnnotation:
        analysisDict["Consequence"] = "Inframe_deletion"
    elif "inframe_insertion" in snpeffAnnotation:
        analysisDict["Consequence"] = "Inframe_insertion"
    elif "frameshift" in snpeffAnnotation:
        if analysisDict["Type"] == "Insertion":
            analysisDict["Consequence"] = "Frameshift_insertion"
        elif analysisDict["Type"] == "Deletion":
            analysisDict["Consequence"] = "Frameshift_deletion"
        else:
            # NCCL2022没有这个判断
            analysisDict["Consequence"] = "Frameshift_variant"
    elif "stop_gained" in snpeffAnnotation:
        analysisDict["Consequence"] = "Nonsense_substitution"
        if analysisDict["Type"] == "Deletion" or analysisDict["Type"] == "Insertion":
            analysisDict["Consequence"] = "Truncation_mutation"
    elif "stop_lost" in snpeffAnnotation:
        analysisDict["Consequence"] = "Elongation_mutation"
    else:
        analysisDict["Consequence"] = "Other"
    return analysisDict

# 获得基因覆盖度
def GeneCov(geneCoverFile):
    geneCover = open(geneCoverFile, "r", encoding="utf-8")
    geneCoverDict = {}
    for line in geneCover:
        if not line.startswith("#"):
            lines = line.replace("\n", "").split("\t")
            gene = lines[0]
            depth = lines[1]
            coverage = lines[2]
            geneCoverDict[gene] = [depth, coverage]
    geneCover.close()
    return geneCoverDict

# 主流程
def main(annovarResultsFile, resultsFile, refTranscriptFile, trioName, silent, omim, syno, pp2, geneCoverage, lite, filterSymbol):
    annovarResults = open(annovarResultsFile, "r")
    results = open(resultsFile, "w")

    # 根据标题生成字典
    lineDict = {}
    lineDictList = []
    for line in annovarResults:
        if line.startswith("Chr\tStart"):
            lines = line.replace("\n", "").split("\t")
            for k in lines:
                lineDict[k] = "-"
                lineDictList.append(k)
    annovarResults.close()

    # 再次打开开始分析
    annovarResults = open(annovarResultsFile, "r")
    headCount = 0
    for line in annovarResults:
        if not line.startswith("Chr\tStart"):
            # 切割snpeff结果，这里注释结果在第212列
            lineSplitList = SplitLine(line, 211)
            
            # lite模式，不进行基因分离
            if lite == "True" or lite == "true":
                lineSplitList = [line]
            
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

                # 使用参考转录本
                snpeff_select = SnpeffRefTranscriptSelect(snpeff_anno, refTranscriptFile)         
                s = snpeff_select.split("|")
                for ss in range(len(s)):
                    if s[ss] == "":
                        s[ss] = "-"

                lineDict["SEAnnotation"] = s[1]
                lineDict["Gene"] = s[3].replace("-CHR_END", "")

                # 过滤MIR/LOC/LINC系列
                if filterSymbol == "True" or filterSymbol == "true":
                    if lineDict["Gene"].startswith("MIR"):
                        continue
                    elif lineDict["Gene"].startswith("LOC"):
                        continue
                    elif lineDict["Gene"].startswith("LINC"):
                        continue
                    else:
                        pass

                lineDict["Feature"] = s[4].replace("-CHR_END", "")
                lineDict["Transcript"] = s[6]
                lineDict["AffectedExon"] = s[8]
                lineDict["cHGVS"] = s[9]
                lineDict["pHGVS"] = s[10]

                # 处理Consequence
                lineDict = Consequence(lineDict, lineDict["SEAnnotation"])

                # 处理ClinvarPathogenic
                if "ClinvarPath" in lineDict.keys():
                    if int(lineDict["ClinvarPath"]) > 0:
                        lineDict["ClinvarPath"] = "Y"
                    else:
                        lineDict["ClinvarPath"] = "N"
                else:
                    lineDict["ClinvarPath"] = "-"

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

                    # 过滤同义突变
                    if syno == "True" or syno == "true":
                        if lineDict["Consequence"] == "Synonymous_substitution":
                            if SynonymousFilter(max(DPList), lineDict["AF_exome"], lineDict["AF_exome_eas"], lineDict["CLNSIG"]):
                                continue

                except:
                    continue

                # 基因覆盖度与深度
                if os.path.exists(geneCoverage):
                    geneCoverageDict = GeneCov(geneCoverage)

                # 列名
                if headCount == 0:
                    header = [
                        "Chr", "Start", "End", "Ref", "Alt", "Gene", "Type", "Transcript", "cHGVS", "pHGVS", "Consequence", "AffectedExon", "cytoBand", "avsnp150",
                        "Lower0.01", "Lower0.05", "AF", "AF_popmax", "AF_male", "AF_female", "AF_raw", "AF_eas", "AF_afr", "AF_sas", "AF_amr", "AF_nfe", "AF_fin", "AF_asj", "AF_oth", "non_topmed_AF_popmax", "non_neuro_AF_popmax", "non_cancer_AF_popmax", "controls_AF_popmax", "AF_exome", "AF_exome_eas", "1000g2015aug_all",
                        "1000g2015aug_eas", "ExAC_ALL", "ExAC_EAS", "esp6500siv2_all",
                        "CLNDN", "CLNSIG", "ClinvarPathHotRegion", "InterVar_automated", "InterVar_sig", "M-CAP_pred", "REVEL_score", "SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationTaster_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE", "SpliceAI",
                        "MAPQ", "VcfInfo"
                    ]
                    if os.path.exists(geneCoverage):
                        header.insert(12, "GeneDepth")
                        header.insert(12, "GeneCoverage")

                    for tn in trioNames:
                        header.append(tn + "_DP")
                        header.append(tn + "_RefAD")
                        header.append(tn + "_AltAD")
                        header.append(tn + "_VAF")

                    header.append("INFO")
                    
                    for tn in trioNames:
                        header.append(tn + "_INFO")

                    for tn in trioNames:
                        header.append(tn + "_GT")

                    results.write("\t".join(header) + "\n")
                    headCount = 1

                # 处理人群频率，仅考虑以下来源
                keyList = ["AF", "AF_eas", "AF_exome", "AF_exome_eas", "1000g2015aug_all",
                    "1000g2015aug_eas", "ExAC_ALL", "ExAC_EAS"
                ]
                AF_list = []
                for k in keyList:
                    if lineDict[k] == ".":
                        lineDict[k] = "-"
                    AF_list.append(lineDict[k])
                Lower001 = checkAF(AF_list, 0.01)
                Lower005 = checkAF(AF_list, 0.05)
                lineDict["Lower0.01"] = Lower001
                lineDict["Lower0.05"] = Lower005

                # InterVar合并
                lineDict = MergeInterVar(lineDict)

                # 整理输出结果
                output = [lineDict["Chr"], lineDict["Start"], lineDict["End"], lineDict["Ref"],
                    lineDict["Alt"], lineDict["Gene"], lineDict["Type"], lineDict["Transcript"],
                    lineDict["cHGVS"], AATranslate(lineDict["pHGVS"]), lineDict["Consequence"],
                    lineDict["AffectedExon"], lineDict["cytoBand"], lineDict["avsnp150"], lineDict["Lower0.01"], lineDict["Lower0.05"],
                    lineDict["AF"], lineDict["AF_popmax"], lineDict["AF_male"], lineDict["AF_female"],
                    lineDict["AF_raw"], lineDict["AF_eas"], lineDict["AF_afr"], lineDict["AF_sas"],
                    lineDict["AF_amr"], lineDict["AF_nfe"], lineDict["AF_fin"], lineDict["AF_asj"], lineDict["AF_oth"],
                    lineDict["non_topmed_AF_popmax"], lineDict["non_neuro_AF_popmax"], lineDict["non_cancer_AF_popmax"],
                    lineDict["controls_AF_popmax"], lineDict["AF_exome"], lineDict["AF_exome_eas"], lineDict["1000g2015aug_all"], lineDict["1000g2015aug_eas"],
                    lineDict["ExAC_ALL"], lineDict["ExAC_EAS"], lineDict["esp6500siv2_all"],
                    lineDict["CLNDN"], lineDict["CLNSIG"], lineDict["ClinvarPath"], lineDict["InterVar_automated"], lineDict["InterVar_sig"],
                    lineDict["M-CAP_pred"], lineDict["REVEL_score"], lineDict["SIFT_pred"],
                    lineDict["Polyphen2_HDIV_pred"], lineDict["Polyphen2_HVAR_pred"], lineDict["LRT_pred"], lineDict["MutationTaster_pred"],
                    lineDict["MutationAssessor_pred"], lineDict["FATHMM_pred"], lineDict["PROVEAN_pred"], lineDict["dbscSNV_ADA_SCORE"],
                    lineDict["dbscSNV_RF_SCORE"], lineDict["SpliceAI"], lineDict["Otherinfo3"], lineDict["Otherinfo11"]
                ]

                if os.path.exists(geneCoverage):
                    if lineDict["Gene"] in geneCoverageDict.keys():
                        geneDepth = geneCoverageDict[lineDict["Gene"]][0]
                        geneCov = geneCoverageDict[lineDict["Gene"]][1]
                    else:
                        geneCov = geneDepth = "-"
                    output.insert(12, geneDepth)
                    output.insert(12, geneCov)

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

    # 补充omim
    if omim == "True" or omim == "true":
        ## 数据库
        omimDB = "/slurm/databases/humandb/OmimGenePhenotype.txt"
        os.rename(resultsFile, resultsFile + ".tmp")
        appendOMIM(resultsFile + ".tmp", resultsFile, omimDB)
        os.remove(resultsFile + ".tmp")

    # 补充PP2
    if pp2 == "True" or pp2 == "true":
        pp2DB = "/slurm/databases/humandb/PP2GeneList.txt"
        os.rename(resultsFile, resultsFile + ".tmp")
        appendPP2Gene(resultsFile + ".tmp", resultsFile, pp2DB, omim, geneCoverage)
        os.remove(resultsFile + ".tmp")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            AnnoGermline
            在默认的-l False模式下，程序会对同一位置比对到不同基因的结果进行拆分，因此最终输出行数会大于vcf原结果行数。

            适用于snpeff+annovar注释后结果整理。
            其中annovar需要注释数据库为
                    refGene,cytoBand,avsnp150,gnomad211_genome,gnomad211_exome,1000g2015aug_all,1000g2015aug_eas,exac03,
                    esp6500siv2_all,clinvar_20220320,dbnsfp42a,dbscsnv11,intervar_20180118,SpliceAI

            建议命令：
            python3 AnnoGermline.py -i sample.hg19_multianno.txt -o sample.anno.txt -t sample
            python3 AnnoGermline.py -i sample.hg19_multianno.txt -o sample.anno.txt -t sample -l True
            python3 AnnoGermline.py -i sample.hg19_multianno.txt -o sample.anno.txt -t sample -gcov sample.region.txt -fs True
        
        """,
        prog="AnnoMT.py",
        usage="python3 AnnoGermline.py [-h] -i <annoFile> -o <annoFixFile> -t <trioNameList>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version",
        version="Version 0.1 20220505")
    parser.add_argument("-i", "--anno", type=str,
        help="注释结果文件")
    parser.add_argument("-o", "--output", type=str,
        help="注释结果整理")
    parser.add_argument("-t", "--trio", type=str,
        help="受检者名称，按照vcf结果中名称顺序填写，以‘,’分隔")
    parser.add_argument("-r", "--reftran", type=str,
        help="参考转录本，默认为/slurm/databases/hg19/LRG/refTranscript_germline.txt", default="/slurm/databases/hg19/LRG/refTranscript_germline.txt")
    parser.add_argument("-omim", "--omim", type=str,
        help="是否注释omim内容，默认是", default="True")
    parser.add_argument("-syno", "--syno", type=str,
        help="是否过滤特定同义突变，默认是", default="True")
    parser.add_argument("-pp2", "--pp2", type=str,
        help="是否注释PP2基因，默认是", default="True")
    parser.add_argument("-gcov", "--gcov", type=str,
        help="导入基因覆盖结果", default="False")
    parser.add_argument("-l", "--lite", type=str,
        help="同一位置不进行多基因分割，默认否", default="False")
    parser.add_argument("-fs", "--fs", type=str,
        help="过滤miRNA、LOC、LINC系列基因，默认否", default="False")
    parser.add_argument("-s", "--silent", type=str,
        help="是否静默运行，默认是", default="True")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(annovarResultsFile=args.anno, resultsFile=args.output, refTranscriptFile=args.reftran, trioName=args.trio, silent=args.silent, omim=args.omim, syno=args.syno, pp2=args.pp2, geneCoverage=args.gcov, lite=args.lite, filterSymbol=args.fs)

# end
        