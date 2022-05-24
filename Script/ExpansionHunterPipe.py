#!/usr/local/bin/python3
# coding=utf-8
# pzw
# 20220519
# ExpansionHunter STR 分析流程
# 包含REViewer


import os
import sys
import cairosvg
import argparse

# ExpansionHunter 分析流程
def ExpansionHunter(bam, reference, catalog, prefix):
    EH = "/home/novelbio/software/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter"
    cmd = """
        {EH} \\
            --reads {bam} \\
            --reference {reference} \\
            --variant-catalog {catalog} \\
            --output-prefix {prefix}
    """.format(EH=EH, bam=bam, reference=reference, catalog=catalog, prefix=prefix)
    os.system(cmd)

# REViewer 作图流程
def REViewer(bam, vcf, reference, catalog, locus, prefix):
    RV = "/home/novelbio/software/ExpansionHunter-v5.0.0-linux_x86_64/bin/REViewer"
    sortbam = bam.replace(".bam", ".sort.bam")
    cmd = """
        samtools sort {bam} -o {sortbam}
        samtools index {sortbam}
        {RV} \\
            --reads {sortbam} \\
            --vcf {vcf} \\
            --reference {reference} \\
            --catalog {catalog} \\
            --locus {locus} \\
            --output-prefix {prefix}
    """.format(bam=bam, sortbam=sortbam, RV=RV, vcf=vcf, reference=reference, catalog=catalog, locus=locus, prefix=prefix)
    os.system(cmd)

# 判定是否致病 hg19 Only
def DiseaseJudger(chrom, pos, STR1, STR2):
    db = [
        ["chr3", 63898357, "ATXN7", "Spinocerebellar Ataxia Type 7", "36-"],
        ["chr4", 3076599, "HTT", "Huntington Disease", "40-"],
        ["chr5", 146258287, "PPP2R2B", "Spinocerebellar Ataxia Type 12", "51-78"],
        ["chr6", 16327862, "ATXN1", "Spinocerebellar Ataxia Type 1", "39-"],
        ["chr6", 170871009, "TBP", "Spinocerebellar Ataxia Type 17", "49-"],
        ["chr9", 71652198, "FXN", "Friedreich's Ataxia", "66-"],
        ["chr12", 7045879, "ATN1", "Dentatorubral-pallidoluysian Atrophy", "48-"],
        ["chr12", 112036750, "ATXN2", "Spinocerebellar Ataxia Type 2", "33-"],
        ["chr13", 70713481, "ATXN8OS", "Spinocerebellar Ataxia Type 8", "51-"],
        ["chr14", 92537350, "ATXN3", "Spinocerebellar Ataxia Type 3", "60-87"],
        ["chr18", 53253386, "TCF4", "Fuchs Corneal Dystrophy", "50-"],
        ["chr19", 13318668, "CACNA1A", "Spinocerebellar Ataxia Type 6", "20-33"],
        ["chr19", 46273458, "DMPK", "Myotonic Dystrophy Type 1", "50-"],
        ["chrX", 146993568, "FMR1", "Fragile X Syndrome", "55-"]
    ]

    if not "chr" in chrom:
        chrom = "chr" + chrom
    
    disease = "-"
    status = "-"
    for d in db:
        if chrom == d[0]:
            if abs(int(pos) - d[1]) <= 10:
                disease = d[3]
                region = d[4].split("-")
                if region[0] == "":
                    if region[1] == "":
                        pass
                    else:
                        if (int(STR1) <= int(region[1])) or (int(STR2) <= int(region[1])):
                            status = "+"
                else:
                    if region[1] == "":
                        if (int(STR1) >= int(region[0])) or (int(STR2) >= int(region[0])):
                            status = "+"
                    else:
                        if (int(STR1) >= int(region[0])) and (int(STR1) <= int(region[1])):
                            status = "+"
                        
                        if (int(STR2) >= int(region[0])) and (int(STR2) <= int(region[1])):
                            status = "+"
    
    return [disease, status]


def VcfReader(vcf, output):
    vcfO = open(vcf, "r", encoding="utf-8")
    out = open(output, "w", encoding="utf-8")
    
    #  header
    out.write("#Chrom\tStart\tEnd\tGene\trefUnit\trefRepeat\tSTR1\tSTR2\tDisease\tJudgement\n")

    geneList = []
    for line in vcfO:
        if not line.startswith("#"):
            lines = line.split("\t")
            chrom = lines[0]
            pos = lines[1]
            filter = lines[6]
            alt = lines[4]
            
            if filter != "PASS":
                continue
        
            if not "STR" in alt:
                continue

            alts = alt.split(",")

            infos = lines[7].split(";")
            end = "-"
            STR1 = STR2 = refRepeat = "-"
            # refLength = "-"
            refUnit = "-"
            gene = "-"
            for i in infos:
                if i.startswith("END="):
                    end = i.split("=")[1]
                elif i.startswith("REF="):
                    refRepeat = i.split("=")[1]
                # elif i.startswith("RL="):
                #     refLength = i.split("=")[1]
                elif i.startswith("RU="):
                    refUnit = i.split("RU=")[1]
                elif i.startswith("REPID="):
                    gene = i.split("REPID=")[1].split("_")[0]
                    geneList.append(gene)
                else:
                    pass
            
            if len(alts) == 1:
                STR1 = refRepeat
                STR2 = alts[0].replace("<STR", "").replace(">", "")
            elif len(alts) == 2:
                STR1 = alts[0].replace("<STR", "").replace(">", "")
                STR2 = alts[1].replace("<STR", "").replace(">", "")
            else:
                pass
            
            diseList = DiseaseJudger(chrom, pos, STR1, STR2)
            
            outLine = [chrom, pos, end, gene, refUnit, refRepeat, STR1, STR2, diseList[0], diseList[1]]
            out.write("\t".join(outLine) + "\n")

    out.close()
    vcfO.close()

    geneList = list(set(geneList))
    return geneList

# Main
def main(bam, reference, outputPrefix, gnomeVersion, outputPics):
    sample = outputPrefix.split("/")[-1]
    outputDir = "/".join(outputPrefix.split("/")[0:-1])
    if outputDir == "":
        outputDir = "./"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    catalogDir = "/home/novelbio/software/ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog"
    if gnomeVersion == "grch38":
        catalog = catalogDir + "/grch38/variant_catalog.json"
    elif gnomeVersion == "hg19":
        catalog = catalogDir + "/hg19/variant_catalog.json"
    elif gnomeVersion == "hg38":
        catalog = catalogDir + "/hg38/variant_catalog.json"
    else:
        catalog = catalogDir + "/grch37/variant_catalog.json"

    ExpansionHunter(bam, reference, catalog, outputPrefix)
    geneList = VcfReader(outputPrefix + ".vcf", outputPrefix + ".EH.txt")
    if not os.path.exists(outputDir + "/pics"):
        os.makedirs(outputDir + "/pics")
    
    if outputPics == "True" or outputPics == "true":
        if len(geneList) >= 1:
            geneString = ",".join(geneList)
            REViewer(outputPrefix + "_realigned.bam", outputPrefix + ".vcf", reference, catalog, geneString, outputDir + "/pics/" + sample)
            for pic in os.listdir(outputDir + "/pics"):
                if pic.endswith(".svg"):    
                    cairosvg.svg2pdf(url=outputDir + "/pics/" + pic, write_to=outputDir + "/pics/" + pic.replace(".svg", ".pdf"))
        else:
            print("Cannot Create Pics without results")
        
# 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            ExpansionHunter Pipeline
        """,
        prog="ExpansionHunterPipe.py",
        usage="python3 ExpansionHunterPipe.py [-h] -i <bamFile> -r <referenceFasta> -o <outputPrefix>",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version",
        version="Version 0.1 20220519")
    parser.add_argument("-i", "--input", type=str,
        help="bam文件")
    parser.add_argument("-o", "--output", type=str,
        help="结果文件")
    parser.add_argument("-r", "--ref", type=str,
        help="参考基因组")
    parser.add_argument("-g", "--genome", type=str,
        help="可选grch37,grch38,hg19,hg38，默认为grch37", default="grch37")
    parser.add_argument("-p", "--pics", type=str,
        help="是否作图，默认是", default="True")
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    main(bam=args.input, reference=args.ref, outputPrefix=args.output, gnomeVersion=args.genome, outputPics=args.pics)

# end

