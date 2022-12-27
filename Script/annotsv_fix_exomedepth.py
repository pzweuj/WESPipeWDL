# coding=utf-8
# pzw
# 20221006
# CNV注释后结果格式调整

import sys

# 输入输出
input = open(sys.argv[1], "r", encoding="utf-8")
output = open(sys.argv[2], "w", encoding="utf-8")


# 输出内容标题
outputHeaderList = ["#Chrom", "Start", "End", "Length","Type", "ISCN", "CytoBand", "Gene", "GeneCount", "BayesFactor", "reads.ratio",
    "ExAC_delZ", "ExAC_dupZ", "ExAC_cnvZ", "ExAC_synZ", "ExAC_misZ", "ExAC_pLI", "GnomAD_pLI",
    "Class", "P_gain_phen", "P_gain_hpo", "P_gain_source", "P_gain_coord", "P_loss_phen", "P_loss_hpo", "P_loss_source", "P_loss_coord",
    "B_gain_source", "B_gain_coord", "B_gain_AFmax", "B_loss_source", "B_loss_coord", "B_loss_AFmax",
    "Decipher_HI", "OMIM_ID", "OMIM_phenotype", "OMIM_inheritance", "OMIM_morbid"
]
outputheader = "\t".join(outputHeaderList)
output.write(outputheader + "\n")

# 整理结果
lineDict = {}
numDict = {}
for line in input:
    if line.startswith("AnnotSV_"):
        header = line.replace("\n", "").split("\t")
        for h in range(len(header)):
            lineDict[header[h]] = ""
            numDict[header[h]] = h

    else:
        lines = line.replace("\n", "").split("\t")
        for i in range(len(lines)):
            for k in numDict:
                if numDict[k] == i:
                    lineDict[k] = lines[i]
        lineDict["ISCN"] = "seq[GRCh37]" + lineDict["SV_type"].lower() + "(" + lineDict["SV_chrom"] + ")(" + lineDict["CytoBand"].replace("-", "") + ")chr"
        lineDict["ISCN"] = lineDict["ISCN"] + lineDict["SV_chrom"] + ":g." + lineDict["SV_start"] + "_" + lineDict["SV_end"] + lineDict["SV_type"].lower()
        
        # Class
        if lineDict["ACMG_class"] == "1":
            lineDict["Class"] = "Benign"
        elif lineDict["ACMG_class"] == "2":
            lineDict["Class"] = "Likely_benign"
        elif lineDict["ACMG_class"] == "4":
            lineDict["Class"] = "Likely_pathogenic"
        elif lineDict["ACMG_class"] == "5":
            lineDict["Class"] = "Pathogenic"
        else:
            lineDict["Class"] = "VUS"

        outputStringList = [lineDict["SV_chrom"], lineDict["SV_start"], lineDict["SV_end"], lineDict["SV_length"], lineDict["SV_type"], lineDict["ISCN"], lineDict["CytoBand"], lineDict["Gene_name"], lineDict["Gene_count"], lineDict["BF"], lineDict["reads.ratio"],
            lineDict["ExAC_delZ"], lineDict["ExAC_dupZ"], lineDict["ExAC_cnvZ"], lineDict["ExAC_synZ"], lineDict["ExAC_misZ"], lineDict["ExAC_pLI"], lineDict["GnomAD_pLI"],
            lineDict["Class"], lineDict["P_gain_phen"], lineDict["P_gain_hpo"], lineDict["P_gain_source"], lineDict["P_gain_coord"], lineDict["P_loss_phen"], lineDict["P_loss_hpo"], lineDict["P_loss_source"], lineDict["P_loss_coord"],
            lineDict["B_gain_source"], lineDict["B_gain_coord"], lineDict["B_gain_AFmax"], lineDict["B_loss_source"], lineDict["B_loss_coord"], lineDict["B_loss_AFmax"],
            lineDict["DDD_HI_percent"], lineDict["OMIM_ID"], lineDict["OMIM_phenotype"], lineDict["OMIM_inheritance"], lineDict["OMIM_morbid"]
        ]
        outputString = "\t".join(outputStringList)
        output.write(outputString + "\n")

output.close()
input.close()




