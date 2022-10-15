version 1.0
# panzhaowen
# 20220914
# 胚系CNV分析
# cnv的输入输出理论上是不固定数量的bam文件
# 分为使用基线及使用对照样本的两种模式


## WES CNV分析





## WGS CNV分析
## CNV-seq分析
### 输入是一个bam文件以及一个作为reference的bam的路径list，后续可以建立基线取代
task Cnvkit {
    input {
        String sample
        File testBam
        File controlBamPathFile
        Int threads
    }

    File reference = "/slurm/databases/b37/human_g1k_v37_decoy.fasta"
    File refFai = "/slurm/databases/b37/human_g1k_v37_decoy.fasta.fai"
    File refFlat = "/slurm/databases/b37/b37_refFlat.txt"

    command <<<
        controlBamList=`cat ~{controlBamPathFile} | tr '\n' ' '`
        cnvkit.py batch ~{testBam} -n $controlBamList \
            -m wgs -f ~{reference} --annotate ~{refFlat} \
            --short-names --target-avg-size 20000 -p ~{threads} \
            --segment-method hmm
        out=`basename ~{testBam} .bam`
        mv ${out}.call.cns ~{sample}.call.cns
    >>>

    output {
        File cnvResults = "~{sample}.call.cns"
    }

    runtime {
        cpus: threads
        docker: "etal/cnvkit:latest"
    }

}

## 单样本分析，使用Control-Freec
task Freec {
    input {
        String sample
        File bam
        File bai
        String sex
        Int threads
    }

    File configFile = "/home/novelbio/pipeline/WESpipeWDL/Config/config_cnv-seq.txt"

    command <<<
        sed "s|{outputDir}|${PWD}|" ~{configFile} > ~{sample}.config
        sed -i "s/{threads}/~{threads}/" ~{sample}.config
        sed -i "s/{sex}/~{sex}/" ~{sample}.config
        sed -i "s|{bamFile}|~{bam}|" ~{sample}.config
        /home/novelbio/software/FREEC-11.6/src/freec -conf ~{sample}.config
        name=`basename ~{bam} .bam`
        mv ${name}.bam_CNVs ~{sample}.CNVs.txt
        mv ${name}.bam_ratio.txt ~{sample}.ratio.txt
        perl /home/novelbio/software/FREEC-11.6/scripts/freec2bed.pl -f ~{sample}.ratio.txt | tr ' ' '\t' > ~{sample}.ratio.bed
        perl /home/novelbio/software/FREEC-11.6/scripts/freec2circos.pl -f ~{sample}.ratio.txt > ~{sample}.ratio.circos
    >>>

    output {
        File cnvResult = "~{sample}.CNVs.txt"
        File ratioResult = "~{sample}.ratio.txt"
    }

    runtime {
        cpu: threads
    }

}

# 画图
task FreecPlot {
    input {
        String sample
        File ratioFile
        File cnvResult
    }

    command <<<
        Rscript /software/FREEC-11.6/scripts/FreecChromPlot.R ~{ratioFile} ~{sample}_chrom_pics
        cat /software/FREEC-11.6/scripts/assess_significance.R | R --slave --args ~{cnvResult} ~{ratioFile}
        cat /software/FREEC-11.6/scripts/makeGraph.R | R --slave --args 2 ~{ratioFile}
        path=`dirname ~{ratioFile}`
        name=`basename ~{ratioFile} .txt`
        mv ${path}/${name}.txt.log2.png ~{sample}.freec.log2.png
        mv ${path}/${name}.txt.png ~{sample}.freec.png
        name=`basename ~{cnvResult} .txt`
        mv ${path}/${name}.txt.p.value.txt ~{sample}.cnv.p.value.txt
    >>>

    output {
        Array[File] outputPics = glob("~{sample}_chrom_pics/~{sample}.*.png")
        File pValueFile = "~{sample}.cnv.p.value.txt"
        File log2Pic = "~{sample}.freec.log2.png"
        File cnPic = "~{sample}.freec.png"
    }

    runtime {
        docker: "pzweuj/freec:v11.6"
    }

}


