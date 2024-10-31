# WESPipeWDL
## 项目介绍
使用WDL流程语言，从零开始搭建一套普适性的遗传全外显子分析流程。

使用GRCh38！

## 路线

fastp + bwa + samtools + markduplicates -> bam
deepvariant + vep -> snp/indel
mutect2 + vep -> mt
cnvkit -> cnv

