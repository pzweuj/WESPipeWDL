version 1.0
# panzhaowen
# 20220427
# WES 注释

# snpeff
task Snpeff {
    input {
        String sample
        File vcf
    }

    String SNPEFF = "/home/novelbio/software/snpEff-5.0d/snpEff.jar"
    File config = "/home/novelbio/software/snpEff-5.0d/snpEff.config"

    command <<<
        java -jar ~{SNPEFF} \
            -c ~{config} \
            hg19 ~{vcf} > ~{sample}.snpeff.vcf
    >>>

    output {
        File annoVcf = "~{sample}.snpeff.vcf"
    }
}

# annovar
task Annovar {
    input {
        String sample
        File vcf
        Int threads
    }

    String humandb = "/home/novelbio/databases/humandb"

    command <<<
        convert2annovar.pl -format vcf4 -allsample -withfreq \
            ~{vcf} \
            --includeinfo > ~{sample}.avinput
        table_annovar.pl ~{sample}.avinput \
            ~{humandb} -buildver hg19 \
            -out ~{sample} -remove \
            -protocol refGene,cytoBand,avsnp150,gnomad211_genome,gnomad211_exome,1000g2015aug_all,1000g2015aug_eas,exac03,esp6500siv2_all,clinvar_20220320,dbnsfp42a,dbscsnv11,intervar_20180118,SpliceAI \
            -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f \
            -nastring - -thread ~{threads} -otherinfo
    >>>

    output {
        File annovarResults = "~{sample}.hg19_multianno.txt"
    }

    runtime {
        cpus: threads
    }

}

# 基因覆盖度分析
task GeneCoverage {
    input {
        String sample
        File regionGz
        File bed
    }

    String script = "/home/novelbio/pipeline/WESpipeWDL/Script/GeneCover.py"

    command <<<
        python3 ~{script} ~{regionGz} ~{bed} ~{sample}.geneCover.txt
    >>>

    output {
        File geneCover = "~{sample}.geneCover.txt"
    }

}


# 线粒体注释
task MTAnnotation {
    input {
        String sample
        File vcf
        Int threads
    }

    String SNPEFF = "/home/novelbio/software/snpEff-5.0d/snpEff.jar"
    File config = "/home/novelbio/software/snpEff-5.0d/snpEff.config"
    String fixScript = "/home/novelbio/pipeline/WESpipeWDL/Script/AnnoMT.py"

    command <<<
        java -jar ~{SNPEFF} -c ~{config} MT ~{vcf} > ~{sample}.MT.snpeff.vcf
        convert2annovar.pl -format vcf4 -allsample -withfreq \
            ~{sample}.MT.snpeff.vcf --includeinfo > ~{sample}.avinput
        table_annovar.pl ~{sample}.avinput \
            /home/novelbio/databases/humandb \
            -buildver hg19 -out ~{sample} -remove \
            -protocol refGene,avsnp150,clinvar_20220227,mitomapB37 \
            -operation g,f,f,f \
            -nastring - -thread ~{threads} -otherinfo
        python3 ~{fixScript} -i ~{sample}.hg19_multianno.txt -o ~{sample}.MT.anno.txt -t ~{sample}
    >>>

    output {
        File mtAnno = "~{sample}.MT.anno.txt"
    }

    runtime {
        cpus: threads
    }

}

# 注释结果整理
task AnnotationFix {
    input {
        String sample
        File annoFile
        File? geneCoverFile
    }

    String fixScript = "/home/novelbio/pipeline/WESpipeWDL/Script/AnnoGermline.py"
    String clinvarPath = "/slurm/databases/humandb/b37_clinvarPathRegion.txt"

    command <<<
        head -n 1 ~{annoFile} | sed 's/$/&\tClinvarPath/' > ~{sample}.anno.header
        tail -n +2 ~{annoFile} | bedtools intersect -a - -b ~{clinvarPath} -c > ~{sample}.anno.sig
        cat ~{sample}.anno.header ~{sample}.anno.sig > ~{sample}.anno.tmp
        python3 ~{fixScript} -i ~{sample}.anno.tmp -o ~{sample}.anno.txt -t ~{sample} -gcov ~{geneCoverFile} -syno False
    >>>

    output {
        File resultsFile = "~{sample}.anno.txt"
    }

}

# VEP
task VEPAnno {
    input {
        String sample
        File vcf
        Int threads
    }

    String vep_cache = "/ykrt/data/backup/databases/VEP_cache"
    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File refFai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

    command <<<
        vep \
            -i ~{vcf} \
            -o ~{sample}.vep.vcf \
            --offline \
            --assembly GRCh37 \
            --cache \
            --dir_cache ~{vep_cache} \
            --fasta ~{reference} \
            --vcf \
            --fork ~{threads}
    >>>

    output {
        File vepVcf = "~{sample}.vep.vcf"
    }

    runtime {
        cpus: threads
    }

}


