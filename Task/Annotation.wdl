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
    }

    String humandb = "/home/novelbio/databases/humandb"

    command <<<
        convert2annovar.pl -format vcf4 -allsample -withfreq \
            ~{vcf} \
            --includeinfo > ~{sample}.avinput
        table_annovar.pl ~{sample}.avinput \
            ~{humandb} -buildver hg19 \
            -out ~{sample} -remove \
            -protocol refGene,cytoBand,avsnp150,gnomad211_genome,gnomad211_exome,1000g2015aug_all,1000g2015aug_eas,exac03,esp6500siv2_all,clinvar_20220812,dbnsfp42a,dbscsnv11,intervar_20180118,SpliceAI \
            -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f \
            -nastring - -otherinfo
    >>>

    output {
        File annovarResults = "~{sample}.hg19_multianno.txt"
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
            -protocol refGene,avsnp150,clinvar_20220812,mitomapB37 \
            -operation g,f,f,f \
            -nastring - -otherinfo
        python3 ~{fixScript} -i ~{sample}.hg19_multianno.txt -o ~{sample}.MT.anno.txt -t ~{sample}
    >>>

    output {
        File mtAnno = "~{sample}.MT.anno.txt"
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

    String vep_cache = "/ykrt/data/backup/databases/vep_data"
    File reference = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta"
    File refFai = "/home/novelbio/databases/b37/human_g1k_v37_decoy.fasta.fai"

    command <<<
        plugins_dir=~{vep_cache}/Plugins
        plu_data_dir=~{vep_cache}/Plu_data
        dbnsfp_str="CADD_phred,SIFT_pred,Polyphen2_HDIV_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,M-CAP_pred,REVEL_score"
        dbnsfp_str=${dbnsfp_str}",clinvar_clnsig"
        fields_str="Uploaded_variation,Location,REF_ALLELE,Allele,Gene,VARIANT_CLASS,CANONICAL,HGVSc,HGVSp,Consequence,EXON,BIOTYPE"
        fields_str=${fields_str}",Existing_variation,gnomAD_EAS_AF,AF,EAS_AF,"
        fields_str=${fields_str}${dbnsfp_str}",ada_score,rf_score,SpliceAI_pred"

        vep \
            -i ~{vcf} \
            -o ~{sample}.vep.vcf \
            --offline --cache \
            --format vcf --refseq --fork ~{threads} \
            --force_overwrite \
            --dir_cache ~{vep_cache} \
            --dir_plugins ${plugins_dir} \
            --plugin dbNSFP,${plu_data_dir}/dbNSFP4.3a_grch37.gz,${dbnsfp_str} \
            --plugin dbscSNV,${plu_data_dir}/dbscSNV1.1_GRCh37.txt.gz \
            --plugin SpliceAI,snv=${plu_data_dir}/spliceai_scores.raw.snv.hg19.vcf.gz,indel=${plu_data_dir}/spliceai_scores.raw.indel.hg19.vcf.gz,cutoff=0.5 \
            --fasta ~{reference} --assembly GRCh37 \
            --shift_3prime 1 --no_escape --show_ref_allele --check_existing \
            --exclude_predicted --canonical --vcf --fields ${fields_str}
    >>>

    output {
        File vepVcf = "~{sample}.vep.vcf"
    }

    runtime {
        cpus: threads
    }

}

# annotSV
# https://github.com/lgmgeo/AnnotSV
task AnnotSV {
    input {
        String sample
        File cnvResult
    }

    String annotSV = "/yk/apps/biosoft/AnnotSV/bin/AnnotSV"
    String annotSVDir = "/yk/apps/biosoft/AnnotSV"
    String script = "/home/novelbio/pipeline/WESpipeWDL/Script/annotsv_fix.py"

    command <<<
        cat ~{cnvResult} | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | sed 's/gain/DUP/' | sed 's/loss/DEL/' > ~{sample}.cnv.bed
        export ANNOTSV=~{annotSVDir}
        ~{annotSV} \
            -SVinputFile ~{sample}.cnv.bed \
            -outputFile ~{sample}.annotsv.tsv \
            -outputDir . \
            -svtBEDcol 4 \
            -annotationMode full \
            -genomeBuild GRCh37 \
            -SVminSize 5000
        python3 ~{script} ~{sample}.annotsv.tsv ~{sample}.annotsv.txt
    >>>

    output {
        File annoResult = "~{sample}.annotsv.txt"
    }

}

# AnnotSv 用于exomedepth结果的版本
task AnnotSVEx {
    input {
        String sample
        File cnvResult
    }

    String annotSV = "/yk/apps/biosoft/AnnotSV/bin/AnnotSV"
    String annotSVDir = "/yk/apps/biosoft/AnnotSV"
    String script = "/home/novelbio/pipeline/WESpipeWDL/Script/annotsv_fix_exomedepth.py"

    command <<<
        cat ~{cnvResult} \
            | awk '{print $7"\t"$5"\t"$6"\t"$3"\t"$9"\t"$10"\t"$11"\t"$12}' \
            | sed 's/duplication/DUP/' \
            | sed 's/deletion/DEL/' \
            | sed 's/chromosome/#chromosome/' > ~{sample}.fix.bed
        export ANNOTSV=~{annotSVDir}
        ~{annotSV} \
            -SVinputFile ~{sample}.fix.bed \
            -outputFile ~{sample}.annotsv.tsv \
            -outputDir . \
            -svtBEDcol 4 \
            -annotationMode full \
            -genomeBuild GRCh37 \
            -SVminSize 5000
        python3 ~{script} ~{sample}.annotsv.tsv ~{sample}.annotsv.txt        
    >>>    

    output {
        File annoResult = "~{sample}.annotsv.txt"
    }
}

