version 1.0
# WES 分析流程
# 质控模块

# fastp过滤数据
task Fastp {
    input {
        String sample_id
        String output_dir
        File read1
        File read2
        Int threads
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        fastp \
            -i ~{read1} \
            -I ~{read2} \
            -o ~{output_dir}/~{sample_id}_R1.fq.gz \
            -O ~{output_dir}/~{sample_id}_R2.fq.gz \
            -w ~{threads} \
            -j ~{output_dir}/~{sample_id}.fastp_stats.json \
            -h ~{output_dir}/~{sample_id}.fastp_stats.html \
            --detect_adapter_for_pe
    >>>

    output {
        File cleanRead1 = "~{output_dir}/~{sample_id}_R1.fq.gz"
        File cleanRead2 = "~{output_dir}/~{sample_id}_R2.fq.gz"
        File jsonReport = "~{output_dir}/~{sample_id}.fastp_stats.json"
        File htmlReport = "~{output_dir}/~{sample_id}.fastp_stats.html"
    }

    runtime {
        cpus: threads
    }
}


# bam QC
task Bamdst {
    input {
        String sample_id
        String output_dir
        File bam
        File bai
        File target_bed
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        bamdst -p ~{target_bed} \
            --cutoffdepth 20 \
            -o ~{output_dir} ~{bam}
    >>>

    output {
        String output_path = "~{output_dir}"
        File targetCovFile = "~{output_dir}/coverage.report"
    }

}


# 插入片段分析
task CollectQCMetrics {
    input {
        String sample_id
        File bam
        File bai
        File bed
        String output_dir
        File reference
        File reference_dict
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        gatk CollectInsertSizeMetrics \
            -I ~{bam} \
            -O ~{output_dir}/~{sample_id}.insertsize.txt \
            -H ~{output_dir}/~{sample_id}.histogram.pdf

        gatk CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -R ~{reference} \
            -O ~{output_dir}/~{sample_id}.metrics.txt

        gatk BedToIntervalList \
            -I ~{bed} \
            -O ~{sample_id}.interval_list \
            -SD ~{reference}

        gatk CollectHsMetrics \
            -BI ~{sample_id}.interval_list \
            -TI ~{sample_id}.interval_list \
            -I ~{bam} \
            -O ~{output_dir}/~{sample_id}.hs.txt

    >>>

    output {
        File insertsizeFile = "~{output_dir}/~{sample_id}.insertsize.txt"
        File hisogramPDF = "~{output_dir}/~{sample_id}.histogram.pdf"
        File summary = "~{output_dir}/~{sample_id}.metrics.txt"
        File hs_metric = "~{output_dir}/~{sample_id}.hs.txt"
    }

}

# 性别检测
task GenderPredict {
    input {
        String sample_id
        File bam
        File bai
        String output_dir
    }

    command <<<
        if [ ! -d ~{output_dir} ]; then
            mkdir -p ~{output_dir}
        fi

        samtools view ~{bam} chrY:2786855-2787682 | wc -l > ~{output_dir}/~{sample_id}.SRY.count.txt
    >>>

    output {
        File SRYcount = "~{output_dir}/~{sample_id}.SRY.count.txt"
    }
}

# end
