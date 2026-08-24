# ---------------------------------------------------------------------------
# Per-sample processing: trim -> STAR -> unique filter -> dedup -> index
# PE and SE handled by parallel rules constrained to the relevant sample set.
# ---------------------------------------------------------------------------

# ---- Trimming ------------------------------------------------------------
rule trim_pe:
    input:
        r1 = lambda wc: meta.loc[wc.sample, "fastq_r1"],
        r2 = lambda wc: meta.loc[wc.sample, "fastq_r2"],
    output:
        r1 = p("TRIMMED", "{sample}_R1.fq.gz"),
        r2 = p("TRIMMED", "{sample}_R2.fq.gz"),
    params:
        tdir    = p("TRIMMED"),
        adapter = "--" + config["trim"]["adapter"],
        q       = config["trim"]["quality"],
        extra   = config["trim"]["extra"],
    threads: min(config["threads"], 8)
    resources:
        mem_mb = 16000, runtime = 180
    wildcard_constraints:
        sample = constr(PE_SAMPLES)
    log:
        p("logs", "trim_{sample}.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        trim_galore {params.adapter} --paired {input.r1} {input.r2} \
            -q {params.q} -j {threads} {params.extra} -o {params.tdir} > {log} 2>&1
        b1=$(basename {input.r1} .fastq.gz); b2=$(basename {input.r2} .fastq.gz)
        mv {params.tdir}/${{b1}}_val_1.fq.gz {output.r1}
        mv {params.tdir}/${{b2}}_val_2.fq.gz {output.r2}
        """

rule trim_se:
    input:
        r1 = lambda wc: meta.loc[wc.sample, "fastq_r1"],
    output:
        r1 = p("TRIMMED", "{sample}_R1.fq.gz"),
    params:
        tdir    = p("TRIMMED"),
        adapter = "--" + config["trim"]["adapter"],
        q       = config["trim"]["quality"],
        extra   = config["trim"]["extra"],
    threads: min(config["threads"], 8)
    resources:
        mem_mb = 16000, runtime = 180
    wildcard_constraints:
        sample = constr(SE_SAMPLES)
    log:
        p("logs", "trim_{sample}.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        trim_galore {params.adapter} {input.r1} \
            -q {params.q} -j {threads} {params.extra} -o {params.tdir} > {log} 2>&1
        b1=$(basename {input.r1} .fastq.gz)
        mv {params.tdir}/${{b1}}_trimmed.fq.gz {output.r1}
        """

# ---- Alignment (STAR) ----------------------------------------------------
_STAR_COMMON = (
    "STAR --runMode alignReads --runThreadN {threads} --readFilesCommand zcat "
    "--genomeDir {params.index} --alignEndsType {params.ends} "
    "--genomeLoad NoSharedMemory --outFilterMultimapNmax {params.mm} "
    "--outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate "
    "--outFilterScoreMin {params.score} --outFilterMatchNmin {params.match} {params.extra} "
)

rule star_pe:
    input:
        r1 = p("TRIMMED", "{sample}_R1.fq.gz"),
        r2 = p("TRIMMED", "{sample}_R2.fq.gz"),
    output:
        bam = p("BAMS", "{sample}Aligned.sortedByCoord.out.bam"),
    params:
        index  = config["star_index"],
        prefix = lambda wc: p("BAMS", wc.sample),
        ends   = config["star"]["align_ends_type"],
        gap    = config["star"]["align_mates_gap_max"],
        mm     = config["star"]["out_filter_multimap_nmax"],
        score  = config["star"]["out_filter_score_min"],
        match  = config["star"]["out_filter_match_nmin"],
        extra  = config["star"]["extra"],
    threads: config["threads"]
    resources:
        mem_mb = 40000, runtime = 240
    wildcard_constraints:
        sample = constr(PE_SAMPLES)
    log:
        p("logs", "star_{sample}.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        _STAR_COMMON +
        "--alignMatesGapMax {params.gap} --readFilesIn {input.r1} {input.r2} > {log} 2>&1"

rule star_se:
    input:
        r1 = p("TRIMMED", "{sample}_R1.fq.gz"),
    output:
        bam = p("BAMS", "{sample}Aligned.sortedByCoord.out.bam"),
    params:
        index  = config["star_index"],
        prefix = lambda wc: p("BAMS", wc.sample),
        ends   = config["star"]["align_ends_type"],
        mm     = config["star"]["out_filter_multimap_nmax"],
        score  = config["star"]["out_filter_score_min"],
        match  = config["star"]["out_filter_match_nmin"],
        extra  = config["star"]["extra"],
    threads: config["threads"]
    resources:
        mem_mb = 40000, runtime = 240
    wildcard_constraints:
        sample = constr(SE_SAMPLES)
    log:
        p("logs", "star_{sample}.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        _STAR_COMMON + "--readFilesIn {input.r1} > {log} 2>&1"

# ---- Unique-mapper filter ------------------------------------------------
rule unique_filter:
    input:
        p("BAMS", "{sample}Aligned.sortedByCoord.out.bam"),
    output:
        temp(p("BAMS", "{sample}_uniq.bam")),
    params:
        q = config["mapq_min"],
    threads: 8
    resources:
        mem_mb = 16000, runtime = 120
    conda:
        "../envs/uvrip.yaml"
    shell:
        "samtools view -@ {threads} -q {params.q} -o {output} {input}"

# ---- Deduplication -------------------------------------------------------
rule dedup:
    input:
        p("BAMS", "{sample}_uniq.bam"),
    output:
        bam = p("BAMS", "nodup_{sample}.bam"),
        bai = p("BAMS", "nodup_{sample}.bam.bai"),
    params:
        paired = lambda wc: is_paired(wc.sample),
        rmdup  = config["remove_duplicates"],
        tmp    = lambda wc: p("BAMS", f"_{wc.sample}"),
    threads: 8
    resources:
        mem_mb = 16000, runtime = 180
    log:
        p("logs", "dedup_{sample}.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        set -euo pipefail
        if [ "{params.rmdup}" = "True" ]; then
            if [ "{params.paired}" = "True" ]; then
                samtools sort   -@ {threads} -n -o {params.tmp}_n.bam {input}
                samtools fixmate -@ {threads} -m {params.tmp}_n.bam {params.tmp}_fm.bam
                samtools sort   -@ {threads} -o {params.tmp}_c.bam {params.tmp}_fm.bam
                samtools markdup -@ {threads} -r -s {params.tmp}_c.bam {output.bam} > {log} 2>&1
                rm -f {params.tmp}_n.bam {params.tmp}_fm.bam {params.tmp}_c.bam
            else
                samtools sort   -@ {threads} -o {params.tmp}_c.bam {input}
                samtools markdup -@ {threads} -r -s {params.tmp}_c.bam {output.bam} > {log} 2>&1
                rm -f {params.tmp}_c.bam
            fi
        else
            samtools sort -@ {threads} -o {output.bam} {input}
        fi
        samtools index -@ {threads} {output.bam}
        """
