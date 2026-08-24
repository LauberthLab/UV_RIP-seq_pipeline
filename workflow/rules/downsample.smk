# ---------------------------------------------------------------------------
# Depth-matching + tracks: plan -> subsample -> merge -> bamCompare/bamCoverage
# ---------------------------------------------------------------------------

# ---- 1. Downsample plan (needs every nodup BAM) --------------------------
rule downsample_plan:
    input:
        expand(p("BAMS", "nodup_{s}.bam"), s=SAMPLES),
    output:
        p("BAMS", "downsample_plan.tsv"),
    params:
        bamdir  = p("BAMS"),
        meta    = config["metadata"],
        mode    = config["downsample"]["mode"],
        target  = config["downsample"]["downsample_reads"],
        seed    = config["downsample"]["seed"],
        groupby = config["downsample"]["group_by"],
        script  = "workflow/scripts/compute_downsample.py",
    resources:
        mem_mb = 8000, runtime = 60
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        python {params.script} \
            --bam-dir {params.bamdir} --meta {params.meta} \
            --mode {params.mode} --target-reads {params.target} \
            --seed {params.seed} --group-by {params.groupby} \
            --out {output}
        """

# ---- 2. Subsample each library to its planned fraction -------------------
rule subsample:
    input:
        bam  = p("BAMS", "nodup_{sample}.bam"),
        plan = p("BAMS", "downsample_plan.tsv"),
    output:
        bam = p("BAMS", "nodup_{sample}_sub.bam"),
        bai = p("BAMS", "nodup_{sample}_sub.bam.bai"),
    threads: 8
    resources:
        mem_mb = 8000, runtime = 90
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        set -euo pipefail
        row=$(awk -F'\t' -v s="{wildcards.sample}" '$1==s {{print; exit}}' {input.plan})
        sarg=$(echo "$row" | cut -f6)
        need=$(echo "$row" | cut -f7)
        if [ "$need" = "true" ]; then
            samtools view -b -@ {threads} -s "$sarg" -o {output.bam} {input.bam}
        else
            cp {input.bam} {output.bam}
        fi
        samtools index -@ {threads} {output.bam}
        """

# ---- 3. Merge replicates per (condition, type) ---------------------------
def merge_inputs(wc):
    return [p("BAMS", f"nodup_{s}_sub.bam") for s in members(wc.condition, wc.type)]

rule merge_group:
    input:
        merge_inputs,
    output:
        bam = p("MERGED", "UVRIP_{condition}_{type}_merged.bam"),
        bai = p("MERGED", "UVRIP_{condition}_{type}_merged.bam.bai"),
    threads: 8
    resources:
        mem_mb = 16000, runtime = 120
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        set -euo pipefail
        n=$(echo {input} | wc -w)
        if [ "$n" -eq 1 ]; then
            cp {input} {output.bam}
        else
            samtools merge -f -@ {threads} {output.bam} {input}
        fi
        samtools index -@ {threads} {output.bam}
        """

# ---- 4. bamCompare: log2(IP / control) per condition ---------------------
rule bigwig_log2:
    input:
        ip   = p("MERGED", "UVRIP_{condition}_IP_merged.bam"),
        ipi  = p("MERGED", "UVRIP_{condition}_IP_merged.bam.bai"),
        ctrl = p("MERGED", "UVRIP_{condition}_" + CONTROL + "_merged.bam"),
        ctli = p("MERGED", "UVRIP_{condition}_" + CONTROL + "_merged.bam.bai"),
    output:
        p("BIGWIG", "UVRIP_{condition}_IP_over_" + CONTROL + "_log2ratio.bw"),
    params:
        norm  = config["bigwig"]["normalize_using"],
        scale = config["bigwig"]["scale_factors_method"],
        bs    = config["bigwig"]["bin_size"],
        extra = config["bigwig"]["extra"],
    threads: config["threads"]
    resources:
        mem_mb = 16000, runtime = 120
    log:
        p("logs", "bamcompare_{condition}.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        bamCompare -b1 {input.ip} -b2 {input.ctrl} \
            --normalizeUsing {params.norm} --scaleFactorsMethod {params.scale} \
            -bs {params.bs} {params.extra} -o {output} > {log} 2>&1
        """

# ---- 4b. Optional per-group coverage tracks ------------------------------
rule bigwig_coverage:
    input:
        bam = p("MERGED", "UVRIP_{condition}_{type}_merged.bam"),
        bai = p("MERGED", "UVRIP_{condition}_{type}_merged.bam.bai"),
    output:
        p("BIGWIG", "UVRIP_{condition}_{type}_" + config["bigwig"]["normalize_using"] + ".bw"),
    params:
        norm  = config["bigwig"]["normalize_using"],
        bs    = config["bigwig"]["bin_size"],
        extra = config["bigwig"]["extra"],
    threads: config["threads"]
    resources:
        mem_mb = 16000, runtime = 120
    log:
        p("logs", "bamcoverage_{condition}_{type}.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        bamCoverage -b {input.bam} --normalizeUsing {params.norm} \
            -bs {params.bs} {params.extra} -o {output} > {log} 2>&1
        """
