# ---------------------------------------------------------------------------
# OPTIONAL differential enrichment (enabled by config["de"]["enabled"]).
# featureCounts over nodup BAMs -> tidy -> DESeq2 (de_analysis.R).
# ---------------------------------------------------------------------------
rule featurecounts:
    input:
        bams = expand(p("BAMS", "nodup_{s}.bam"), s=SAMPLES),
    output:
        raw = p("counts", "raw_counts.txt"),
    params:
        gtf     = config["gtf"],
        strand  = config["de"]["strandedness"],
        paired  = is_paired(SAMPLES[0]) if SAMPLES else False,
    threads: config["threads"]
    resources:
        mem_mb = 16000, runtime = 180
    log:
        p("logs", "featurecounts.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        set -euo pipefail
        pflag=""
        if [ "{params.paired}" = "True" ]; then pflag="-p --countReadPairs"; fi
        featureCounts -T {threads} -s {params.strand} $pflag \
            -a {params.gtf} -o {output.raw} {input.bams} > {log} 2>&1
        """

rule tidy_counts:
    input:
        raw = p("counts", "raw_counts.txt"),
    output:
        tsv = p("counts", "counts.tsv"),
    resources:
        mem_mb = 8000, runtime = 30
    conda:
        "../envs/uvrip.yaml"
    run:
        import os
        import pandas as pd
        df = pd.read_csv(input.raw, sep="\t", comment="#")
        meta_cols = ["Chr", "Start", "End", "Strand", "Length"]
        count_cols = [c for c in df.columns if c not in ["Geneid"] + meta_cols]
        rename = {c: os.path.basename(c).replace("nodup_", "").replace(".bam", "")
                  for c in count_cols}
        df = df.rename(columns=rename)
        df[["Geneid"] + [rename[c] for c in count_cols]].to_csv(
            output.tsv, sep="\t", index=False)

rule deseq2:
    input:
        counts = p("counts", "counts.tsv"),
    output:
        touch(p("de", "deseq2.done")),
    params:
        meta    = config["metadata"],
        outdir  = p("de"),
        control = config["de"]["control_type"],
        min_cpm = config["de"]["min_cpm"],
        lfc     = config["de"]["lfc_threshold"],
        padj    = config["de"]["padj_threshold"],
        script  = "workflow/scripts/de_analysis.R",
    resources:
        mem_mb = 16000, runtime = 60
    log:
        p("logs", "deseq2.log")
    conda:
        "../envs/deseq2.yaml"
    shell:
        r"""
        Rscript {params.script} \
            --counts {input.counts} --metadata {params.meta} \
            --outdir {params.outdir} --control {params.control} \
            --min-cpm {params.min_cpm} --lfc {params.lfc} --padj {params.padj} \
            > {log} 2>&1
        """
