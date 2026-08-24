# ---------------------------------------------------------------------------
# Aggregated QC. Depends on all nodup BAMs so it runs after processing.
# ---------------------------------------------------------------------------
rule multiqc:
    input:
        expand(p("BAMS", "nodup_{s}.bam"), s=SAMPLES),
    output:
        p("qc", "multiqc_report.html"),
    params:
        scan   = lambda wc: " ".join([p("TRIMMED"), p("BAMS"), p("logs")]),
        outdir = p("qc"),
    resources:
        mem_mb = 8000, runtime = 60
    log:
        p("logs", "multiqc.log")
    conda:
        "../envs/uvrip.yaml"
    shell:
        r"""
        multiqc {params.scan} -o {params.outdir} -n multiqc_report.html -f > {log} 2>&1
        """
