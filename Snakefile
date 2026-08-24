# ===========================================================================
# UV-RIP-seq Snakemake workflow
# Generalization of UV-RIP_rDNA.sh + UVRIP_downsample_BW.sh
# ===========================================================================
import os
import re
import pandas as pd

configfile: "config.yaml"

# ---- Metadata ------------------------------------------------------------
meta = (
    pd.read_csv(config["metadata"], sep="\t", dtype=str)
      .fillna("")
      .apply(lambda c: c.str.strip())
)
meta = meta.set_index("sample", drop=False)

SAMPLES    = list(meta["sample"])
PE_SAMPLES = [s for s in SAMPLES if meta.loc[s, "fastq_r2"]]
SE_SAMPLES = [s for s in SAMPLES if not meta.loc[s, "fastq_r2"]]
CONDITIONS = list(dict.fromkeys(meta["condition"]))
GROUPS     = list(dict.fromkeys(zip(meta["condition"], meta["type"])))  # (condition, type)
CONTROL    = config["control_type"]
OUT        = config["outdir"]


def members(condition, typ):
    """Sample IDs in a (condition, type) group, ordered by replicate."""
    sub = meta[(meta["condition"] == condition) & (meta["type"] == typ)].copy()
    sub["_rep"] = pd.to_numeric(sub["replicate"], errors="coerce")
    return list(sub.sort_values("_rep")["sample"])


def is_paired(sample):
    return bool(meta.loc[sample, "fastq_r2"])


def constr(names):
    """Regex matching exactly this set of sample names (never-match if empty)."""
    return "|".join(re.escape(n) for n in names) if names else r"(?!x)x"


def p(*parts):
    return os.path.join(OUT, *parts)


# Conditions where an IP-over-control track can be built.
BW_CONDS = [c for c in CONDITIONS if members(c, "IP") and members(c, CONTROL)]

wildcard_constraints:
    sample    = constr(SAMPLES),
    condition = "|".join(re.escape(c) for c in CONDITIONS) or r"(?!x)x",
    type      = "IP|IgG|Input",

# ---- Targets -------------------------------------------------------------
def all_targets():
    t = [p("BIGWIG", f"UVRIP_{c}_IP_over_{CONTROL}_log2ratio.bw") for c in BW_CONDS]
    if config["bigwig"]["make_coverage_tracks"]:
        norm = config["bigwig"]["normalize_using"]
        t += [p("BIGWIG", f"UVRIP_{c}_{ty}_{norm}.bw") for (c, ty) in GROUPS]
    t += [p("qc", "multiqc_report.html")]
    if config["de"]["enabled"]:
        t += [p("de", "deseq2.done")]
    return t

rule all:
    input:
        all_targets()

# ---- Rule modules --------------------------------------------------------
include: "workflow/rules/process.smk"
include: "workflow/rules/downsample.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/diffenrich.smk"
