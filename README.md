# UV-RIP-seq Pipeline (Snakemake)

A Snakemake workflow for UV-crosslinked RNA immunoprecipitation sequencing (UV-RIP-seq), from raw FASTQ through depth-matched, replicate-merged **log2(IP / control) bigWig tracks**, with an optional DESeq2 differential-enrichment step.

It is a generalization of two lab scripts (`UV-RIP_rDNA.sh` and `UVRIP_downsample_BW.sh`): the hardcoded sample arrays, hand-computed downsample fractions, merge groupings, and absolute paths all move into `config.yaml` + `metadata.tsv`.

> **Design note for Quest.** This workflow is set up to run inside **one activated conda env, without `--use-conda`**. The SLURM profile passes `--export=ALL`, so compute nodes inherit that env's `PATH` and every rule finds its tools. This deliberately avoids the failure modes we've hit before: `--use-conda` silently failing in non-interactive `/usr/bin/bash` subshells (`conda: command not found`), and `module load snakemake` running inside a Singularity container without SLURM access (exit 127). The `conda:` directives are still present on every rule, so `--use-conda` remains available if you ever want per-rule envs.

---

## Overview

| Stage | Rules | Tools | Purpose |
|-------|-------|-------|---------|
| Process | `trim_{pe,se}`, `star_{pe,se}`, `unique_filter`, `dedup` | Trim Galore, STAR, samtools | Trim → STAR align to hg38+rDNA (unique mappers, `EndToEnd`) → MAPQ filter → PE `markdup` dedup → index |
| Depth-match + tracks | `downsample_plan`, `subsample`, `merge_group`, `bigwig_log2`, `bigwig_coverage` | samtools, deepTools | Auto-compute downsample fractions to a common depth → subsample → merge reps per (condition, type) → `bamCompare` log2(IP/control) (+ optional coverage) |
| QC | `multiqc` | MultiQC | Aggregate FastQC/trim/STAR/dedup reports |
| DE *(optional)* | `featurecounts`, `tidy_counts`, `deseq2` | Subread, DESeq2 | IP-vs-control enrichment per condition; interaction test + delta-log2FC scatter when two conditions exist (e.g. ISD vs DMSO) |

The downsampling step reproduces automatically what was previously done by hand (`samtools view -s 42.60172 …`): it counts mapped/primary/deduplicated reads per library and subsamples each to a shared depth so the IP/control ratio isn't driven by library size. Fractions and the exact `-s` arguments are written to `BAMS/downsample_plan.tsv` as an audit trail.

---

## Repository structure

```
uvrip-seq-snakemake/
├── Snakefile                       # loads config/metadata, builds targets
├── config.yaml                     # all params/toggles — edit this
├── metadata.tsv                    # sample sheet — edit this
├── Run_pipeline.sh                 # robust entry point (portable conda + SLURM)
├── profiles/slurm/config.yaml      # SLURM profile (edit account/partition)
└── workflow/
    ├── rules/
    │   ├── process.smk             # trim / align / filter / dedup
    │   ├── downsample.smk          # plan / subsample / merge / bigwig
    │   ├── qc.smk                  # multiqc
    │   └── diffenrich.smk          # optional featureCounts + DESeq2
    ├── envs/
    │   ├── uvrip.yaml              # Snakemake + all NGS tools (one env)
    │   └── deseq2.yaml             # optional R/DESeq2 env
    └── scripts/
        ├── compute_downsample.py   # fractions to a common depth
        └── de_analysis.R           # DESeq2 differential enrichment
```

---

## Install

```bash
mamba env create -f workflow/envs/uvrip.yaml
#   or into a project allocation:
mamba env create -f workflow/envs/uvrip.yaml -p /projects/b1042/LauberthLab/conda-envs/uvrip
```

This env includes Snakemake **and** every tool the rules call, so a single `conda activate uvrip` is all that's needed. If you'd rather keep Snakemake in your existing `snake` env, delete the `snakemake-minimal` line from `uvrip.yaml`, run with `--use-conda` (the `conda:` directives will build per-rule envs), and let `Run_pipeline.sh` handle the conda-on-Quest wiring.

The `deseq2.yaml` env is only needed for the optional DE step.

---

## Sample sheet (`metadata.tsv`)

Tab-separated, one row per library:

| Column | Required | Notes |
|--------|----------|-------|
| `sample` | yes | Unique ID; final BAM is `nodup_<sample>.bam` |
| `condition` | yes | e.g. `TOP1`, or `DMSO`/`ISD` for a drug contrast |
| `type` | yes | `IP`, `IgG`, or `Input`; the bamCompare denominator is `control_type` |
| `replicate` | yes | Integer; sets merge order |
| `fastq_r1` | yes | Absolute path to R1 |
| `fastq_r2` | no | Path to R2; **blank ⇒ single-end** (SE rule variants kick in automatically) |
| `downsample_frac` | no | Manual override of the computed fraction (0–1); blank ⇒ auto |

Pre-filled with the four samples from your original run (TOP1 IP rep1/2, IgG rep1/2). Add a second `condition` (e.g. ISD/DMSO) and the workflow automatically builds a track per condition and runs the interaction test.

---

## Configuration (`config.yaml`)

Highlights (STAR/trim defaults reproduce `UV-RIP_rDNA.sh` exactly):

- `outdir` — output root (default `results`).
- `star_index`, `gtf` — references; the `star=2.7.5a` pin in `uvrip.yaml` must match whatever built the index.
- `remove_duplicates` — `true` runs the PE `fixmate -m → sort → markdup -r` path.
- `downsample.mode` — `min` (equalize to smallest library) or `target` (`downsample_reads`); `group_by: all|type`; `seed` (default 42).
- `control_type` — bamCompare denominator (`IgG`/`Input`).
- `bigwig.*` — `normalize_using` (RPKM), `scale_factors_method` (`None`, trusting the explicit downsampling as in the original), `bin_size` (10), `make_coverage_tracks`.
- `de.*` — the optional step, `enabled: false` by default.

SLURM account/partition live in **`profiles/slurm/config.yaml`** (`default-resources`), not in `config.yaml` — edit them there.

---

## Usage

```bash
conda activate uvrip

# Dry runs
bash Run_pipeline.sh --dry-run             # local DAG preview
bash Run_pipeline.sh --slurm --dry-run     # SLURM DAG preview (uses the profile)

# Real runs
bash Run_pipeline.sh --slurm               # submit to SLURM (recommended on Quest)
bash Run_pipeline.sh                        # run locally with $CORES (default 20)

# If a previous run left a lock:
bash Run_pipeline.sh --unlock

# Pass extra Snakemake args through after `--`:
bash Run_pipeline.sh --slurm -- --until dedup       # stop after dedup
bash Run_pipeline.sh --slurm -- results/BIGWIG/UVRIP_TOP1_IP_over_IgG_log2ratio.bw   # one target
```

`Run_pipeline.sh` sources conda portably (via `$CONDA_EXE`, with fallbacks), `module purge`s to stop system tools shadowing the env, then runs Snakemake locally or via the SLURM profile.

You can also call Snakemake directly once the env is active:

```bash
snakemake --profile profiles/slurm --rerun-incomplete
```

---

## Outputs

```
results/
├── TRIMMED/                                    # trimmed FASTQs + FastQC/trim reports
├── BAMS/
│   ├── nodup_<sample>.bam(.bai)                # unique, deduplicated
│   ├── nodup_<sample>_sub.bam(.bai)            # depth-matched (copy if frac==1)
│   └── downsample_plan.tsv                     # reads, target, fraction, -s arg
├── MERGED/  UVRIP_<condition>_<type>_merged.bam(.bai)
├── BIGWIG/
│   ├── UVRIP_<condition>_IP_over_<control>_log2ratio.bw   # main track
│   └── UVRIP_<condition>_<type>_RPKM.bw                   # per-group coverage (optional)
├── qc/  multiqc_report.html
├── counts/  de/                                # optional DE step
└── logs/  (+ logs/slurm/)
```

---

## Differential enrichment (optional)

Set `de.enabled: true`, provide a `gtf`, and (if not using the bundled env) create `deseq2.yaml`. `de_analysis.R`:

- **Per condition:** DESeq2 `~type` → IP-vs-control enrichment (`<cond>_IP_vs_<control>_{all,significant}.tsv`).
- **Two conditions (e.g. ISD vs DMSO):** an interaction model `~condition + type + condition:type` whose interaction term is the *differential* enrichment, plus `delta_log2FC_scatter.pdf`.

Confirm `de.strandedness` (default `2`, reverse) against your library prep.

---

## Troubleshooting

- **`conda: command not found` in a SLURM job.** Use `Run_pipeline.sh` (it sources `conda.sh` from `$CONDA_EXE`) and keep `--export=ALL` in the profile so the env's PATH reaches compute nodes. With the single-env setup you should not need `--use-conda` at all.
- **`module load snakemake` fails with exit 127.** Don't use the Quest Snakemake module (Singularity, no SLURM access). Use the conda `uvrip` env.
- **System deeptools/Python shadows the env.** `Run_pipeline.sh` runs `module purge` after activating conda; if calling Snakemake directly, `module purge` first.
- **`samtools markdup` missing / `fixmate` rejects `-m`.** Needs samtools ≥ 1.7 (env pins 1.19.2).
- **STAR OOM / "genome version mismatch."** hg38+rDNA index needs ~32–40 GB (the `star` rules request 40 GB); the `star=2.7.5a` pin must match the index build.
- **A track is missing for a condition.** `bigwig_log2` only runs where both an `IP` and a `control_type` merged BAM exist — check the sample sheet.
- **Reserved-name errors when editing rules.** `params` names like `remove`, `input`, `output`, etc. are reserved by Snakemake (that's why the dedup rule uses `rmdup`).
