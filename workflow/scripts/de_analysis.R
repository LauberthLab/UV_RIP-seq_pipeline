#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# UV-RIP differential enrichment.
#
#  * Per condition, models IP vs control (IgG/Input) enrichment with DESeq2.
#  * If two conditions are present (e.g. ISD vs DMSO), also fits an interaction
#    model (~ condition * type) whose interaction term is the differential
#    enrichment between conditions, and draws a delta-log2FC scatter
#    (log2FC_condB enrichment vs log2FC_condA enrichment).
#
# Counts matrix: first column "Geneid", remaining columns = sample IDs matching
# metadata$sample. Input-type samples are ignored.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse); library(DESeq2); library(edgeR)
  library(ggplot2); library(ggrepel); library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts",   type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--outdir",   type = "character"),
  make_option("--control",  type = "character", default = "IgG"),
  make_option("--min-cpm",  type = "double",    default = 0.3,    dest = "min_cpm"),
  make_option("--lfc",      type = "double",    default = 0.5849),
  make_option("--padj",     type = "double",    default = 0.05)
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## ---- Load ----------------------------------------------------------------
cnt  <- as.data.frame(fread(opt$counts))
rownames(cnt) <- cnt$Geneid; cnt$Geneid <- NULL
meta <- as.data.frame(fread(opt$metadata))
meta <- meta[meta$type != "Input", ]                 # drop Input
meta <- meta[meta$sample %in% colnames(cnt), ]
cnt  <- cnt[, meta$sample, drop = FALSE]
meta$type      <- factor(meta$type, levels = c(opt$control,
                     setdiff(unique(meta$type), opt$control)))  # control = ref
meta$condition <- factor(meta$condition)

cpm_ok <- function(mat, min_cpm, k) rowSums(edgeR::cpm(mat) >= min_cpm) >= k

run_deseq <- function(sub_meta, design, coef = NULL, coef_pattern = NULL, tag) {
  m <- cnt[, sub_meta$sample, drop = FALSE]
  keep <- cpm_ok(m, opt$min_cpm, max(2, floor(ncol(m) / 2)))
  m <- m[keep, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(round(m), sub_meta, design)
  dds <- DESeq(dds, quiet = TRUE)
  if (is.null(coef) && !is.null(coef_pattern)) {
    coef <- grep(coef_pattern, resultsNames(dds), value = TRUE)[1]
    if (is.na(coef)) stop("no coefficient matching: ", coef_pattern)
  }
  res <- results(dds, name = coef)
  out <- as.data.frame(res); out$Geneid <- rownames(out)
  out <- out[order(out$padj), c("Geneid", setdiff(colnames(out), "Geneid"))]
  fwrite(out, file.path(opt$outdir, sprintf("%s_all.tsv", tag)), sep = "\t")
  sig <- subset(out, !is.na(padj) & padj < opt$padj & abs(log2FoldChange) > opt$lfc)
  fwrite(sig, file.path(opt$outdir, sprintf("%s_significant.tsv", tag)), sep = "\t")
  message(sprintf("%s: %d tested, %d significant", tag, nrow(out), nrow(sig)))
  out
}

## ---- Per-condition IP vs control ----------------------------------------
enrich <- list()
for (cond in levels(meta$condition)) {
  sm <- meta[meta$condition == cond, ]
  if (!all(c(opt$control, "IP") %in% sm$type)) next
  ip_lvl <- setdiff(levels(meta$type), opt$control)[1]
  enrich[[cond]] <- run_deseq(
    sm, design = ~type,
    coef = paste0("type_", ip_lvl, "_vs_", opt$control),
    tag  = sprintf("%s_IP_vs_%s", cond, opt$control))
}

## ---- Two conditions: interaction + delta-log2FC scatter ------------------
conds <- levels(meta$condition)
if (length(conds) == 2 && all(c(opt$control, "IP") %in% meta$type)) {
  a <- conds[1]; b <- conds[2]
  meta$condition <- relevel(meta$condition, ref = a)
  inter <- run_deseq(
    meta, design = ~condition + type + condition:type,
    coef_pattern = "condition.*\\.type",   # the condition:type interaction term
    tag = sprintf("interaction_%s_vs_%s_IPenrichment", b, a))

  ea <- enrich[[a]]; eb <- enrich[[b]]
  if (!is.null(ea) && !is.null(eb)) {
    m <- merge(ea[, c("Geneid","log2FoldChange")],
               eb[, c("Geneid","log2FoldChange")],
               by = "Geneid", suffixes = c(paste0("_",a), paste0("_",b)))
    colnames(m) <- c("Geneid","lfc_a","lfc_b")
    sig_ids <- inter$Geneid[!is.na(inter$padj) & inter$padj < opt$padj]
    m$sig <- m$Geneid %in% sig_ids
    p <- ggplot(m, aes(lfc_a, lfc_b, color = sig)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
      geom_point(size = 1, alpha = 0.6) +
      scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
      labs(x = sprintf("log2 IP/%s enrichment (%s)", opt$control, a),
           y = sprintf("log2 IP/%s enrichment (%s)", opt$control, b),
           title = "UV-RIP differential enrichment") +
      theme_bw()
    ggsave(file.path(opt$outdir, "delta_log2FC_scatter.pdf"), p, width = 5.5, height = 5)
    message("wrote delta_log2FC_scatter.pdf")
  }
}
message("DESeq2 differential enrichment complete.")
