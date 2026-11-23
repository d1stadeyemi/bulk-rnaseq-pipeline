#!/usr/bin/env Rscript
###############################################################################
# deseq2_04_results.R   --  STAGE 4 of 4:  CONTRAST, SHRINKAGE, DEG TABLES
#
# Extracts the contrast from the fitted model and turns it into tables and
# figures. This is the ONLY stage that reads --alpha / --lfc-report, so you can
# re-run it alone as many times as you like without repeating the fit:
#
#     Rscript deseq2_04_results.R <same flags> --alpha 0.01 --lfc-report 1
#
#  ==========================================================================
#  TWO THRESHOLDS, TWO DIFFERENT JOBS   (the single most common confusion)
#  ==========================================================================
#  --lfc-report    0.585   DEG-CALLING cutoff. A feature is a DEG when
#                          padj < --alpha AND its fold change clears this.
#                          0.585 = log2(1.5) = a 1.5-fold change.   <-- DEFAULT
#  --lfc-threshold 0       TESTING threshold, i.e. the null handed to
#                          results(lfcThreshold=). 0 is DESeq2's own default
#                          and means the familiar H0: beta = 0.       <-- DEFAULT
#  --alpha         0.05    FDR cutoff on padj.                        <-- DEFAULT
#
#  Setting --lfc-threshold > 0 changes what the p-value MEANS: it becomes a test
#  of |beta| > x rather than beta != 0. That is statistically stricter and the
#  vignette prefers it to post-hoc filtering, but the resulting padj is NOT
#  comparable to any other DESeq2 output your collaborators will have. So the
#  default is the familiar one -- test against zero, report at 1.5-fold -- and
#  the strict route is one flag away. Whichever you pick, everything downstream
#  (shrinkage, MA, volcano, DEG calls, s-values) stays consistent with it.
#
#  Both thresholds are swept regardless of your choice, in
#  22_threshold_sensitivity_*, so the robustness of the DEG count is always on
#  disk with the active setting flagged.
#  CAVEAT: independent filtering is optimised for the --alpha actually passed
#  to results() (vignette: filtering optimizes "the number of genes which will
#  have an adjusted p value below a given FDR cutoff, alpha"). Rows at OTHER
#  alpha values reuse that padj and are therefore close approximations, not
#  exact. For an exact count at another alpha, re-run with that --alpha.
#
#  ==========================================================================
#  --alt-hypothesis AND WHY IT CONSTRAINS --lfc-threshold
#  ==========================================================================
#  Vignette definitions, where x is lfcThreshold:
#      greaterAbs  |beta| > x   two-tailed                        <-- DEFAULT
#      lessAbs     |beta| < x   p = max of the upper and lower tests
#      greater      beta  > x
#      less         beta  < -x
#
#  lessAbs is an EQUIVALENCE test -- "show me genes that did NOT change much" --
#  and it is meaningless without a threshold, because |beta| < 0 is the empty
#  set. This stage therefore HARD-STOPS on `--alt-hypothesis lessAbs` when
#  --lfc-threshold is 0, rather than writing an empty or nonsensical DEG table.
#
#  greater/less are one-sided and ARE valid at lfcThreshold = 0. The
#  altHypothesis is forwarded to results() whenever it is not greaterAbs, so the
#  reported FDR is the correct one-sided FDR for the selected tail. Filtering a
#  tail out of a two-sided p-value afterwards would not give that.
#
#  ==========================================================================
#  SHRINKAGE:  --shrink, --rank-on, AND s-VALUES
#  ==========================================================================
#  --shrink <apeglm|ashr|normal|none>   Default apeglm.
#      apeglm and ashr use wide-tailed priors that preserve large LFCs; normal
#      is the original 2014 estimator and over-shrinks them.
#      apeglm needs a COEFFICIENT (not an arbitrary contrast); ashr accepts a
#      contrast; normal cannot shrink interaction terms. All three accept a
#      positive lfcThreshold (vignette, "Extended section on shrinkage
#      estimators"); with a threshold, apeglm and ashr return FSOS s-values.
#      Fallback chain:  <--shrink> -> apeglm -> ashr -> normal
#      normal is last because it can over-shrink large LFCs and is not
#      compatible with s-values. A fallback is always logged with the reason.
#
#  --rank-on <mle|shrunk>   Default mle. Which log2FC calls and ranks DEGs and
#      draws the volcano. Whichever you choose is recorded in the output as
#      `lfc_used_source`, so the table can never be ambiguous about it.
#
#  s-VALUES. With a positive lfcThreshold, apeglm and ashr return s-values
#  giving the "false sign or small" (FSOS) rate -- the rate at which a gene's
#  LFC is either the wrong sign or smaller than the threshold. In exactly that
#  configuration the s-value is the more defensible statistic, so it is carried
#  into the results table as `svalue` rather than discarded. It is reported
#  alongside padj, not silently substituted for it.
#
#  CONSISTENCY DETAILS THAT ARE EASY TO GET WRONG, AND ARE HANDLED HERE
#  --------------------------------------------------------------------
#  * lfcShrink is called with res = res. lfcShrink does not compute its own
#    p-values -- it recycles them -- and without `res` it rebuilds a results
#    table internally at DESeq2's DEFAULT alpha of 0.1. Passing res forwards
#    your independent-filtering and alpha choices into the shrunken table.
#  * plotMA takes alpha explicitly. plotMA.DESeqResults defaults to 0.1
#    regardless of what you passed to results(), so without this the raw MA plot
#    and the shrunken MA plot in the same folder would colour points at two
#    different cutoffs with nothing saying so.
#
#  ==========================================================================
#  neg_log10_pvalue / neg_log10_padj  --  WHY padj = 0 IS NOT AN ERROR
#  ==========================================================================
#  DESeq2 does not cap or round. p = 2*pnorm(-|stat|) for a Wald statistic of
#  -53 is roughly 1e-615, and the smallest number an IEEE double can hold is
#  ~2.2e-308, so 0 is the only representable answer. `padj = 0` in the results
#  table means "below double precision", not zero.
#
#  The quantity itself is recoverable: pnorm(log.p = TRUE) returns log p without
#  underflowing, and Benjamini-Hochberg is a product, so it is exact in log
#  space too. Two extra columns therefore carry the true values:
#      neg_log10_pvalue   -log10(p),    finite however large
#      neg_log10_padj     -log10(padj), finite however large
#  They are self-checked against DESeq2's own padj on every gene DESeq2 could
#  compute, and the maximum log-unit deviation is logged. If it exceeds 0.05 the
#  columns are left NA rather than trusted.
#
#  Available only for a two-sided Wald test (altHypothesis greaterAbs) at
#  lfcThreshold = 0 -- the LRT statistic is chi-square, a threshold test uses a
#  different null, and greater/less are one-sided, so `stat` would not map to
#  the p-value this way. Recorded in the run summary as
#  exact_neglog10_available.
#
#  These columns also fix a RANKING bug: every underflowed gene has padj == 0,
#  so ordering by padj left the most significant genes tied and broke the tie by
#  |LFC| -- which is what selects the top-N heatmap, the count plots and the
#  volcano labels.
#
#  ==========================================================================
#  OTHER PARAMETERS
#  ==========================================================================
#  --no-independent-filtering   Disable results() independent filtering.
#  --top-n <INT>                Default 30. Count plots and top-DEG heatmap.
#  --volcano-y <padj|pvalue>    Default padj. The y-axis is then the same
#                               quantity significance is called on, so the
#                               dotted line IS the boundary and nothing above it
#                               can be grey; genes with padj = NA (dropped by
#                               independent filtering, never ranked) are omitted
#                               and the count is logged. `pvalue` gives the
#                               limma/EnhancedVolcano convention: finer
#                               resolution at the top and no BH step-function
#                               banding, but the line is then only the raw p at
#                               which padj crossed alpha, and filtered genes can
#                               appear above it.
###############################################################################

.args_raw <- commandArgs(trailingOnly = TRUE)
.self <- (function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE) else NA_character_
})()
source(file.path(dirname(.self), "deseq2_common.R"))

opt <- parse_args(.args_raw)
.unknown <- setdiff(names(opt), .known_keys)
if (length(.unknown) > 0) { message("ERROR: unknown option(s): ",
  paste0("--", .unknown, collapse = ", ")); quit(status = 2) }
g <- mk_getters(opt)

outdir <- g$chr("outdir", "deseq2_out")
log_open(outdir, "4/results")
logbanner("STAGE 4: CONTRAST, SHRINKAGE, DEG TABLES")

## ---------------------------------------------------------------------------
## 4.1  Thresholds and their coherence
##      Flag-only checks run BEFORE the upstream state is loaded, so a typo is
##      reported as a typo rather than as a signature mismatch.
## ---------------------------------------------------------------------------
alpha    <- g$num("alpha", 0.05, min = 1e-10, max = 0.999999)
lfc_rep  <- g$num("lfc-report", 0.585, min = 0)
lfc_thr  <- g$num("lfc-threshold", 0, min = 0)
alt_hyp  <- g$chr("alt-hypothesis", "greaterAbs")
rank_on  <- tolower(g$chr("rank-on", "mle"))
shrink   <- g$chr("shrink", "apeglm")
indep_f  <- !g$bool("no-independent-filtering")
top_n    <- g$int("top-n", 30L, min = 1)

# padj by default. The y-axis is then the same quantity significance is called
# on, so the cutoff line IS the boundary and nothing above it can be grey.
# Genes with padj = NA (removed by independent filtering, never ranked) are
# dropped rather than drawn at a raw-p height they were never judged at.
# `pvalue` reproduces the limma / EnhancedVolcano convention: more resolution at
# the top and no BH step-function banding, but the line is then only the raw p
# at which padj crossed alpha, and filtered genes can appear above it.
volcano_y <- tolower(g$chr("volcano-y", "padj"))
if (!volcano_y %in% c("padj","pvalue")) die("--volcano-y must be padj or pvalue.")

if (!alt_hyp %in% c("greaterAbs","lessAbs","greater","less")) die("--alt-hypothesis invalid: ", alt_hyp)
if (!rank_on %in% c("mle","shrunk"))                          die("--rank-on invalid: ", rank_on)
if (!shrink  %in% c("apeglm","ashr","normal","none"))         die("--shrink invalid: ", shrink)
if (rank_on == "shrunk" && shrink == "none")
  die("--rank-on shrunk requires a --shrink estimator (you set --shrink none).")

# lessAbs is an equivalence test and is undefined without a threshold.
if (alt_hyp == "lessAbs" && lfc_thr <= 0)
  die("--alt-hypothesis lessAbs tests |LFC| < x, which is the EMPTY SET when\n",
      "       --lfc-threshold is 0 (the default).\n",
      "       lessAbs is an equivalence test: 'show me genes that did NOT change by\n",
      "       more than x'. It needs an explicit x, e.g.\n",
      "         --alt-hypothesis lessAbs --lfc-threshold 0.585\n",
      "       If you meant the usual 'genes that DID change', drop --alt-hypothesis.")

force <- g$bool("force")
st1 <- load_state(outdir, 1, opt, force = force)
st2 <- load_state(outdir, 2, opt, force = force)
st3 <- load_state(outdir, 3, opt, force = force)

control <- st1$control; treat <- st1$treat; contrast_nm <- st1$contrast_nm
design_formula <- st1$design_formula; test_type <- st3$test_type

if (test_type == "lrt") {
  if (alt_hyp == "lessAbs")
    die("--alt-hypothesis lessAbs cannot be combined with --test lrt.\n",
        "       lessAbs is an equivalence test on a COEFFICIENT (|LFC| < x). The LRT\n",
        "       tests whether a set of terms improves the model at all, so there is no\n",
        "       coefficient to bound and lfcThreshold does not apply.\n",
        "       Use --test wald --alt-hypothesis lessAbs --lfc-threshold <x> instead.")
  if (alt_hyp != "greaterAbs") {
    logmsg("NOTE: --alt-hypothesis ", alt_hyp, " is a Wald-test concept and does not ",
           "apply to the LRT; ignoring it for the test. The DEG call still uses ",
           "--lfc-report ", lfc_rep, " on the reported contrast.")
  }
  if (lfc_thr > 0) {
    logmsg("NOTE: --lfc-threshold does not apply to the LRT (it compares nested models ",
           "rather than testing a coefficient against a threshold). Setting it to 0; ",
           "--lfc-report ", lfc_rep, " still governs the DEG call.")
    lfc_thr <- 0
  }
}
if (lfc_thr > 0) {
  if (alt_hyp == "lessAbs") {
    if (lfc_rep > lfc_thr || lfc_rep == 0) {
      logmsg("NOTE: lessAbs tests |LFC| < ", lfc_thr, "; coercing --lfc-report ", lfc_rep, " -> ", lfc_thr)
      lfc_rep <- lfc_thr }
  } else if (lfc_rep < lfc_thr) {
    logmsg("NOTE: --lfc-report (", lfc_rep, ") is below --lfc-threshold (", lfc_thr,
           "); coercing up to ", lfc_thr, " (genes below the testing null cannot be significant).")
    lfc_rep <- lfc_thr
  }
}
lfc_pass <- make_lfc_pass(alt_hyp, lfc_rep)

rlib <- g$chr("rlib")
if (!is.null(rlib)) { if (!dir.exists(rlib)) die("--rlib not found: ", rlib); .libPaths(c(rlib, .libPaths())) }
suppressWarnings(suppressMessages({ library(DESeq2); if (has("ggplot2")) library(ggplot2) }))
set.seed(st1$seed)

logmsg("thresholds: alpha=", alpha, "  lfc-report=", lfc_rep, " (", sprintf("%.2f", 2^lfc_rep),
       "-fold)  lfc-threshold(test)=", lfc_thr, "  altHypothesis=", alt_hyp, "  rank-on=", rank_on)

## ===========================================================================
## RESULTS ENGINE -- once per level
## ===========================================================================
results_level <- function(L) {
  dds <- L$dds; ann <- L$ann; coldata <- L$coldata
  qc_obj <- L$qc_obj; qc_lab <- L$qc_lab; vsd <- L$vsd; norm_cts <- L$norm_cts
  io <- make_io(L$level_dir)
  logmsg("---- level: ", L$level_lab, " ----")

  ## --- the test -----------------------------------------------------------
  ra <- list(object = dds, contrast = c("condition", treat, control),
             alpha = alpha, independentFiltering = indep_f)
  # Forward altHypothesis whenever it is not the default. greater/less are valid
  # at lfcThreshold = 0 and give the correct ONE-SIDED FDR, which post-hoc
  # tail-filtering of a two-sided p-value would not.
  if (test_type == "lrt") {
    # Nothing to forward: the LRT p-value comes from the model comparison, not
    # from a coefficient, so lfcThreshold/altHypothesis do not apply.
  } else if (lfc_thr > 0) {
    ra$lfcThreshold <- lfc_thr; ra$altHypothesis <- alt_hyp
  } else if (alt_hyp != "greaterAbs") {
    ra$altHypothesis <- alt_hyp
    logmsg("one-sided test requested (altHypothesis=", alt_hyp,
           ") at lfcThreshold=0; results() returns the one-sided FDR.")
  }
  res <- do.call(results, ra)

  ## --- shrinkage ----------------------------------------------------------
  rn <- resultsNames(dds)
  want_coef <- paste0("condition_", treat, "_vs_", control)
  coef_idx <- match(want_coef, rn)

  # Fallback order: requested estimator first, then the wide-tailed priors
  # (apeglm, ashr), then normal last. Vignette: "All three methods allow for a
  # positive lfcThreshold to be specified"; with a threshold, apeglm and ashr
  # return s-values (FSOS). normal is last because it can over-shrink large
  # LFCs and "is not really compatible with computing s-values".
  order_try <- unique(c(shrink, "apeglm", "ashr", "normal"))

  resShrink <- NULL; shrink_used <- "none"
  if (shrink != "none") {
    attempt <- function(type) {
      a <- list(dds = dds, type = type, res = res)   # res= keeps alpha + filtering
      if (lfc_thr > 0) { a$lfcThreshold <- lfc_thr; a$svalue <- TRUE }
      if (type == "apeglm") {
        if (is.na(coef_idx))
          stop("coefficient '", want_coef, "' is not in resultsNames(): ",
               paste(rn, collapse = ", "), "; apeglm shrinks a coefficient, not a contrast")
        a$coef <- coef_idx
      } else a$contrast <- c("condition", treat, control)
      suppressWarnings(do.call(lfcShrink, a))
    }
    for (ty in order_try) {
      resShrink <- tryCatch(attempt(ty), error = function(e) {
        logmsg("NOTE: shrinkage '", ty, "' unavailable (", conditionMessage(e), ")"); NULL })
      if (!is.null(resShrink)) { shrink_used <- ty; break }
    }
    if (is.null(resShrink)) logmsg("WARNING: all shrinkage estimators failed; MLE only.")
    else if (shrink_used != shrink)
      logmsg("NOTE: fell back from '", shrink, "' to '", shrink_used, "' shrinkage.")
    else logmsg("shrinkage: ", shrink_used, if (lfc_thr > 0) "  (with s-values)" else "")
  }
  if (rank_on == "shrunk" && is.null(resShrink))
    die("--rank-on shrunk requested but no shrinkage estimator could be fit.")

  ## --- assemble ------------------------------------------------------------
  res_df <- as.data.frame(res); res_df$gene_id <- rownames(res_df)
  if (!is.null(resShrink)) {
    sh <- as.data.frame(resShrink)
    add <- data.frame(gene_id = rownames(sh),
                      log2FoldChange_shrunk = sh$log2FoldChange, stringsAsFactors = FALSE)
    # s-values exist only when a threshold was used; carry them, do not discard.
    if ("svalue" %in% colnames(sh)) add$svalue <- sh$svalue
    res_df <- merge(res_df, add, by = "gene_id", all.x = TRUE)
  }
  if (!is.null(ann)) res_df <- merge(res_df, ann, by = "gene_id", all.x = TRUE)

  eff_col <- if (rank_on == "shrunk") "log2FoldChange_shrunk" else "log2FoldChange"
  if (!eff_col %in% colnames(res_df)) eff_col <- "log2FoldChange"
  res_df$lfc_used <- res_df[[eff_col]]; res_df$lfc_used_source <- eff_col

  ## --- exact -log10 p and padj, immune to double underflow -----------------
  # DESeq2 is NOT capping or rounding. p = 2*pnorm(-|stat|) for stat = -53 is
  # about 1e-615, and the smallest number an IEEE double can hold is ~2.2e-308,
  # so 0 is the only representable answer. The QUANTITY is fine, though:
  # pnorm(log.p = TRUE) returns log p directly without underflowing, and BH is a
  # product, so it is exact in log space too. Recovered here for the results
  # TABLE, the ordering, and the volcano.

  res_df$neg_log10_pvalue <- NA_real_
  res_df$neg_log10_padj   <- NA_real_
  exact_ok <- FALSE
  n_uf <- sum(res_df$padj == 0, na.rm = TRUE)
  # Two-sided Wald only: log(2) + pnorm(-|stat|) is the greaterAbs p-value.
  # greater/less are one-sided tests (vignette: "greater - beta > x",
  # "less - beta < -x"), so this formula would not reproduce their p-values.
  if (identical(test_type, "wald") && lfc_thr == 0 && alt_hyp == "greaterAbs" &&
      "stat" %in% colnames(res_df)) {
    tested <- !is.na(res_df$padj) & is.finite(res_df$stat)
    if (any(tested)) {
      lp <- log(2) + pnorm(-abs(res_df$stat[tested]), log.p = TRUE)
      m  <- sum(tested); o <- order(lp); ro <- order(o)
      lq <- pmin(0, cummin((lp[o] + log(m) - log(seq_len(m)))[m:1])[m:1])[ro]
      pj  <- res_df$padj[tested]; chk <- pj > 0
      # Must reproduce DESeq2's own padj wherever DESeq2 could compute one.
      dev <- if (any(chk)) max(abs(log(pj[chk]) - lq[chk])) else Inf
      if (is.finite(dev) && dev <= 0.05) {
        res_df$neg_log10_pvalue[tested] <- -lp / log(10)
        res_df$neg_log10_padj[tested]   <- -lq / log(10)
        exact_ok <- TRUE
        # A gene can carry a padj while its Wald statistic is NA (an outlier
        # replaced by DESeq2, say). Leaving its exact column NA would drop it
        # from the volcano silently AND sort it to the bottom of the table, so
        # fall back to the standard value where one exists.
        fb <- !is.na(res_df$padj) & is.na(res_df$neg_log10_padj) & res_df$padj > 0
        if (any(fb)) {
          res_df$neg_log10_padj[fb]   <- -log10(res_df$padj[fb])
          res_df$neg_log10_pvalue[fb] <- -log10(res_df$pvalue[fb])
          logmsg("NOTE: ", sum(fb), " gene(s) have a padj but no usable Wald statistic; ",
                 "their exact columns fall back to the standard value.")
        }
        logmsg("exact -log10 p/padj recovered in log space; agreement with DESeq2 on ",
               sum(chk), " gene(s): max ", signif(dev, 2), " log units.",
               if (n_uf > 0) paste0(" ", n_uf, " gene(s) underflow to padj = 0 in the ",
                 "standard column; the exact column reaches ",
                 round(max(res_df$neg_log10_padj, na.rm = TRUE)), ".") else "")
      } else
        logmsg("NOTE: exact -log10 recovery disagreed with DESeq2 by ", signif(dev, 3),
               " log units; the neg_log10_* columns are left NA.")
    }
  } else if (n_uf > 0)
    logmsg("NOTE: ", n_uf, " gene(s) have padj = 0 (IEEE underflow, not a computed zero). ",
           "Exact recovery requires a two-sided Wald test (greaterAbs) at lfcThreshold = 0.")

  # Rank on the exact value when we have it, so ties at padj = 0 do not decide
  # which genes appear in the top-N figures.
  ord_key <- if (exact_ok) -res_df$neg_log10_padj else res_df$padj
  res_df <- res_df[order(ord_key, -abs(res_df$lfc_used)), ]
  front <- intersect(c("gene_id","gene_name","biotype","chrom","baseMean","log2FoldChange",
                       "log2FoldChange_shrunk","lfc_used","lfc_used_source","lfcSE","stat",
                       "pvalue","padj","neg_log10_pvalue","neg_log10_padj","svalue"),
                     colnames(res_df))
  res_df <- res_df[, c(front, setdiff(colnames(res_df), front))]
  io$table(res_df, paste0("20_results_", contrast_nm, "_full.csv"))

  is_deg <- !is.na(res_df$padj) & res_df$padj < alpha & lfc_pass(res_df$lfc_used)
  sig <- res_df[is_deg, , drop = FALSE]
  io$table(sig, paste0("20_results_", contrast_nm, "_significant.csv"))
  io$table(data.frame(gene_id = rownames(norm_cts), as.data.frame(norm_cts), check.names = FALSE),
           "21_normalized_counts.csv")
  if (!is.null(vsd))
    io$table(data.frame(gene_id = rownames(assay(vsd)), as.data.frame(assay(vsd)), check.names = FALSE),
             "21_vst_matrix.csv")

  n_up <- sum(is_deg & res_df$lfc_used > 0); n_dn <- sum(is_deg & res_df$lfc_used < 0)
  summ <- data.frame(
    metric = c("features_tested","DEGs_total","DEG_up","DEG_down","NA_padj",
               "alpha","lfc_report","lfc_report_fold","lfc_threshold_test","alt_hypothesis",
               "effect_size_column","shrink_used","svalues_available",
               "exact_neglog10_available","n_padj_underflowed","test"),
    value = c(sum(!is.na(res$pvalue)), nrow(sig), n_up, n_dn, sum(is.na(res_df$padj)),
              alpha, lfc_rep, sprintf("%.3f", 2^lfc_rep), lfc_thr, alt_hyp, eff_col, shrink_used,
              "svalue" %in% colnames(res_df), exact_ok, n_uf,
              toupper(test_type)), stringsAsFactors = FALSE)
  io$table(summ, paste0("20_results_", contrast_nm, "_summary.csv"))
  logmsg("DEGs @ padj<", alpha, " & ", alt_hyp, " ", lfc_rep, " (", sprintf("%.2f", 2^lfc_rep),
         "-fold): ", nrow(sig), "  (up=", n_up, ", down=", n_dn, ")")
  if ("svalue" %in% colnames(res_df))
    logmsg("s-values (FSOS: wrong sign OR smaller than ", lfc_thr, ") carried in column `svalue`; ",
           "at s<", alpha, " that set has ", sum(!is.na(res_df$svalue) & res_df$svalue < alpha),
           " features. Reported alongside padj, not substituted for it.")

  ## --- sensitivity ---------------------------------------------------------
  grid <- expand.grid(alpha = c(0.01,0.05,0.10), lfc = c(0,0.585,1,1.5,2), KEEP.OUT.ATTRS = FALSE)
  grid$fold <- round(2^grid$lfc, 3)
  grid$n_DEG <- vapply(seq_len(nrow(grid)), function(i) {
    p <- make_lfc_pass(alt_hyp, grid$lfc[i])
    sum(!is.na(res_df$padj) & res_df$padj < grid$alpha[i] & p(res_df$lfc_used)) }, numeric(1))
  grid$n_up <- vapply(seq_len(nrow(grid)), function(i) {
    p <- make_lfc_pass(alt_hyp, grid$lfc[i])
    sum(!is.na(res_df$padj) & res_df$padj < grid$alpha[i] & p(res_df$lfc_used) &
        res_df$lfc_used > 0) }, numeric(1))
  grid$n_down <- grid$n_DEG - grid$n_up
  grid$is_active <- abs(grid$alpha - alpha) < 1e-12 & abs(grid$lfc - lfc_rep) < 1e-12
  io$table(grid, paste0("22_threshold_sensitivity_", contrast_nm, ".csv"))
  if (has("pheatmap")) {
    # expand.grid varies its FIRST argument fastest, so the vector fills a matrix
    # column-wise with alpha down the rows. Getting this backwards transposes the figure.
    gm <- matrix(grid$n_DEG, nrow = length(unique(grid$alpha)),
                 dimnames = list(paste0("alpha=", unique(grid$alpha)),
                                 paste0("lfc ", alt_hyp, " ", unique(grid$lfc))))
    io$base(paste0("22_threshold_sensitivity_", contrast_nm), function()
      pheatmap::pheatmap(gm, display_numbers = matrix(sprintf("%d", gm), nrow = nrow(gm)),
        number_color = "grey15", fontsize_number = 9, cluster_rows = FALSE, cluster_cols = FALSE,
        main = paste0("DEG count vs thresholds: ", contrast_nm, "\n(effect size = ", eff_col, ")")),
      width = 6, height = 5)
  }

  ## --- DEG by biotype ------------------------------------------------------
  if ("biotype" %in% colnames(res_df) && any(!is.na(res_df$biotype))) {
    tst <- res_df[!is.na(res_df$padj), , drop = FALSE]
    if (nrow(tst) > 0) {
      tst$biotype[is.na(tst$biotype) | !nzchar(tst$biotype)] <- "unknown"
      deg <- tst$padj < alpha & lfc_pass(tst$lfc_used)
      is_up <- deg & tst$lfc_used > 0; is_dn <- deg & tst$lfc_used < 0
      bts <- sort(unique(tst$biotype))
      cnt <- function(v) as.integer(vapply(bts, function(b) sum(v & tst$biotype == b), numeric(1)))
      nt <- as.integer(vapply(bts, function(b) sum(tst$biotype == b), numeric(1)))
      nu <- cnt(is_up); nd <- cnt(is_dn)
      dbt <- data.frame(biotype=bts, n_tested=nt, n_up=nu, n_down=nd,
                        pct_of_tested_DE=round(100*(nu+nd)/pmax(nt,1),2),
                        pct_of_all_up=round(100*nu/max(sum(is_up),1),2),
                        pct_of_all_down=round(100*nd/max(sum(is_dn),1),2), stringsAsFactors=FALSE)
      dbt <- dbt[order(-(dbt$n_up + dbt$n_down)), , drop=FALSE]
      io$table(dbt, paste0("23_DEG_by_biotype_", contrast_nm, ".csv"))
      npc <- c("snRNA","snoRNA","scaRNA","misc_RNA","rRNA","rRNA_pseudogene","Mt_rRNA",
               "vault_RNA","sRNA","ribozyme")
      # Distinct from stage 2's non-polyA %, which is a share of COUNTS and
      # diagnoses library prep. This is a share of UP DEGs and diagnoses the
      # result. Same word, different quantity -- do not merge the two.
      NONPOLYA_DEG_WARN_PCT <- 20
      fr <- 100 * sum(dbt$n_up[dbt$biotype %in% npc]) / max(sum(is_up), 1)
      logmsg("DEG-by-biotype: ", sprintf("%.1f%%", fr), " of UP features are non-polyA ncRNA classes",
             " (threshold ", NONPOLYA_DEG_WARN_PCT, "%)")
      if (fr > NONPOLYA_DEG_WARN_PCT)
        logmsg("WARNING: >20% of UP features are non-polyA ncRNA classes -- consistent with a ",
               "poly-A vs total-RNA library difference rather than biology. ",
               "See tables/02c_library_type_qc.")
    }
  }

  ## --- figures -------------------------------------------------------------
  io$base("30_pvalue_histogram", function() {
    ft <- metadata(res)$filterThreshold
    pv <- if (!is.null(ft) && length(ft)==1 && is.finite(ft))
      res$pvalue[!is.na(res$baseMean) & res$baseMean > ft] else res$pvalue
    hist(pv[!is.na(pv)], breaks=0:20/20, col="grey80", xlab="p-value",
         main=paste0("p-value histogram (passing filter)",
                     if (lfc_thr > 0) "\nNOTE: threshold test -- a spike near 1 is EXPECTED" else "")) })

  if (!is.null(metadata(res)$filterNumRej))
    io$base("31_independent_filtering", function() {
      plot(metadata(res)$filterNumRej, type="b", ylab="# rejections",
           xlab="quantiles of mean-normalized-count filter", main="Independent filtering")
      if (!is.null(metadata(res)$lo.fit)) lines(metadata(res)$lo.fit, col="red")
      if (!is.null(metadata(res)$filterTheta)) abline(v=metadata(res)$filterTheta, lty=2) })

  ma_ylim <- function(x) { m <- max(4, ceiling(quantile(abs(x[is.finite(x)]), .999, na.rm=TRUE))); c(-m,m) }
  # alpha passed EXPLICITLY: plotMA.DESeqResults defaults to 0.1 regardless of
  # what results() was given, which would colour these two panels differently.
  io$base("32_MA_raw", function() plotMA(res, alpha=alpha, ylim=ma_ylim(res$log2FoldChange),
    main=paste0("MA (MLE): ", contrast_nm, "   [padj<", alpha, "]")))
  if (!is.null(resShrink))
    io$base("32_MA_shrunk", function() plotMA(resShrink, alpha=alpha,
      ylim=ma_ylim(resShrink$log2FoldChange),
      main=paste0("MA (", shrink_used, "-shrunk): ", contrast_nm, "   [",
                  if ("svalue" %in% colnames(as.data.frame(resShrink))) "s" else "padj", "<", alpha, "]")))

  if (volcano_y == "padj") {
    volc <- res_df[!is.na(res_df$padj) & !is.na(res_df$lfc_used), , drop = FALSE]
    n_drop <- sum(!is.na(res_df$pvalue) & is.na(res_df$padj))
    if (n_drop > 0)
      logmsg("volcano: ", n_drop, " gene(s) omitted (padj = NA -- removed by ",
             "independent filtering before multiple-testing correction, so never ranked).")
    # Exact values from the assemble step when available. Falling back to a
    # ceiling only when recovery was impossible (LRT, or a threshold test).
    y_cap <- -log10(.Machine$double.xmin)
    volc$at_cap <- FALSE
    if (exact_ok) {
      volc$y_val <- volc$neg_log10_padj
    } else {
      volc$y_val <- -log10(volc$padj)
      volc$at_cap <- volc$padj == 0
      if (any(volc$at_cap)) {
        volc$y_val[volc$at_cap] <- y_cap
        logmsg("volcano: ", sum(volc$at_cap), " gene(s) with padj = 0 drawn on a ceiling at ",
               round(y_cap, 1), "; the exact value could not be recovered for this test type.")
      }
    }
    y_lab <- "-log10(padj)"
    y_fdr <- -log10(alpha)          # exact boundary on this axis
  } else {
    volc <- res_df[!is.na(res_df$pvalue) & !is.na(res_df$lfc_used), , drop = FALSE]
    volc$y_val <- -log10(volc$pvalue)
    volc$at_cap <- FALSE
    y_lab <- "-log10(pvalue)"
    p_at_fdr <- suppressWarnings(max(volc$pvalue[!is.na(volc$padj) & volc$padj < alpha],
                                     na.rm = TRUE))
    y_fdr <- if (is.finite(p_at_fdr)) -log10(p_at_fdr) else NA_real_
  }
  v_sig <- !is.na(volc$padj) & volc$padj < alpha & lfc_pass(volc$lfc_used)
  volc$regulation <- factor(ifelse(!v_sig, "NS", ifelse(volc$lfc_used > 0, "Up", "Down")),
                            levels=c("Up","Down","NS"))
  vu <- sum(volc$regulation=="Up"); vd <- sum(volc$regulation=="Down")
  vcols <- c(Up="firebrick", Down="royalblue3", NS="grey75")
  vline <- if (alt_hyp %in% c("greaterAbs","lessAbs")) c(-lfc_rep, lfc_rep)
           else if (alt_hyp == "greater") lfc_rep else -lfc_rep
  lab_col <- if ("gene_name" %in% colnames(volc)) "gene_name" else "gene_id"

  if (has("ggplot2")) {
    # Same tie problem as the table ordering: padj = 0 makes the top genes
    # indistinguishable, so label by the exact value when we have it.
    ord <- volc[order(if (exact_ok) -volc$neg_log10_padj else volc$padj), ]
    toplab <- rbind(head(ord[ord$regulation=="Up", , drop=FALSE], 15),
                    head(ord[ord$regulation=="Down", , drop=FALSE], 15))
    gg <- ggplot2::ggplot(volc, ggplot2::aes(lfc_used, y_val, color=regulation)) +
      ggplot2::geom_point(alpha=.6, size=1) +
      ggplot2::scale_color_manual(values=vcols, drop=FALSE,
        name=paste0("padj<", alpha, " & ", alt_hyp, " ", lfc_rep),
        labels=c(paste0("Up (",vu,")"), paste0("Down (",vd,")"), "Not sig.")) +
      ggplot2::geom_vline(xintercept=vline, lty=2, col="grey40") +
      # Open frame: x and y axis lines only. A full box plus a grid competes with
      # the point cloud, which is the whole content of a volcano.
      ggplot2::theme_classic() +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey20"),
                     panel.grid = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank()) +
      ggplot2::labs(title=paste0("Volcano: ", contrast_nm),
                    subtitle=paste0("Up ", vu, " | Down ", vd, "   (padj<", alpha, ", ", alt_hyp,
                                    " ", lfc_rep, " = ", sprintf("%.2f", 2^lfc_rep),
                                    "-fold; x = ", eff_col,
                                    if (is.na(y_fdr)) paste0("; no gene reached padj<", alpha)
                                    else if (volcano_y == "padj") paste0("; dotted = padj ", alpha)
                                    else paste0("; dotted = raw p at padj<", alpha), ")"),
                    x=eff_col, y=y_lab)
    if (!is.na(y_fdr))
      gg <- gg + ggplot2::geom_hline(yintercept=y_fdr, lty=3, col="grey60")
    # Only reached when recomputation was not possible; a faint rule reads as a
    # limit without competing with the points.
    if (any(volc$at_cap))
      gg <- gg + ggplot2::geom_hline(yintercept = y_cap, lty = 2, col = "grey55",
                                     linewidth = 0.3)
    if (has("ggrepel") && nrow(toplab) > 0)
      gg <- gg + ggrepel::geom_text_repel(data=toplab, ggplot2::aes(label=.data[[lab_col]]),
                                          size=3, max.overlaps=20, show.legend=FALSE)
    io$gg("33_volcano", gg, width = 9)
  } else io$base("33_volcano", function() {
    plot(volc$lfc_used, volc$y_val, pch=19, cex=.4,
         col=vcols[as.character(volc$regulation)], main=paste0("Volcano: ", contrast_nm),
         xlab=eff_col, ylab=y_lab)
    abline(v=vline, lty=2)
    if (!is.na(y_fdr)) abline(h=y_fdr, lty=3)
    if (any(volc$at_cap)) abline(h=y_cap, lty=2, col="grey55") })

  top_genes <- head(sig$gene_id, min(top_n, 12))
  if (!length(top_genes)) top_genes <- head(res_df$gene_id[!is.na(res_df$padj)], min(top_n, 12))
  if (length(top_genes) && has("ggplot2")) {
    dd <- do.call(rbind, lapply(top_genes, function(gid) {
      d <- plotCounts(dds, gene=gid, intgroup="condition", returnData=TRUE)
      nm <- if (!is.null(ann)) ann$gene_name[match(gid, ann$gene_id)] else NA
      d$gene <- if (length(nm) && !is.na(nm)) as.character(nm)[1] else gid; d }))
    dd$gene <- factor(as.character(dd$gene), levels=unique(as.character(dd$gene)))
    io$gg("34_top_gene_counts",
      ggplot2::ggplot(dd, ggplot2::aes(condition, count, color=condition)) +
        ggplot2::geom_point(position=ggplot2::position_jitter(width=.15)) +
        ggplot2::scale_y_log10() + ggplot2::facet_wrap(~ gene, scales="free_y") +
        ggplot2::theme_bw() + ggplot2::ggtitle("Normalized counts: top features") +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)),
      width=10, height=8)
  }

  sig_ids <- head(sig$gene_id, top_n)
  sig_ids <- sig_ids[sig_ids %in% rownames(assay(qc_obj))]
  if (length(sig_ids) >= 2 && has("pheatmap")) {
    md <- assay(qc_obj)[sig_ids, , drop=FALSE]
    md <- md[matrixStats::rowSds(md) > 0, , drop=FALSE]   # scale() -> NaN otherwise
    if (nrow(md) >= 2) {
      md <- t(scale(t(md)))
      rn2 <- if (!is.null(ann)) { m2 <- ann$gene_name[match(rownames(md), ann$gene_id)]
                                  ifelse(is.na(m2), rownames(md), m2) } else rownames(md)
      rownames(md) <- make.unique(rn2)
      anndf <- data.frame(condition=coldata$condition); rownames(anndf) <- colnames(qc_obj)
      io$base("35_top_DEG_heatmap", function() pheatmap::pheatmap(md, annotation_col=anndf,
        show_rownames=(nrow(md) <= 60),
        main=paste0("Top ", nrow(md), " DEGs (", qc_lab, " z-score)")), width=8, height=10)
    }
  }

  ## --- manifest ------------------------------------------------------------
  gi <- git_info(.self)
  writeLines(c(
    paste0("rnaseq DESeq2 pipeline ", PIPELINE_VERSION, " -- run manifest"),
    paste0("timestamp: ", format(Sys.time())),
    paste0("host: ", Sys.info()[["nodename"]]),
    paste0("git_commit: ", gi$commit,
           if (isTRUE(gi$dirty)) "  [DIRTY - uncommitted changes]" else
           if (identical(gi$dirty, FALSE)) "  [clean]" else ""),
    paste0("R: ", R.version.string),
    paste0("DESeq2: ", as.character(packageVersion("DESeq2"))),
    paste0("level: ", L$level_lab),
    paste0("ref: ", st1$ref, "   gtf: ", st2$gtf_path),
    paste0("contrast: ", treat, " vs ", control, " (", contrast_nm, ")"),
    paste0("design: ", paste(deparse(design_formula), collapse="")),
    paste0("samples_in_model: ", if (is.null(st1$model_mode)) "contrast-only" else st1$model_mode),
    paste0("test: ", toupper(test_type),
           if (test_type=="lrt") paste0("  reduced: ", paste(deparse(st3$reduced_formula), collapse="")) else ""),
    paste0("prefilter: >=", st3$min_count, " counts in >=", L$min_samples_used, " samples"),
    paste0("alpha(FDR): ", alpha),
    paste0("lfc_report(DEG call): ", lfc_rep, " = ", sprintf("%.3f", 2^lfc_rep), "-fold"),
    paste0("lfc_threshold(test null): ", lfc_thr, "   alt_hypothesis: ", alt_hyp),
    paste0("effect_size_used: ", eff_col, "   shrink: ", shrink_used,
           "   s-values: ", "svalue" %in% colnames(res_df)),
    paste0("independent_filtering: ", indep_f),
    paste0("volcano_y_axis: ", volcano_y),
    paste0("n_samples: ", ncol(dds), "   n_features_tested: ", sum(!is.na(res$pvalue))),
    paste0("n_DEG: ", nrow(sig), " (up=", n_up, ", down=", n_dn, ")"),
    paste0("strandedness: ", st2$strand_verdict),
    paste0("library_type: ", st2$libtype_verdict),
    paste0("args: ", paste(.args_raw, collapse=" ")),
    "outputs:", paste0("  tables:  ", io$tab_dir), paste0("  figures: ", io$fig_dir)),
    file.path(L$level_dir, "MANIFEST.txt"))
  writeLines(capture.output(sessionInfo()), file.path(L$level_dir, "sessionInfo.txt"))
  logmsg("level ", L$level_lab, " done: ", nrow(sig), " DEGs -> ", L$level_dir)
  invisible(nrow(sig))
}

for (nm in names(st3$levels)) results_level(st3$levels[[nm]])

logmsg("STAGE 4 COMPLETE")
cat("\n==> Finished. See ", outdir, "\n", sep = "")
log_close()
