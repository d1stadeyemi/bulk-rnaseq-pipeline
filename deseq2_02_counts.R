#!/usr/bin/env Rscript
###############################################################################
# deseq2_02_counts.R   --  STAGE 2 of 4:  COUNTS, ANNOTATION, LIBRARY QC
#
# Builds the integer count matrix, attaches gene symbols/biotypes, and asks one
# question that decides whether the whole analysis is interpretable:
#
#     ARE THESE LIBRARIES COMPARABLE TO EACH OTHER?
#
# Still no statistics. Everything here is descriptive, and every diagnostic is
# written to disk whether or not it looks bad.
#
# WRITES
#   tables/02_count_assignment_qc      per-sample read fate from featureCounts
#   tables/02b_strandedness_qc         read back from *_library_qc.tsv
#   tables/02c_library_type_qc         poly-A vs total-RNA call per sample
#   tables/05_annotation_report        how much of the matrix got annotated
#   tables/05_biotype_composition_per_sample
#   figures/02_count_assignment_bars, 03_library_size, 03_genes_detected,
#           04_count_distribution_raw, 05_biotype_composition_heatmap
#   _state/02_state.rds
#
#  ==========================================================================
#  PARAMETERS OWNED BY THIS STAGE
#  ==========================================================================
#  --gtf <FILE>            Override the annotation GTF chosen by --ref.
#  --id-map <FILE>         CSV/TSV mapping gene_id -> gene_name[,biotype].
#                          Overrides GTF symbols. This is how you get symbols
#                          for hg19 UCSC knownGene: supply a kgXref table with
#                          columns gene_id,gene_name (kgID -> geneSymbol).
#  --no-annotation         Skip symbols entirely; results show gene IDs.
#  --strip-gene-version    Strip ".NN" from Ensembl IDs (counts AND GTF).
#                          Refuses if stripping would create duplicate IDs.
#  --biotype-filter <STR>  Keep only this biotype before testing, e.g.
#                          protein_coding. Requires annotation.
#  --exclude-chrom <LIST>  Drop all genes on these chromosomes, comma-separated
#                          and matched against the annotation's `chrom` column,
#                          e.g. "chrM" or "chrM,chrY". Composes with
#                          --biotype-filter: biotype is applied first, then this.
#                          Names must match the annotation exactly -- "chrM" for
#                          GENCODE, "MT" for Ensembl -- and a name matching
#                          nothing is an ERROR, so a typo cannot silently
#                          do nothing.
#                          The mitochondrial genes are protein_coding, so
#                          --biotype-filter protein_coding does NOT remove them.
#                          Use this when a handful of very highly expressed
#                          chrM transcripts dominate a principal component and
#                          you want a sensitivity analysis without them. Report
#                          both versions: removing genes because they dominate a
#                          PC is defensible only if the exclusion is stated.
#  --collapse-to-gene      Also run a gene-level analysis in stage 3/4 by
#                          summing transcript-model counts per symbol.
#
#  ==========================================================================
#  ORDER MATTERS: QC RUNS **BEFORE** --biotype-filter
#  ==========================================================================
#  The non-polyA marker set is entirely histones, sno/scaRNA, RN7SL, RMRP and
#  RPPH1 -- none of them protein_coding. If the biotype filter ran first,
#  `--biotype-filter protein_coding` would delete every marker and the
#  library-type diagnostic would silently report nothing. That is precisely the
#  situation in which you need it, because filtering to protein_coding is the
#  recommended WORKAROUND for a poly-A vs total-RNA confound. So all
#  composition QC is computed on the unfiltered matrix, and the filter is
#  applied afterwards.
#
#  ==========================================================================
#  THE TWO CONFOUNDS THIS STAGE IS DESIGNED TO CATCH
#  ==========================================================================
#  1. LIBRARY TYPE (poly-A vs ribo-depleted). Poly-A selection removes
#     replication-dependent histones, Pol III transcripts and sno/scaRNAs;
#     total-RNA libraries keep them. If the CALL differs BY CONDITION, then
#     prep is confounded with biology and those gene classes will look
#     differential for purely technical reasons. Reported as a hard warning.
#  2. STRANDEDNESS. Read back from the upstream *_library_qc.tsv, including
#     `stranded_used` -- what featureCounts -s actually produced the counts.
#     If samples were COUNTED differently from each other, the matrix is not
#     internally comparable and nothing downstream can fix it.
#
#  AMBIGUITY THRESHOLD  (--ambiguous-warn-pct, default 25)
#  -------------------
#  Flagged when ambiguous fragments exceed the threshold, measured against the
#  COUNTABLE denominator (total minus alignment RECORDS and unmapped), not the
#  raw featureCounts total. 12-18% is normal for unstranded (-s 0) counting on
#  GENCODE, where strand cannot separate a gene from its antisense neighbour.
#  The value used is recorded in the qc_thresholds column. On hg19 UCSC
#  knownGene, `-g gene_id` groups by TRANSCRIPT accession, so multi-isoform
#  genes fragment and reads over shared exons are discarded as ambiguous.
#  --collapse-to-gene recovers the isoform-split signal but cannot recover
#  reads featureCounts already threw away -- only re-counting with a gene-level
#  annotation can. grch38/GENCODE does not have this problem.
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
log_open(outdir, "2/counts")
logbanner("STAGE 2: COUNTS, ANNOTATION, LIBRARY QC")

st1 <- load_state(outdir, 1, opt, force = g$bool("force"))
coldata       <- st1$coldata
ref           <- st1$ref
counts_suffix <- st1$counts_suffix
io <- make_io(outdir)

# ---- stage-2 parameters ----------------------------------------------------
default_gtf <- c(
  hg19   = "/lustre/isaac24/proj/UTK0009/RNA-Seq/reference_genomes/hg19/HomoSapiens_hg19.gtf",
  grch38 = "/lustre/isaac24/proj/UTK0009/RNA-Seq/reference_genomes/grch38/gencode.v49.annotation.gtf")
gtf_path      <- g$chr("gtf", default_gtf[[ref]])
id_map        <- g$chr("id-map")
do_annotate   <- !g$bool("no-annotation")
strip_version <- g$bool("strip-gene-version")
biotype_filt  <- g$chr("biotype-filter")
exclude_chrom <- g$chr("exclude-chrom")
collapse_gene <- g$bool("collapse-to-gene")

rlib <- g$chr("rlib")
if (!is.null(rlib)) { if (!dir.exists(rlib)) die("--rlib not found: ", rlib); .libPaths(c(rlib, .libPaths())) }

# Settable AND recorded, matching QC_MAX_AMBIGUOUS in rnaseq_fastq_to_counts.sh
# so both files call the same thing a problem. A hardcoded gate makes an old QC
# table uninterpretable: "OK" means nothing unless you know what it was judged
# against, so the value travels in the qc_thresholds column.
# 25 rather than 15 because the denominator is COUNTABLE, not the raw
# featureCounts total, which roughly doubles the reported value (6.4% -> 11.9%
# on this data). 12-18% is normal for unstranded counting on GENCODE; 15 would
# fire on every sample. It still catches what this check exists for: hg19 UCSC
# knownGene, where -g gene_id groups by transcript accession and ambiguity
# runs 30-50%.
AMBIGUOUS_WARN_PCT <- g$num("ambiguous-warn-pct", 25, min = 0, max = 100)

# Library-type call boundaries. Named, not literals buried in an ifelse.
# These define what "poly-A" MEANS rather than setting a pass/fail tolerance,
# so they are constants rather than flags -- but they are recorded below so the
# call in 02c_library_type_qc can be re-derived years later.
POLYA_MAX_PCT    <- 2    # below this -> poly-A selected
TOTALRNA_MIN_PCT <- 5    # above this -> total RNA / ribo-depleted

## ===========================================================================
## 2.1  READ COUNT FILES
## ===========================================================================
# Accepts htseq-style 2-column files AND native featureCounts output, so
# pointing --counts-suffix at *_featurecounts.txt gives a result instead of a
# parse error. STAR ReadsPerGene.out.tab is refused with instructions.
read_counts <- function(path) {
  first <- readLines(path, n = 2L, warn = FALSE)
  if (length(first) == 0) stop("count file is empty: ", path)
  if (any(grepl("^N_unmapped", first)))
    stop("this is a STAR ReadsPerGene.out.tab file: ", path,
         "\n  It has one count column per strandedness. Extract the right column first.")
  is_fc <- startsWith(first[1], "# Program:featureCounts") ||
           (length(first) >= 2 && startsWith(first[2], "Geneid\t"))
  if (is_fc) {
    df0 <- read.table(path, sep = "\t", header = TRUE, skip = 1, stringsAsFactors = FALSE,
                      quote = "", comment.char = "", check.names = FALSE)
    if (ncol(df0) < 7) stop("featureCounts file has ", ncol(df0), " columns (expected >=7): ", path)
    if (ncol(df0) > 7) stop("featureCounts file has ", ncol(df0)-6, " count columns (expected 1): ", path)
    df <- data.frame(gene_id = as.character(df0[[1]]), count = as.numeric(df0[[7]]),
                     stringsAsFactors = FALSE)
    sp <- setNames(numeric(0), character(0))
    sfile <- paste0(path, ".summary")
    if (file.exists(sfile)) {
      s <- read.table(sfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "", check.names = FALSE)
      key <- c(Unassigned_NoFeatures = "__no_feature", Unassigned_Ambiguity = "__ambiguous",
               Unassigned_MappingQuality = "__too_low_aQual", Unassigned_Unmapped = "__not_aligned",
               Unassigned_MultiMapping = "__alignment_not_unique")
      idx <- match(names(key), s[[1]]); okk <- !is.na(idx)
      sp <- setNames(as.numeric(s[[2]][idx[okk]]), unname(key[okk]))
    }
    return(list(genes = df, special = sp))
  }
  df <- read.table(path, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("gene_id","count"), quote = "", comment.char = "")
  df$gene_id <- as.character(df$gene_id); df$count <- as.numeric(df$count)
  sp <- startsWith(df$gene_id, "__")
  list(genes = df[!sp, , drop = FALSE], special = setNames(df$count[sp], df$gene_id[sp]))
}

logmsg("reading ", nrow(coldata), " count files ...")
reads <- lapply(coldata$count_file, function(p)
  tryCatch(read_counts(p), error = function(e) { logmsg("while reading ", p); die(conditionMessage(e)) }))

## ---- assignment QC ---------------------------------------------------------
# CAVEAT recorded in the table: STAR writes only ALIGNED reads to its BAM, so
# __not_aligned is ~0 and `total` here means "fragments that reached
# featureCounts", NOT fragments sequenced. For the true mapping rate use
# <sample>_Log.final.out (summarised in the upstream *_library_qc.tsv).
special_names <- c("__no_feature","__ambiguous","__too_low_aQual",
                   "__not_aligned","__alignment_not_unique")
qc <- data.frame(sample = coldata$sample, condition = as.character(coldata$condition),
                 stringsAsFactors = FALSE)
qc$assigned <- vapply(reads, function(r) sum(r$genes$count), numeric(1))
for (sn in special_names)
  qc[[sub("^__","",sn)]] <- vapply(reads, function(r) {
    v <- r$special[[sn]]; if (is.null(v) || is.na(v)) 0 else as.numeric(v) }, numeric(1))
qc$total <- rowSums(qc[, c("assigned", sub("^__","",special_names))])

# COUNTABLE denominator, defined by SUBTRACTION rather than by summing named
# categories. featureCounts reports Unassigned_MultiMapping in ALIGNMENT
# RECORDS (STAR writes up to --outFilterMultimapNmax per read), so including it
# divides fragments by records: the ratio then tracks the library's multimapping
# rate, not its assignment rate, and is not comparable between samples.
# Unmapped is excluded because a fragment that never aligned could not have been
# assigned to a gene. Everything else -- NoFeatures, Ambiguity, MappingQuality,
# and any category a future --fc-extra switches on -- stays IN automatically.
qc$countable <- qc$total - qc$alignment_not_unique - qc$not_aligned
# Distinct names on purpose. Reusing `assignment_rate` for the new denominator
# would silently change the meaning of a column that already exists in older
# output. These names also match rnaseq_fastq_to_counts.sh, so the same quantity
# is called the same thing everywhere in the pipeline.
qc$assignment_rate_all_records <- round(qc$assigned / pmax(qc$total, 1), 4)
qc$assignment_rate_countable   <- round(qc$assigned / pmax(qc$countable, 1), 4)
qc$pct_ambiguous   <- round(100 * qc$ambiguous / pmax(qc$countable, 1), 2)
qc$pct_multimap_records <- round(100 * qc$alignment_not_unique / pmax(qc$total, 1), 2)
qc$denominator_note <- paste0(
  "assignment_rate_countable/pct_ambiguous use COUNTABLE = total - alignment_not_unique - ",
  "not_aligned (fragments on both sides). assignment_rate_all_records uses the raw ",
  "featureCounts total, which mixes fragments with alignment RECORDS -- do not QC on it. ",
  "pct_multimap_records scales with --outFilterMultimapNmax; diagnostic only.")
qc$qc_thresholds <- paste0("ambiguous_warn_pct=", AMBIGUOUS_WARN_PCT)
io$table(qc, "02_count_assignment_qc.csv")
logmsg("assignment (countable denominator): ",
       paste(sprintf("%s=%.1f%%", qc$sample, 100 * qc$assignment_rate_countable), collapse = "  "))

amb_bad <- qc$pct_ambiguous > AMBIGUOUS_WARN_PCT
if (any(amb_bad)) {
  logmsg("WARNING: ", sum(amb_bad), " sample(s) have >", AMBIGUOUS_WARN_PCT,
         "% ambiguous fragments (max ", max(qc$pct_ambiguous),
         "%, of the COUNTABLE denominator): ", paste(qc$sample[amb_bad], collapse = ", "))
  logmsg("         On hg19 UCSC knownGene this is expected: `-g gene_id` groups by ",
         "TRANSCRIPT accession, so multi-isoform genes fragment and shared-exon reads ",
         "are discarded. --collapse-to-gene recovers the isoform SPLIT but NOT the ",
         "discarded reads -- only re-counting against a gene-level annotation does. ",
         "On grch38/GENCODE, investigate instead (wrong -s, or overlapping features).")
} else logmsg("ambiguity: max ", max(qc$pct_ambiguous), "% (threshold ", AMBIGUOUS_WARN_PCT, "%) -- OK")

## ---- build the matrix ------------------------------------------------------
gene_lists   <- lapply(reads, function(r) r$genes$gene_id)
common_genes <- Reduce(intersect, gene_lists)
if (length(common_genes) == 0) die("no gene IDs shared across samples.")
if (any(vapply(gene_lists, length, integer(1)) != length(common_genes)))
  logmsg("WARNING: samples had differing gene sets; using the ", length(common_genes),
         " IDs common to all.")

cts <- vapply(reads, function(r) r$genes$count[match(common_genes, r$genes$gene_id)],
              numeric(length(common_genes)))
if (!is.matrix(cts)) cts <- matrix(cts, nrow = length(common_genes))
rownames(cts) <- common_genes; colnames(cts) <- coldata$sample

if (strip_version) {
  newrn <- sub("\\.\\d+$", "", rownames(cts))
  if (anyDuplicated(newrn)) {
    dup <- unique(newrn[duplicated(newrn)])
    die("--strip-gene-version created ", length(dup), " duplicate ID(s), e.g. ",
        paste(head(dup, 5), collapse = ", "), ". Re-run without it.")
  }
  rownames(cts) <- newrn
}
if (any(is.na(cts))) die("count matrix contains NA after aligning on common gene IDs.")
if (any(cts < 0))    die("count matrix contains negative values -- these are not raw counts.")
if (any((cts - floor(cts)) > 1e-8)) {
  logmsg("WARNING: non-integer counts (Salmon/RSEM estimates?). DESeq2 needs raw integer ",
         "counts; rounding. For transcript estimates use tximport instead.")
  cts <- round(cts)
}
storage.mode(cts) <- "integer"
logmsg("count matrix: ", nrow(cts), " features x ", ncol(cts), " samples")

cond_int <- as.integer(coldata$condition)
io$base("02_count_assignment_bars", function() {
  mat <- t(as.matrix(qc[, c("assigned", sub("^__","",special_names))])); colnames(mat) <- qc$sample
  cols <- c(assigned="#3B8C4E", no_feature="#F2A65A", ambiguous="#F2CC8F",
            too_low_aQual="#BDBDBD", not_aligned="#7F9DB9", alignment_not_unique="#C1483B")[rownames(mat)]
  par(mar = c(9,5,3,9), xpd = TRUE)
  barplot(mat, col = cols, border = NA, las = 2, ylab = "fragments", main = "Fragment assignment per sample")
  legend("topright", inset = c(-0.30,0), legend = rownames(mat), fill = cols, border = NA, bty = "n", cex = .75)
}, width = 9, height = 6)
io$base("03_library_size", function() { par(mar=c(9,5,3,2))
  barplot(colSums(cts), las=2, ylab="assigned counts", main="Library size", col=cond_int) }, width=8, height=5)
detected <- colSums(cts > 0)
io$table(data.frame(sample = names(detected), features_detected = as.integer(detected)), "03_genes_detected.csv")
io$base("03_genes_detected", function() { par(mar=c(9,5,3,2))
  barplot(detected, las=2, ylab="# features with count > 0", main="Detected features", col=cond_int) }, width=8, height=5)
io$base("04_count_distribution_raw", function() { par(mar=c(9,5,3,2))
  boxplot(log2(cts+1), las=2, ylab="log2(count+1)", main="Raw count distribution", col=cond_int) }, width=8, height=5)

## ===========================================================================
## 2.2  ANNOTATION
## ===========================================================================
annotate_genes <- function(gtf, strip_version) {
  if (!file.exists(gtf)) { logmsg("WARNING: GTF not found (", gtf, "); skipping annotation."); return(NULL) }
  ann <- NULL
  if (has("rtracklayer")) {
    ann <- tryCatch({
      gr <- rtracklayer::import(gtf); md <- as.data.frame(S4Vectors::mcols(gr))
      if ("type" %in% colnames(md) && any(md$type == "gene")) {
        keep <- md$type == "gene"; gr <- gr[keep]; md <- md[keep, , drop = FALSE]
      }
      gname <- if ("gene_name" %in% colnames(md)) md$gene_name else md$gene_id
      btype <- if ("gene_type" %in% colnames(md)) md$gene_type
               else if ("gene_biotype" %in% colnames(md)) md$gene_biotype else NA
      data.frame(gene_id = as.character(md$gene_id),
                 gene_name = ifelse(is.na(gname), as.character(md$gene_id), as.character(gname)),
                 biotype = btype, chrom = as.character(GenomeInfoDb::seqnames(gr)),
                 stringsAsFactors = FALSE)
    }, error = function(e) { logmsg("rtracklayer failed (", conditionMessage(e), "); manual parse."); NULL })
  }
  if (is.null(ann)) {
    grab <- function(line, key) {
      m <- regmatches(line, regexpr(paste0(key, ' "[^"]+"'), line))
      if (length(m)) sub(paste0(key, ' "([^"]+)"'), "\\1", m) else NA_character_
    }
    stream_parse <- function(feature) {
      con <- file(gtf, "r"); on.exit(close(con))
      ids <- nms <- bts <- chr <- character(0)
      repeat {
        chunk <- readLines(con, n = 100000L); if (length(chunk) == 0) break
        chunk <- chunk[!startsWith(chunk, "#")]; if (length(chunk) == 0) next
        parts <- strsplit(chunk, "\t", fixed = TRUE)
        f1 <- vapply(parts, function(x) if (length(x)>=1) x[1] else "", character(1))
        f3 <- vapply(parts, function(x) if (length(x)>=3) x[3] else "", character(1))
        if (!is.null(feature)) { sel <- f3 == feature; chunk <- chunk[sel]; f1 <- f1[sel] }
        if (length(chunk) == 0) next
        ids <- c(ids, vapply(chunk, grab, character(1), key="gene_id"))
        nms <- c(nms, vapply(chunk, grab, character(1), key="gene_name"))
        bts <- c(bts, vapply(chunk, function(l){b<-grab(l,"gene_type"); if(is.na(b)) grab(l,"gene_biotype") else b}, character(1)))
        chr <- c(chr, f1)
      }
      data.frame(gene_id=ids, gene_name=nms, biotype=bts, chrom=chr, stringsAsFactors=FALSE)
    }
    ann <- tryCatch({
      gg <- stream_parse("gene"); if (nrow(gg) == 0) gg <- stream_parse(NULL)
      gg$gene_name <- ifelse(is.na(gg$gene_name), gg$gene_id, gg$gene_name); gg
    }, error = function(e) { logmsg("manual GTF parse failed: ", conditionMessage(e)); NULL })
  }
  if (is.null(ann) || nrow(ann) == 0) return(NULL)
  if (strip_version) ann$gene_id <- sub("\\.\\d+$", "", ann$gene_id)
  ann[!duplicated(ann$gene_id), , drop = FALSE]
}

map_ucsc_to_symbol <- function(ann) {
  if (is.null(ann)) return(ann)
  if (mean(grepl("^uc[0-9]", ann$gene_id)) < 0.5) return(ann)
  need <- c("TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","AnnotationDbi")
  if (!all(vapply(need, has, logical(1)))) {
    logmsg("hg19 knownGene IDs detected but ", paste(need[!vapply(need,has,logical(1))], collapse=", "),
           " not installed. Supply --id-map (UCSC kgXref) for symbols."); return(ann)
  }
  tryCatch({
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    ent <- AnnotationDbi::mapIds(txdb, keys=ann$gene_id, column="GENEID", keytype="TXNAME", multiVals="first")
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=as.character(ent),
                                 column="SYMBOL", keytype="ENTREZID", multiVals="first")
    ann$gene_name <- ifelse(is.na(sym), ann$gene_name, sym)
    logmsg("mapped UCSC knownGene -> symbols (", sum(!is.na(sym)), "/", nrow(ann), ")"); ann
  }, error = function(e) { logmsg("UCSC->symbol mapping failed: ", conditionMessage(e)); ann })
}

ann <- if (do_annotate) annotate_genes(gtf_path, strip_version) else NULL
if (do_annotate && !is.null(id_map)) {
  if (!file.exists(id_map)) die("--id-map not found: ", id_map)
  m <- read_delim_auto(id_map)
  if (!all(c("gene_id","gene_name") %in% colnames(m)))
    die("--id-map must have 'gene_id' and 'gene_name' columns.")
  m$gene_id <- as.character(m$gene_id); m$gene_name <- as.character(m$gene_name)
  if (strip_version) m$gene_id <- sub("\\.\\d+$", "", m$gene_id)
  m <- m[!duplicated(m$gene_id), , drop = FALSE]
  if (is.null(ann)) {
    ann <- data.frame(gene_id=m$gene_id, gene_name=m$gene_name,
                      biotype=if ("biotype" %in% colnames(m)) m$biotype else NA_character_,
                      chrom=NA_character_, stringsAsFactors=FALSE)
  } else {
    idx <- match(ann$gene_id, m$gene_id); hit <- !is.na(idx)
    ann$gene_name[hit] <- m$gene_name[idx[hit]]
    if ("biotype" %in% colnames(m)) ann$biotype[hit] <- m$biotype[idx[hit]]
  }
  logmsg("applied --id-map: ", id_map, " (", nrow(m), " entries)")
} else if (do_annotate) ann <- map_ucsc_to_symbol(ann)

if (!is.null(ann)) {
  atomize <- function(x) if (is.list(x))
    vapply(x, function(v) if (length(v)) as.character(v)[[1]] else NA_character_, character(1)) else as.character(x)
  for (.c in intersect(c("gene_id","gene_name","biotype","chrom"), names(ann))) ann[[.c]] <- atomize(ann[[.c]])
  ann$gene_name <- ifelse(is.na(ann$gene_name) | !nzchar(ann$gene_name), ann$gene_id, ann$gene_name)
  ann <- ann[!duplicated(ann$gene_id), , drop = FALSE]

  n_match <- sum(rownames(cts) %in% ann$gene_id)
  matched <- ann[ann$gene_id %in% rownames(cts), , drop = FALSE]
  sym_frac <- if (nrow(matched) > 0) mean(!is.na(matched$gene_name) & matched$gene_name != matched$gene_id) else 0
  io$table(data.frame(gtf = gtf_path, id_map = if (is.null(id_map)) "<none>" else id_map,
                      n_features = nrow(cts), n_annotated = n_match,
                      frac_annotated = round(n_match/nrow(cts), 3),
                      symbols_available = sym_frac > 0.1,
                      frac_with_symbol = round(sym_frac, 3), stringsAsFactors = FALSE),
           "05_annotation_report.csv")
  logmsg("annotation: ", n_match, "/", nrow(cts), " IDs matched; symbols on ", round(100*sym_frac), "%")
  if (n_match/nrow(cts) < 0.5)
    logmsg("WARNING: <50% of IDs matched the annotation -- check --ref / --strip-gene-version / --gtf.")
} else if (do_annotate) logmsg("annotation unavailable; results will use gene IDs only.")

## ===========================================================================
## 2.3  LIBRARY QC  (BEFORE the biotype filter -- see header)
## ===========================================================================
strand_verdict <- "not available (no *_library_qc.tsv beside the count files)"
lib_qc_files <- vapply(coldata$count_file, function(p) {
  s <- sub(paste0(counts_suffix, "$"), "", basename(p))
  f <- file.path(dirname(p), paste0(s, "_library_qc.tsv")); if (file.exists(f)) f else NA_character_
}, character(1))
if (any(!is.na(lib_qc_files))) {
  rows <- lapply(lib_qc_files[!is.na(lib_qc_files)], function(f)
    tryCatch(read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  # rbind() ERRORS on data.frames with different columns, and a bare tryCatch
  # would swallow that to NULL -- silently disabling the strandedness check in
  # exactly the mixed-schema case it exists to catch. Bind on shared columns and
  # say so instead.
  sq <- NULL
  if (length(rows)) {
    common <- Reduce(intersect, lapply(rows, colnames))
    allcol <- unique(unlist(lapply(rows, colnames)))
    if (length(common) < length(allcol)) {
      logmsg("WARNING: the *_library_qc.tsv files do not share a column layout ",
             "(", length(common), " of ", length(allcol), " columns common). Some samples ",
             "were produced by a DIFFERENT version of rnaseq_fastq_to_counts.sh. ",
             "Dropped column(s): ", paste(setdiff(allcol, common), collapse = ", "),
             ".  Re-emit the QC rows so every sample is described the same way.")
    }
    if (length(common) > 0)
      sq <- tryCatch(do.call(rbind, lapply(rows, function(d) d[, common, drop = FALSE])),
                     error = function(e) { logmsg("WARNING: could not combine ",
                       "*_library_qc.tsv files: ", conditionMessage(e)); NULL })
  }
  if (!is.null(sq) && "strand_verdict" %in% colnames(sq)) {
    io$table(sq, "02b_strandedness_qc.csv")
    v <- unique(sq$strand_verdict)
    if (length(v) == 1) { strand_verdict <- v; logmsg("strandedness (all samples): ", v) }
    else {
      strand_verdict <- paste0("MIXED: ", paste(v, collapse=", "))
      logmsg("WARNING: samples disagree on detected strandedness (", paste(v, collapse=", "),
             "). Different preps => a batch effect possibly confounded with condition.")
    }
    if ("stranded_used" %in% colnames(sq)) {
      su <- unique(sq$stranded_used)
      logmsg("counted upstream with -s: ", paste(su, collapse=", "))
      if (length(su) > 1)
        logmsg("WARNING: samples were COUNTED with DIFFERENT -s settings. The count matrix ",
               "is not internally comparable and no downstream step can repair that. ",
               "Re-count everything the same way.")
    } else logmsg("NOTE: *_library_qc.tsv has no 'stranded_used' column, so this run cannot ",
                  "record which featureCounts -s produced the counts.")
  }
}

nonpolya_pattern <- paste(c("^H1-[0-9]+$","^H2AC[0-9]+","^H2BC[0-9]+","^H3C[0-9]+","^H4C[0-9]+",
  "^HIST[0-9]","^RMRP$","^RPPH1$","^RN7SL[0-9]","^RN7SK$","^RNU[0-9]","^SNORD","^SNORA",
  "^SCARNA","^VTRNA","^RNY[0-9]","^RNA5S","^RNA5-8S"), collapse="|")
libtype_verdict <- "not assessed (no gene symbols available)"
if (!is.null(ann)) {
  sym <- ann$gene_name[match(rownames(cts), ann$gene_id)]
  sym[is.na(sym)] <- rownames(cts)[is.na(sym)]
  is_np <- grepl(nonpolya_pattern, sym)
  if (sum(is_np) > 0) {
    pct_np <- 100 * colSums(cts[is_np, , drop=FALSE]) / pmax(colSums(cts), 1)
    call1 <- ifelse(pct_np < POLYA_MAX_PCT, "poly-A selected (likely)",
             ifelse(pct_np > TOTALRNA_MIN_PCT, "total RNA / ribo-depleted (likely)", "ambiguous"))
    lt <- data.frame(sample=colnames(cts), condition=as.character(coldata$condition),
                     n_marker_genes=sum(is_np), pct_counts_nonpolyA=round(as.numeric(pct_np),3),
                     call=call1,
                     call_thresholds=paste0("polyA<", POLYA_MAX_PCT,
                                            "; totalRNA>", TOTALRNA_MIN_PCT),
                     stringsAsFactors=FALSE)
    io$table(lt, "02c_library_type_qc.csv")
    for (i in seq_len(nrow(lt)))
      logmsg("library type: ", lt$sample[i], " non-polyA = ",
             sprintf("%.2f%%", lt$pct_counts_nonpolyA[i]), " -> ", lt$call[i])
    uc <- unique(lt$call)
    libtype_verdict <- if (length(uc)==1) uc else paste0("MIXED: ", paste(uc, collapse=" | "))
    bycond <- tapply(lt$call, lt$condition, function(x) paste(unique(x), collapse="/"))
    if (length(unique(bycond)) > 1) {
      libtype_verdict <- paste0("CONFOUNDED with condition -> ",
                                paste(sprintf("%s: %s", names(bycond), bycond), collapse=" | "))
      logmsg("WARNING: library TYPE differs BY CONDITION (",
             paste(sprintf("%s=%s", names(bycond), bycond), collapse=", "),
             "). Poly-A vs total-RNA prep is confounded with the biology: histones, ",
             "RMRP/RPPH1, sno/scaRNA and RN7SL will appear differential for purely ",
             "technical reasons. Restrict interpretation to poly-A mRNAs (e.g. ",
             "--biotype-filter protein_coding) and confirm kit names with the facility.")
    }
  }
}

if (!is.null(ann) && any(!is.na(ann$biotype))) {
  bt <- ann$biotype[match(rownames(cts), ann$gene_id)]
  bt[is.na(bt) | !nzchar(bt)] <- "unknown"
  comp <- rowsum(cts, group = bt)
  pct <- sweep(comp, 2, pmax(colSums(comp), 1), "/") * 100
  io$table(data.frame(biotype=rownames(pct), round(as.data.frame(pct),4), check.names=FALSE),
           "05_biotype_composition_per_sample.csv")
  if (has("pheatmap") && nrow(pct) > 1) {
    ordbt <- order(apply(pct,1,max), decreasing=TRUE); keepn <- min(15L, nrow(pct))
    top <- pct[ordbt[seq_len(keepn)], , drop=FALSE]
    if (nrow(pct) > keepn) top <- rbind(top, other = colSums(pct[ordbt[-seq_len(keepn)], , drop=FALSE]))
    nums <- matrix(sprintf("%.2f", top), nrow=nrow(top), dimnames=dimnames(top))
    anncol <- data.frame(condition=coldata$condition); rownames(anncol) <- colnames(cts)
    io$base("05_biotype_composition_heatmap", function()
      pheatmap::pheatmap(sqrt(top), display_numbers=nums, number_color="grey15", fontsize_number=8,
        cluster_rows=FALSE, annotation_col=anncol,
        main="% of assigned counts per biotype (numbers = %, colour = sqrt)"), width=8, height=7)
  }
}

## ===========================================================================
## 2.4  BIOTYPE FILTER  (applied last)
## ===========================================================================
if (!is.null(biotype_filt)) {
  if (is.null(ann)) die("--biotype-filter requires annotation.")
  if (!any(!is.na(ann$biotype)))
    die("--biotype-filter requested but the annotation carries no biotype (e.g. hg19 ",
        "UCSC knownGene). Use grch38, or drop --biotype-filter.")
  keep_bt <- rownames(cts) %in% ann$gene_id[!is.na(ann$biotype) & ann$biotype == biotype_filt]
  if (sum(keep_bt) == 0) die("--biotype-filter '", biotype_filt, "' matched 0 features.")
  logmsg("biotype filter '", biotype_filt, "': keeping ", sum(keep_bt), "/", nrow(cts), " features")
  cts <- cts[keep_bt, , drop = FALSE]
}

## ---------------------------------------------------------------------------
## 2.5  CHROMOSOME EXCLUSION  (after the biotype filter, so the two compose)
## ---------------------------------------------------------------------------
# Mitochondrial genes are protein_coding, so --biotype-filter protein_coding
# leaves all 13 of them in place. They are also among the most highly expressed
# transcripts in the matrix, which is why a handful of them can dominate a
# principal component. Removing them is a legitimate sensitivity analysis and
# an illegitimate way to make a PCA look tidier -- so it is opt-in, logged, and
# written into the run config.
if (!is.null(exclude_chrom)) {
  if (is.null(ann)) die("--exclude-chrom requires annotation (it matches the `chrom` column).")
  if (!"chrom" %in% colnames(ann) || all(is.na(ann$chrom)))
    die("--exclude-chrom requested but the annotation carries no chromosome column.\n",
        "       The GTF parse may have failed; check 05_annotation_report.")
  want <- trimws(strsplit(exclude_chrom, ",", fixed = TRUE)[[1]])
  want <- want[nzchar(want)]
  if (!length(want)) die("--exclude-chrom was given an empty list.")
  # A name that matches nothing is an error, not a no-op: "MT" against a GENCODE
  # annotation that says "chrM" would otherwise silently remove nothing.
  present <- unique(ann$chrom[!is.na(ann$chrom)])
  bad <- setdiff(want, present)
  if (length(bad))
    die("--exclude-chrom name(s) not found in the annotation: ", paste(bad, collapse = ", "),
        "\n       Chromosome names in this annotation look like: ",
        paste(head(sort(present), 6), collapse = ", "), " ...\n",
        "       (GENCODE uses 'chrM'; Ensembl uses 'MT'.)")
  drop_ids <- ann$gene_id[!is.na(ann$chrom) & ann$chrom %in% want]
  keep_ch  <- !(rownames(cts) %in% drop_ids)
  n_drop   <- sum(!keep_ch)
  if (n_drop == 0)
    logmsg("NOTE: --exclude-chrom ", exclude_chrom, " matched no gene in the count matrix ",
           "(they may already have been removed by --biotype-filter or prefiltering).")
  else {
    dropped <- ann$gene_name[match(rownames(cts)[!keep_ch], ann$gene_id)]
    dropped <- dropped[!is.na(dropped)]
    logmsg("chromosome filter '", exclude_chrom, "': removing ", n_drop, "/", nrow(cts),
           " features", if (length(dropped) && length(dropped) <= 20)
             paste0(" (", paste(sort(dropped), collapse = ", "), ")") else "")
    if (sum(keep_ch) < 100)
      die("--exclude-chrom left only ", sum(keep_ch), " features. Check the names given.")
    cts <- cts[keep_ch, , drop = FALSE]
  }
}

save_state(list(cts = cts, ann = ann, qc = qc,
                strand_verdict = strand_verdict, libtype_verdict = libtype_verdict,
                gtf_path = gtf_path, collapse_gene = collapse_gene,
                biotype_filt = biotype_filt, exclude_chrom = exclude_chrom),
           outdir, 2, opt)
logmsg("STAGE 2 COMPLETE  (", nrow(cts), " features x ", ncol(cts), " samples)")
log_close()
