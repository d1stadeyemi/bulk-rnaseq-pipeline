#!/bin/bash
#SBATCH -J deseq2_dge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=00-23:00:00
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --error=/lustre/isaac24/proj/UTK0009/RNA-Seq/logs/errors/%x-%j.err
#SBATCH --output=/lustre/isaac24/proj/UTK0009/RNA-Seq/logs/jobs/%x-%j.out

set -euo pipefail
umask 0002        # group-writable outputs so colleagues can rerun/debug

###############################################################################
# submit_rnaseq_deseq2_dge.sh   (v3.1)  --  runs the 4-stage DESeq2 pipeline
#
#  ==========================================================================
#  QUICK START: MANY COMPARISONS IN ONE GO  (--contrasts FILE.csv)
#  ==========================================================================
#  1. Make a CSV, one comparison per row. In Excel: File > Save As >
#     "CSV UTF-8 (Comma delimited)". Three columns; the header row is optional.
#
#         treatment,control,fit_group
#         HUVEC_2Gy_24hr,HUVEC_untreated,HUVEC
#         HUVEC_600mGy_24hr,HUVEC_untreated,HUVEC
#         IMR90_100mGy_24hr,IMR90_untreated,
#         IMR90_1Gy_5day,IMR90_100mGy_5day,ALL
#
#     col 1  treatment  fold-change NUMERATOR   (a value in the sheet's `condition`)
#     col 2  control    fold-change DENOMINATOR (a value in the sheet's `condition`)
#     col 3  fit_group  which samples enter the model (sent to --fit-all-groups):
#                         a pattern, e.g. HUVEC -> every condition containing it
#                         ALL                   -> every sample in the sheet
#                         left blank            -> only these two groups (default)
#     Lines starting with # and empty lines are ignored.
#
#  2. Run it with plain `bash` on the LOGIN node -- NOT sbatch. It submits the
#     jobs itself:
#       bash submit_rnaseq_deseq2_dge.sh --contrasts my_contrasts.csv \
#            --ref grch38 --counts-dir DIR --samplesheet samplesheet.csv \
#            --outdir results/deseq2_batch1   [any other flags]
#
#  3. What happens:
#     - EVERY row is validated first (a stage-1 dry run, a few seconds each).
#       If ANY row fails, NOTHING is submitted, and the failing rows are listed
#       with their error. Full logs: <outdir>/_contrasts_validation/
#     - If all rows pass, ONE SLURM job per row is submitted, each writing to
#       its own folder:   <outdir>/<treatment>_vs_<control>/
#     - Job IDs are recorded in <outdir>/contrasts_submitted.tsv
#     - Add --dry-run to validate only (nothing is submitted).
#     - Every other flag (--batch, --alpha, --stage, ...) applies to EVERY row.
#       Do NOT also pass --control / --treatment / --fit-all-groups /
#       --contrast-name: the CSV supplies them, one row at a time.
#
#  Without --contrasts the script behaves exactly as before: one comparison.
#
# The analysis is split into four R scripts. This wrapper runs them IN ORDER,
# in ONE job, forwarding your flags unchanged to each. All four APPEND to a
# single <outdir>/run.log, tagged [1/setup] [2/counts] [3/model] [4/results].
#
#   deseq2_01_setup.R    samples, groups, design, batch, full-rank check
#   deseq2_02_counts.R   count matrix, annotation, library-composition QC
#   deseq2_03_model.R    prefilter, DESeq() fit, transforms, sample QC
#   deseq2_04_results.R  contrast, shrinkage, thresholds, DEG tables, figures
#
# Each R script documents its OWN parameters in its header, so to judge whether
# a threshold is right you read one header, not 1200 lines. Nothing about the
# split changes the statistics: stage N reads stage N-1's .rds and refuses to
# run if that file was produced with different arguments.
#
#  ==========================================================================
#  USAGE
#  ==========================================================================
#  # Validate the design without touching count data (ALWAYS DO THIS FIRST):
#    sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 \
#        --counts-dir /lustre/.../rnaseq_output \
#        --samplesheet /lustre/.../rnaseq_output/samplesheet.csv \
#        --control Parental --treatment Bottom10 --batch batch \
#        --outdir results/exploratory/deseq2_B10_vs_Parental --dry-run
#
#  # Full run (defaults: padj<0.05 AND >=1.5-fold):
#    sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 \
#        --counts-dir /lustre/.../rnaseq_output \
#        --samplesheet .../samplesheet.csv \
#        --control Parental --treatment Bottom10 --batch batch \
#        --outdir results/exploratory/deseq2_B10_vs_Parental
#
#  # Keyword scan instead of a sample sheet, batches with no file at all:
#    sbatch submit_rnaseq_deseq2_dge.sh --ref grch38 --counts-dir DIR \
#        --group1 Bottom10 --group2 Parental --control Parental \
#        --batch-map "rep1=A,rep2=B,rep3=C" --outdir OUT
#
#  # Change only the thresholds -- re-runs stage 4 in seconds, not the fit:
#    sbatch submit_rnaseq_deseq2_dge.sh <same flags> --stage 4 \
#        --alpha 0.01 --lfc-report 1
#
#  # Strict threshold TEST instead of post-hoc filtering:
#    ... --lfc-threshold 0.585        (p-values then test |LFC| > 0.585)
#
#  # Equivalence test ("which genes did NOT change?"):
#    ... --alt-hypothesis lessAbs --lfc-threshold 0.585
#
#  # Three conditions:
#    ... --control DMSO --treatment UNC1999          # UNC vs DMSO, all groups fitted
#    ... --control DMSO --treatment UNC1999 --subset-to-contrast   # drop the rest
#    ... --test lrt                                   # does condition matter at all
#
#  # Same library across two lanes (technical, NOT biological, replicates):
#    ... --collapse-technical biological_sample
#
#  ==========================================================================
#  SUBMITTER-ONLY FLAGS
#  ==========================================================================
#  --stage <1|2|3|4|all|N-|N-M>   Which stages to run. Default: all
#                                 4      only stage 4
#                                 3-     stages 3 and 4
#                                 2-3    stages 2 and 3
#  --force                        Reuse upstream state even if it was produced
#                                 with different arguments. Forwarded to R.
#  --contrasts <FILE.csv>         Batch mode: one job per CSV row. See QUICK
#                                 START at the top. Run with bash, not sbatch.
#  Everything else is passed to the R scripts unchanged; run --help for the list.
#
#  Prefer interactive? The wrapper is optional:
#    salloc -A ACF-UTK0009 -p condo-rmccord -q condo -c 4 --mem 32G -t 2:00:00
#    conda activate multiomics-r
#    Rscript deseq2_01_setup.R --ref grch38 --counts-dir ... [flags]
#    Rscript deseq2_02_counts.R  <same flags>     # etc.
###############################################################################

# ---------------------------------------------------------------------------
# Locate the R scripts (they sit next to this wrapper).
# Under `sbatch`, SLURM runs a SPOOLED COPY, so ${BASH_SOURCE[0]} points at the
# spool dir. Ask scontrol for the real submitted path when in a job.
# ---------------------------------------------------------------------------
SELF=""
if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v scontrol >/dev/null 2>&1; then
    SELF="$(scontrol show job "${SLURM_JOB_ID}" 2>/dev/null \
            | awk -F'Command=' 'NF>1{print $2; exit}' | awk '{print $1}')"
fi
[[ -z "${SELF:-}" || ! -f "${SELF:-}" ]] && SELF="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="${DESEQ2_DGE_DIR:-$(cd "$(dirname "${SELF}")" >/dev/null 2>&1 && pwd)}"

STAGE_FILES=( "deseq2_01_setup.R" "deseq2_02_counts.R" "deseq2_03_model.R" "deseq2_04_results.R" )
STAGE_NAMES=( "1/setup" "2/counts" "3/model" "4/results" )

# ---------------------------------------------------------------------------
# Split --stage off; everything else is forwarded verbatim
# ---------------------------------------------------------------------------
STAGE_SPEC="all"
CONTRASTS_FILE=""
PASS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stage)   [[ $# -ge 2 ]] || { echo "ERROR: --stage needs a value."; exit 1; }
                   STAGE_SPEC="$2"; shift 2 ;;
        --stage=*) STAGE_SPEC="${1#--stage=}"; shift ;;
        --contrasts)   [[ $# -ge 2 ]] || { echo "ERROR: --contrasts needs a CSV file."; exit 1; }
                       CONTRASTS_FILE="$2"; shift 2 ;;
        --contrasts=*) CONTRASTS_FILE="${1#--contrasts=}"; shift ;;
        *)         PASS+=( "$1" ); shift ;;
    esac
done

show_help() {
    echo "Wrapper for the 4-stage DESeq2 pipeline. All flags pass through to the R scripts."
    echo "Submitter-only: --stage <1|2|3|4|all|N-|N-M>"
    echo
    if [[ -f "${SCRIPT_DIR}/${STAGE_FILES[0]}" ]]; then
        Rscript "${SCRIPT_DIR}/${STAGE_FILES[0]}" --help 2>/dev/null \
          || sed -n '1,80p' "${SCRIPT_DIR}/${STAGE_FILES[0]}"
    fi
}
if [[ ${#PASS[@]} -eq 0 && -z "${CONTRASTS_FILE}" ]] || printf '%s\n' "${PASS[@]}" | grep -qxE -- '--help|-h'; then
    show_help; exit 0
fi

for f in "${STAGE_FILES[@]}" "deseq2_common.R"; do
    if [[ ! -f "${SCRIPT_DIR}/${f}" ]]; then
        echo "ERROR: cannot find ${f}"
        echo "       looked in: ${SCRIPT_DIR}"
        echo "       Keep all five R files together, or set DESEQ2_DGE_DIR=/path/to/them"
        exit 1
    fi
done

# --stage -> FIRST..LAST
case "${STAGE_SPEC}" in
    all)              FIRST=1; LAST=4 ;;
    [1-4])            FIRST="${STAGE_SPEC}"; LAST="${STAGE_SPEC}" ;;
    [1-4]-)           FIRST="${STAGE_SPEC:0:1}"; LAST=4 ;;
    [1-4]-[1-4])      FIRST="${STAGE_SPEC:0:1}"; LAST="${STAGE_SPEC:2:1}" ;;
    *) echo "ERROR: --stage must be 1|2|3|4|all|N-|N-M (got '${STAGE_SPEC}')"; exit 1 ;;
esac
(( FIRST <= LAST )) || { echo "ERROR: --stage range is backwards ('${STAGE_SPEC}')"; exit 1; }

# --dry-run stops inside stage 1 anyway; do not run later stages against no state.
if printf '%s\n' "${PASS[@]}" | grep -qx -- '--dry-run'; then
    (( FIRST == 1 )) || { echo "ERROR: --dry-run only makes sense from stage 1."; exit 1; }
    LAST=1
    echo "NOTE: --dry-run set -> running stage 1 only (design validation)."
fi


# ---------------------------------------------------------------------------
# Conda environment. A function so that batch mode can also use it on the
# login node to validate every row before submitting anything.
# ---------------------------------------------------------------------------
CONDA_ENV=/lustre/isaac24/proj/UTK0009/Muhammad/melanoma_pertubation_rnaseq_cut-run_hi-c_multiomics_integration/conda-envs/multiomics-r
activate_env() {
    if [[ -n "${DESEQ2_SKIP_ENV:-}" ]]; then return 0; fi
    module purge
    module load anaconda3/2024.06
    unset PYTHONPATH                  # BEFORE activation
    export PYTHONNOUSERSITE=1
    # shellcheck disable=SC1090  # path provided by the module system at runtime
    source "$ANACONDA3_SH"
    conda activate "${CONDA_ENV}"
}

# ===========================================================================
# BATCH MODE (--contrasts): validate every row, then submit one job per row.
# Each child job is THIS script run normally for one comparison, so the
# single-comparison path below is the only code that ever runs an analysis.
# ===========================================================================
if [[ -n "${CONTRASTS_FILE}" ]]; then
    [[ -f "${CONTRASTS_FILE}" ]] || { echo "ERROR: --contrasts file not found: ${CONTRASTS_FILE}"; exit 1; }
    command -v sbatch >/dev/null 2>&1 || {
        echo "ERROR: sbatch not found. Run --contrasts on the cluster login node."; exit 1; }
    SELF_ABS="$(cd "$(dirname "${SELF}")" && pwd)/$(basename "${SELF}")"

    # Split off --outdir (the BASE folder) and --dry-run; refuse the per-row flags.
    BASE_OUT=""; BATCH_DRY=false; COMMON=()
    i=0
    while (( i < ${#PASS[@]} )); do
        a="${PASS[$i]}"
        case "$a" in
            --outdir)   BASE_OUT="${PASS[$((i+1))]:-}"; i=$((i+2)); continue ;;
            --outdir=*) BASE_OUT="${a#--outdir=}";      i=$((i+1)); continue ;;
            --dry-run)  BATCH_DRY=true;                 i=$((i+1)); continue ;;
            --control|--control=*|--treatment|--treatment=*|--fit-all-groups|--fit-all-groups=*|--contrast-name|--contrast-name=*)
                echo "ERROR: ${a%%=*} cannot be used with --contrasts -- the CSV supplies it per row."
                exit 1 ;;
        esac
        COMMON+=( "$a" ); i=$((i+1))
    done
    [[ -n "${BASE_OUT}" ]] || {
        echo "ERROR: --contrasts needs --outdir <BASE>. Row results go to <BASE>/<treatment>_vs_<control>/"
        exit 1; }
    mkdir -p "${BASE_OUT}"; BASE_OUT="$(cd "${BASE_OUT}" && pwd)"

    trim() {   # strip surrounding whitespace, then one pair of surrounding quotes
        local s="$1"
        s="${s#"${s%%[![:space:]]*}"}"; s="${s%"${s##*[![:space:]]}"}"
        s="${s#\"}"; s="${s%\"}"
        printf '%s' "$s"
    }

    # --- read the CSV (tolerates Excel's BOM, CRLF line endings, trailing commas)
    ROW_T=(); ROW_C=(); ROW_G=(); ROW_N=(); ROW_L=()
    lineno=0; nerr=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        lineno=$((lineno + 1))
        line="${line//$'\r'/}"
        if (( lineno == 1 )); then line="${line#$'\xEF\xBB\xBF'}"; fi
        if [[ -z "${line//[[:space:],]/}" || "$line" =~ ^[[:space:]]*# ]]; then continue; fi
        IFS=',' read -r t c gcol extra <<< "$line"
        t="$(trim "${t:-}")"; c="$(trim "${c:-}")"; gcol="$(trim "${gcol:-}")"
        # optional header row
        if (( ${#ROW_T[@]} == 0 )) && [[ "${t,,}" == "treatment" ]]; then continue; fi
        if [[ -n "${extra//[[:space:],]/}" ]]; then
            echo "ERROR: line ${lineno}: more than 3 columns: '${line}'"; nerr=$((nerr + 1)); continue; fi
        if [[ -z "$t" || -z "$c" ]]; then
            echo "ERROR: line ${lineno}: treatment and control are both required: '${line}'"; nerr=$((nerr + 1)); continue; fi
        if [[ "$t" == "$c" ]]; then
            echo "ERROR: line ${lineno}: treatment and control are the same ('${t}')"; nerr=$((nerr + 1)); continue; fi
        n="$(printf '%s_vs_%s' "$t" "$c" | tr -c 'A-Za-z0-9._-' '_')"
        for prev in ${ROW_N[@]+"${ROW_N[@]}"}; do
            if [[ "$prev" == "$n" ]]; then
                echo "ERROR: line ${lineno}: duplicate comparison '${n}' (each row needs its own folder)"
                nerr=$((nerr + 1)); continue 2
            fi
        done
        ROW_T+=( "$t" ); ROW_C+=( "$c" ); ROW_G+=( "$gcol" ); ROW_N+=( "$n" ); ROW_L+=( "$lineno" )
    done < "${CONTRASTS_FILE}"
    (( nerr == 0 )) || { echo "ERROR: ${nerr} problem(s) in ${CONTRASTS_FILE}; nothing submitted."; exit 1; }
    (( ${#ROW_T[@]} > 0 )) || { echo "ERROR: no comparisons found in ${CONTRASTS_FILE}."; exit 1; }

    row_args() {   # fills ARGS for row $1
        local k="$1"
        ARGS=( ${COMMON[@]+"${COMMON[@]}"} --outdir "${BASE_OUT}/${ROW_N[$k]}"
               --treatment "${ROW_T[$k]}" --control "${ROW_C[$k]}" )
        case "${ROW_G[$k]}" in
            "")            ;;                                   # contrast-only (default)
            [Aa][Ll][Ll])  ARGS+=( --fit-all-groups ) ;;        # bare: every sample
            *)             ARGS+=( --fit-all-groups "${ROW_G[$k]}" ) ;;
        esac
    }

    # --- validate EVERY row before submitting ANY -----------------------------
    echo "============================================================"
    echo "Batch mode: ${#ROW_T[@]} comparison(s) from ${CONTRASTS_FILE}"
    echo "Base outdir: ${BASE_OUT}"
    echo "============================================================"
    VLOG="${BASE_OUT}/_contrasts_validation"; mkdir -p "${VLOG}"
    nbad=0
    if (( FIRST == 1 )); then
        activate_env
        for k in "${!ROW_T[@]}"; do
            row_args "$k"
            if Rscript "${SCRIPT_DIR}/${STAGE_FILES[0]}" "${ARGS[@]}" --dry-run \
                   > "${VLOG}/${ROW_N[$k]}.log" 2>&1; then
                printf '  OK    %s\n' "${ROW_N[$k]}"
            else
                nbad=$((nbad + 1))
                printf '  FAIL  %s   (csv line %s)\n' "${ROW_N[$k]}" "${ROW_L[$k]}"
                grep -m 3 'ERROR' "${VLOG}/${ROW_N[$k]}.log" | sed 's/^/          /' || true
            fi
        done
    else
        # Resuming (--stage 2+): stage 1 is not re-run, so only check that each
        # row's earlier output exists. (A dry run here would reset its run.log.)
        for k in "${!ROW_T[@]}"; do
            if [[ -d "${BASE_OUT}/${ROW_N[$k]}/_state" ]]; then
                printf '  OK    %s\n' "${ROW_N[$k]}"
            else
                nbad=$((nbad + 1))
                printf '  FAIL  %s   (no earlier run in %s)\n' "${ROW_N[$k]}" "${BASE_OUT}/${ROW_N[$k]}"
            fi
        done
    fi
    if (( nbad > 0 )); then
        echo "ERROR: ${nbad} of ${#ROW_T[@]} row(s) failed validation. NOTHING was submitted."
        echo "       Details: ${VLOG}/"
        exit 1
    fi
    if $BATCH_DRY; then
        echo "--dry-run: all ${#ROW_T[@]} row(s) are valid. Nothing submitted."
        exit 0
    fi

    # --- submit: one job per row ---------------------------------------------
    MANIFEST="${BASE_OUT}/contrasts_submitted.tsv"
    printf 'job_id\tcomparison\ttreatment\tcontrol\tfit_group\toutdir\tsubmitted\n' > "${MANIFEST}"
    for k in "${!ROW_T[@]}"; do
        row_args "$k"
        jid="$(sbatch --parsable --job-name="dge_${ROW_N[$k]}" "${SELF_ABS}" \
                      --stage "${STAGE_SPEC}" "${ARGS[@]}")"
        jid="${jid%%;*}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${jid}" "${ROW_N[$k]}" "${ROW_T[$k]}" \
            "${ROW_C[$k]}" "${ROW_G[$k]:-<contrast-only>}" "${BASE_OUT}/${ROW_N[$k]}" "$(date '+%F %T')" \
            >> "${MANIFEST}"
        printf '  submitted %-10s %s\n' "${jid}" "${ROW_N[$k]}"
    done
    echo "Manifest : ${MANIFEST}"
    echo "Track    : squeue -u \$USER"
    exit 0
fi

echo "============================================================"
echo "submit_rnaseq_deseq2_dge.sh (v3, 4 stages)"
echo "Started  : $(date)"
echo "Node     : $(hostname)"
echo "CPUs     : ${SLURM_CPUS_PER_TASK:-NA}   Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "Scripts  : ${SCRIPT_DIR}"
echo "Stages   : ${FIRST}..${LAST}"
echo "Args     : ${PASS[*]}"
echo "============================================================"

# --- Environment hygiene (matches the rest of the project) ------------------
activate_env
export XDG_CACHE_HOME="${SLURM_SUBMIT_DIR:-$PWD}/.cache"
mkdir -p "${XDG_CACHE_HOME}"

# --- Run the stages in order ------------------------------------------------
# Stage 1 truncates run.log; stages 2-4 append. A failure stops the chain and
# names the stage, so you know whether the problem is the sample sheet, the
# annotation, the fit, or the contrast.
START_ALL=$SECONDS
for (( s=FIRST; s<=LAST; s++ )); do
    idx=$(( s - 1 ))
    echo
    echo "------------------------------------------------------------"
    echo ">> STAGE ${STAGE_NAMES[$idx]}  (${STAGE_FILES[$idx]})"
    echo "------------------------------------------------------------"
    t0=$SECONDS
    rc=0
    Rscript "${SCRIPT_DIR}/${STAGE_FILES[$idx]}" ${PASS[@]+"${PASS[@]}"} || rc=$?
    if (( rc != 0 )); then
        echo
        echo "============================================================"
        echo "FAILED at stage ${STAGE_NAMES[$idx]} (exit ${rc}) after $(( SECONDS - t0 ))s"
        echo "  The error above names the cause. See <outdir>/run.log for the"
        echo "  full record; earlier stages' output is intact, so once fixed you"
        echo "  can resume with:  --stage ${s}-"
        echo "============================================================"
        exit "${rc}"
    fi
    echo ">> stage ${STAGE_NAMES[$idx]} finished in $(( SECONDS - t0 ))s"
done

echo
echo "============================================================"
echo "All stages complete in $(( SECONDS - START_ALL ))s"
echo "Finished : $(date)"
echo "============================================================"
