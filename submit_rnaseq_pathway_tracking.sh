#!/bin/bash
#SBATCH -J gene_tracking
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=00-23:00:00
#SBATCH -A ACF-UTK0009
#SBATCH --partition=condo-rmccord
#SBATCH --qos=condo
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out

module load anaconda3/2024.06
source $ANACONDA3_SH
conda activate multiomics-r

R=/lustre/isaac24/proj/UTK0009/Muhammad/huvec_imr90_nl20_radiation_project
G=$R/data/refs/msigdb/h.all.v2026.1.Hs.symbols.gmt
FOCUS="CORE_DNA_REPAIR"
H=$R/results/HUVEC_IMR90_NL20_deseq2_outputs/HUVEC_2Gy24hr_vs_untreated
I=$R/results/HUVEC_IMR90_NL20_deseq2_outputs/IMR90_1Gy24hr_vs_untreated
N=/lustre/isaac24/proj/UTK0009/Muhammad/RNAseq/T2_NL20_600mGy24h_vs_untreated
OUT=$R/results/pathway_tracking_outputs
SCRIPT=$R/../RNASeq_pipeline_20AUG2026/rnaseq_pathway_tracking.R

# samples to drop, per cell line (empty = keep all)
declare -A EXCL=(
  [IMR90]="IMR90-100mGy-24hr-R2,IMR90-600mGy-5day-R3,IMR90-1Gy-5day-R3"
  [NL20]=""
)

# 1. reference (HUVEC): unchanged
Rscript $SCRIPT --deseq-dir $H --cell-name HUVEC --outdir $OUT/pathway_tracking_HUVEC \
    --genes-per-set 0 --gmt $G --focus-sets "$FOCUS" \
    --common-genes "$H,$I,$N" --export-order $OUT/pathway_tracking_row_order.tsv

# 2-3. followers: exclusion added only where a list is given
for X in "IMR90 $I" "NL20 $N"; do set -- $X
  EX=(); [[ -n "${EXCL[$1]}" ]] && EX=(--exclude-samples "${EXCL[$1]}")
  Rscript $SCRIPT --deseq-dir $2 --cell-name $1 --outdir $OUT/pathway_tracking_$1 \
      --genes-per-set 0 --gmt $G --focus-sets "$FOCUS" \
      --row-order $OUT/pathway_tracking_row_order.tsv "${EX[@]}"
done
