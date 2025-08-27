#!/bin/bash
#SBATCH --job-name=7_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/7_analysis.err"
#SBATCH --output="./logs/7_analysis.out"

module load bedtools/2.30.0
module load deeptools/3.5.1
module load R/4.2.0

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
PEAKS_DIR="${BASE_DIR}/results/05_peaks"
BIGWIG_DIR="${BASE_DIR}/results/06_bigwig"
OUTPUT_DIR="${BASE_DIR}/results/07_analysis"

# 1. Find overlapping peaks between TES and TEAD1
echo "Finding overlapping peaks between TES and TEAD1..."
bedtools intersect \
    -a ${PEAKS_DIR}/TES_peaks.narrowPeak \
    -b ${PEAKS_DIR}/TEAD1_peaks.narrowPeak \
    -u > ${OUTPUT_DIR}/TES_TEAD1_overlap.bed

bedtools intersect \
    -a ${PEAKS_DIR}/TES_peaks.narrowPeak \
    -b ${PEAKS_DIR}/TEAD1_peaks.narrowPeak \
    -v > ${OUTPUT_DIR}/TES_only.bed

bedtools intersect \
    -a ${PEAKS_DIR}/TEAD1_peaks.narrowPeak \
    -b ${PEAKS_DIR}/TES_peaks.narrowPeak \
    -v > ${OUTPUT_DIR}/TEAD1_only.bed

# 2. Compare TES vs TESmut
echo "Comparing TES vs TESmut..."
bedtools intersect \
    -a ${PEAKS_DIR}/TES_peaks.narrowPeak \
    -b ${PEAKS_DIR}/TESmut_peaks.narrowPeak \
    -v > ${OUTPUT_DIR}/TES_not_in_TESmut.bed

# 3. Generate peak overlap statistics
wc -l ${PEAKS_DIR}/TES_peaks.narrowPeak | awk '{print "TES peaks: "$1}' > ${OUTPUT_DIR}/peak_stats.txt
wc -l ${PEAKS_DIR}/TESmut_peaks.narrowPeak | awk '{print "TESmut peaks: "$1}' >> ${OUTPUT_DIR}/peak_stats.txt
wc -l ${PEAKS_DIR}/TEAD1_peaks.narrowPeak | awk '{print "TEAD1 peaks: "$1}' >> ${OUTPUT_DIR}/peak_stats.txt
wc -l ${OUTPUT_DIR}/TES_TEAD1_overlap.bed | awk '{print "TES-TEAD1 overlapping peaks: "$1}' >> ${OUTPUT_DIR}/peak_stats.txt
wc -l ${OUTPUT_DIR}/TES_only.bed | awk '{print "TES-specific peaks: "$1}' >> ${OUTPUT_DIR}/peak_stats.txt
wc -l ${OUTPUT_DIR}/TEAD1_only.bed | awk '{print "TEAD1-specific peaks: "$1}' >> ${OUTPUT_DIR}/peak_stats.txt
wc -l ${OUTPUT_DIR}/TES_not_in_TESmut.bed | awk '{print "TES peaks not in TESmut: "$1}' >> ${OUTPUT_DIR}/peak_stats.txt

# 4. Create heatmaps and profiles
echo "Creating heatmaps..."

# Compute matrix for TES peaks
computeMatrix reference-point \
    --referencePoint center \
    -b 3000 -a 3000 \
    -R ${PEAKS_DIR}/TES_peaks.narrowPeak \
    -S ${BIGWIG_DIR}/TES-1_CPM.bw \
       ${BIGWIG_DIR}/TES-2_CPM.bw \
       ${BIGWIG_DIR}/TES-3_CPM.bw \
       ${BIGWIG_DIR}/TESmut-1_CPM.bw \
       ${BIGWIG_DIR}/TESmut-2_CPM.bw \
       ${BIGWIG_DIR}/TESmut-3_CPM.bw \
       ${BIGWIG_DIR}/TEAD1-1_CPM.bw \
       ${BIGWIG_DIR}/TEAD1-2_CPM.bw \
       ${BIGWIG_DIR}/TEAD1-3_CPM.bw \
    --skipZeros \
    -o ${OUTPUT_DIR}/TES_peaks_matrix.gz \
    -p 16

# Plot heatmap
plotHeatmap \
    -m ${OUTPUT_DIR}/TES_peaks_matrix.gz \
    -out ${OUTPUT_DIR}/TES_peaks_heatmap.pdf \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar' \
    --samplesLabel TES-1 TES-2 TES-3 TESmut-1 TESmut-2 TESmut-3 TEAD1-1 TEAD1-2 TEAD1-3

# Plot profile
plotProfile \
    -m ${OUTPUT_DIR}/TES_peaks_matrix.gz \
    -out ${OUTPUT_DIR}/TES_peaks_profile.pdf \
    --samplesLabel TES-1 TES-2 TES-3 TESmut-1 TESmut-2 TESmut-3 TEAD1-1 TEAD1-2 TEAD1-3 \
    --plotTitle "Signal at TES peaks"

echo "Analysis complete!"