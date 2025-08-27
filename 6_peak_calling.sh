#!/bin/bash
#SBATCH --job-name=5_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5_peaks.err"
#SBATCH --output="./logs/5_peaks.out"

module load macs2/2.2.7.1
module load bedtools/2.30.0

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
BAM_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/05_peaks"

# Peak calling for each condition with appropriate controls
echo "Calling peaks for TES samples..."
macs2 callpeak \
    -t ${BAM_DIR}/TES-1_filtered.bam \
       ${BAM_DIR}/TES-2_filtered.bam \
       ${BAM_DIR}/TES-3_filtered.bam \
    -c ${BAM_DIR}/IggMs_filtered.bam \
    -f BAMPE \
    -g hs \
    -n TES \
    --outdir ${OUTPUT_DIR} \
    -q 0.01 \
    --keep-dup all \
    2> ${OUTPUT_DIR}/TES_macs2.log

echo "Calling peaks for TESmut samples..."
macs2 callpeak \
    -t ${BAM_DIR}/TESmut-1_filtered.bam \
       ${BAM_DIR}/TESmut-2_filtered.bam \
       ${BAM_DIR}/TESmut-3_filtered.bam \
    -c ${BAM_DIR}/IggMs_filtered.bam \
    -f BAMPE \
    -g hs \
    -n TESmut \
    --outdir ${OUTPUT_DIR} \
    -q 0.01 \
    --keep-dup all \
    2> ${OUTPUT_DIR}/TESmut_macs2.log

echo "Calling peaks for TEAD1 samples..."
macs2 callpeak \
    -t ${BAM_DIR}/TEAD1-1_filtered.bam \
       ${BAM_DIR}/TEAD1-2_filtered.bam \
       ${BAM_DIR}/TEAD1-3_filtered.bam \
    -c ${BAM_DIR}/IggRb_filtered.bam \
    -f BAMPE \
    -g hs \
    -n TEAD1 \
    --outdir ${OUTPUT_DIR} \
    -q 0.01 \
    --keep-dup all \
    2> ${OUTPUT_DIR}/TEAD1_macs2.log

# Also call peaks for individual replicates for reproducibility analysis
for SAMPLE in TES-1 TES-2 TES-3 TESmut-1 TESmut-2 TESmut-3; do
    echo "Calling peaks for ${SAMPLE}..."
    macs2 callpeak \
        -t ${BAM_DIR}/${SAMPLE}_filtered.bam \
        -c ${BAM_DIR}/IggMs_filtered.bam \
        -f BAMPE \
        -g hs \
        -n ${SAMPLE} \
        --outdir ${OUTPUT_DIR} \
        -q 0.01 \
        --keep-dup all \
        2> ${OUTPUT_DIR}/${SAMPLE}_macs2.log
done

for SAMPLE in TEAD1-1 TEAD1-2 TEAD1-3; do
    echo "Calling peaks for ${SAMPLE}..."
    macs2 callpeak \
        -t ${BAM_DIR}/${SAMPLE}_filtered.bam \
        -c ${BAM_DIR}/IggRb_filtered.bam \
        -f BAMPE \
        -g hs \
        -n ${SAMPLE} \
        --outdir ${OUTPUT_DIR} \
        -q 0.01 \
        --keep-dup all \
        2> ${OUTPUT_DIR}/${SAMPLE}_macs2.log
done

echo "Peak calling complete!"