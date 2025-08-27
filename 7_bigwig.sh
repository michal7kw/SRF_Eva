#!/bin/bash
#SBATCH --job-name=6_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6_bigwig_%a.err"
#SBATCH --output="./logs/6_bigwig_%a.out"

module load deeptools/3.5.1

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
BAM_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/06_bigwig"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Creating BigWig for ${SAMPLE}"

# Generate normalized BigWig files (CPM normalization)
bamCoverage \
    -b ${BAM_DIR}/${SAMPLE}_filtered.bam \
    -o ${OUTPUT_DIR}/${SAMPLE}_CPM.bw \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 8 \
    --extendReads

echo "BigWig generation complete for ${SAMPLE}"