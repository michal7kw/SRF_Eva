#!/bin/bash
#SBATCH --job-name=2_trim
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/2_trim_%a.err"
#SBATCH --output="./logs/2_trim_%a.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate trim

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
FASTQ_DIR="${BASE_DIR}/90-1222471453/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/02_trimmed"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Trimming sample: ${SAMPLE}"

# Trim Galore for paired-end Cut&Tag data
# Using standard Illumina adapters and stringent quality trimming
trim_galore \
    --paired \
    --cores 8 \
    --quality 20 \
    --stringency 5 \
    --length 20 \
    --fastqc \
    --gzip \
    --output_dir ${OUTPUT_DIR} \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz

echo "Trimming complete for ${SAMPLE}"