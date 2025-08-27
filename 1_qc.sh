#!/bin/bash
#SBATCH --job-name=1_fastqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/1_fastqc_%a.err"
#SBATCH --output="./logs/1_fastqc_%a.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate quality

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
FASTQ_DIR="${BASE_DIR}/90-1222471453/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/01_fastqc"

# Array of samples
SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: ${SAMPLE}"

# Run FastQC on both R1 and R2
fastqc -t 4 -o ${OUTPUT_DIR} \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz

echo "FastQC complete for ${SAMPLE}"