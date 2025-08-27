#!/bin/bash
#SBATCH --job-name=3_align
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/3_align_%a.err"
#SBATCH --output="./logs/3_align_%a.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda  /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
TRIMMED_DIR="${BASE_DIR}/results/02_trimmed"
OUTPUT_DIR="${BASE_DIR}/results/03_aligned"

# Reference genome - adjust path as needed
GENOME_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/GRCh38"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Aligning sample: ${SAMPLE}"

# Bowtie2 alignment with parameters optimized for Cut&Tag
# --very-sensitive for Cut&Tag, --no-mixed --no-discordant for proper pairs only
bowtie2 \
    -p 32 \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    --rg-id "${SAMPLE}" \
    --rg "SM:${SAMPLE}" \
    --rg "LB:${SAMPLE}" \
    --rg "PL:ILLUMINA" \
    --rg "PU:${SAMPLE}" \
    -x ${GENOME_INDEX} \
    -1 ${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz \
    -2 ${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz \
    2> ${OUTPUT_DIR}/${SAMPLE}_bowtie2.log | \
    samtools view -@ 32 -bS - | \
    samtools sort -@ 32 -o ${OUTPUT_DIR}/${SAMPLE}_sorted.bam -

# Index the BAM file
samtools index -@ 32 ${OUTPUT_DIR}/${SAMPLE}_sorted.bam

# Generate alignment statistics
samtools flagstat ${OUTPUT_DIR}/${SAMPLE}_sorted.bam > ${OUTPUT_DIR}/${SAMPLE}_flagstat.txt

echo "Alignment complete for ${SAMPLE}"