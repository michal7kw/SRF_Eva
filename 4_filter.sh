#!/bin/bash
#SBATCH --job-name=4_filter
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/4_filter_%a.err"
#SBATCH --output="./logs/4_filter_%a.out"

# Function to handle errors
handle_error() {
    echo "Error on line $1"
    exit 1
}
trap 'handle_error $LINENO' ERR

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment

# Wait a bit to avoid conda lock conflicts
sleep $((RANDOM % 30))

# Install Picard if not already present (idempotent)
mamba install -y picard

# Find Picard.jar and set PICARD variable
PICARD_PATH=$(find /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment -name "picard.jar" | head -n 1)
if [ -z "$PICARD_PATH" ]; then
    echo "Error: picard.jar not found in the alignment conda environment."
    exit 1
fi
PICARD="$PICARD_PATH"
echo "PICARD set to: $PICARD"

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
ALIGNED_DIR="${BASE_DIR}/results/03_aligned"
OUTPUT_DIR="${BASE_DIR}/results/04_filtered"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Filtering sample: ${SAMPLE}"

# Check if input file exists
if [ ! -f "${ALIGNED_DIR}/${SAMPLE}_sorted.bam" ]; then
    echo "Error: Input file ${ALIGNED_DIR}/${SAMPLE}_sorted.bam not found"
    exit 1
fi

# Add read groups to BAM file (required for Picard MarkDuplicates)
echo "Adding read groups to BAM file..."
java -Xmx30g -jar $PICARD AddOrReplaceReadGroups \
    I=${ALIGNED_DIR}/${SAMPLE}_sorted.bam \
    O=${OUTPUT_DIR}/${SAMPLE}_with_rg.bam \
    RGID=${SAMPLE} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${SAMPLE} \
    CREATE_INDEX=true

# Check if read group addition was successful
if [ ! -f "${OUTPUT_DIR}/${SAMPLE}_with_rg.bam" ]; then
    echo "Error: Failed to add read groups"
    exit 1
fi

# Remove duplicates with Picard
echo "Removing duplicates..."
java -Xmx30g -jar $PICARD MarkDuplicates \
    I=${OUTPUT_DIR}/${SAMPLE}_with_rg.bam \
    O=${OUTPUT_DIR}/${SAMPLE}_dedup.bam \
    M=${OUTPUT_DIR}/${SAMPLE}_dedup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    CREATE_INDEX=true

# Check if deduplication was successful
if [ ! -f "${OUTPUT_DIR}/${SAMPLE}_dedup.bam" ]; then
    echo "Error: Deduplication failed"
    exit 1
fi

# Filter for properly paired, high-quality alignments
# Remove mitochondrial reads and blacklisted regions
echo "Filtering alignments..."
samtools view -@ 16 -b -q 30 -f 2 -F 1804 \
    ${OUTPUT_DIR}/${SAMPLE}_dedup.bam \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
    chr21 chr22 chrX chrY | \
    samtools sort -@ 16 -o ${OUTPUT_DIR}/${SAMPLE}_filtered.bam -

# Check if filtering was successful
if [ ! -f "${OUTPUT_DIR}/${SAMPLE}_filtered.bam" ]; then
    echo "Error: Filtering failed"
    exit 1
fi

# Index filtered BAM
echo "Indexing filtered BAM..."
samtools index -@ 16 ${OUTPUT_DIR}/${SAMPLE}_filtered.bam

# Generate final statistics
echo "Generating statistics..."
samtools flagstat ${OUTPUT_DIR}/${SAMPLE}_filtered.bam > ${OUTPUT_DIR}/${SAMPLE}_filtered_flagstat.txt

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm -f ${OUTPUT_DIR}/${SAMPLE}_with_rg.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_with_rg.bai
rm -f ${OUTPUT_DIR}/${SAMPLE}_dedup.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_dedup.bai

echo "Filtering complete for ${SAMPLE}"
echo "Final output: ${OUTPUT_DIR}/${SAMPLE}_filtered.bam"