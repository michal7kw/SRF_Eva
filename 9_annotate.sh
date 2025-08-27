#!/bin/bash
#SBATCH --job-name=8_annotate
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/8_annotate.err"
#SBATCH --output="./logs/8_annotate.out"

module load R/4.2.0
module load homer/4.11

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
PEAKS_DIR="${BASE_DIR}/results/05_peaks"
ANALYSIS_DIR="${BASE_DIR}/results/07_analysis"

# R script for peak annotation using ChIPseeker
cat > ${BASE_DIR}/scripts/annotate_peaks.R << 'EOF'
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva")

# Load TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peak files
tes_peaks <- readPeakFile("results/05_peaks/TES_peaks.narrowPeak")
tesmut_peaks <- readPeakFile("results/05_peaks/TESmut_peaks.narrowPeak")
tead1_peaks <- readPeakFile("results/05_peaks/TEAD1_peaks.narrowPeak")

# Annotate peaks
tes_anno <- annotatePeak(tes_peaks, tssRegion=c(-3000, 3000), TxDb=txdb)
tesmut_anno <- annotatePeak(tesmut_peaks, tssRegion=c(-3000, 3000), TxDb=txdb)
tead1_anno <- annotatePeak(tead1_peaks, tssRegion=c(-3000, 3000), TxDb=txdb)

# Plot annotation statistics
pdf("results/07_analysis/peak_annotation_pie.pdf", width=12, height=4)
par(mfrow=c(1,3))
plotAnnoPie(tes_anno, main="TES Peak Distribution")
plotAnnoPie(tesmut_anno, main="TESmut Peak Distribution")
plotAnnoPie(tead1_anno, main="TEAD1 Peak Distribution")
dev.off()

# Get gene lists
tes_genes <- as.data.frame(tes_anno)$SYMBOL
tead1_genes <- as.data.frame(tead1_anno)$SYMBOL

# Find overlapping genes
common_genes <- intersect(tes_genes, tead1_genes)
tes_specific <- setdiff(tes_genes, tead1_genes)
tead1_specific <- setdiff(tead1_genes, tes_genes)

# Save gene lists
write.table(common_genes, "results/07_analysis/TES_TEAD1_common_genes.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(tes_specific, "results/07_analysis/TES_specific_genes.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# GO enrichment for TES target genes
ego <- enrichGO(gene = na.omit(tes_genes),
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

pdf("results/07_analysis/TES_GO_enrichment.pdf", width=10, height=8)
dotplot(ego, showCategory=20)
dev.off()

# Save results
write.csv(as.data.frame(tes_anno), "results/07_analysis/TES_peaks_annotated.csv")
write.csv(as.data.frame(tead1_anno), "results/07_analysis/TEAD1_peaks_annotated.csv")
EOF

# Run R script
Rscript ${BASE_DIR}/scripts/annotate_peaks.R

# HOMER motif analysis for TES peaks
echo "Running motif analysis..."
findMotifsGenome.pl \
    ${PEAKS_DIR}/TES_peaks.narrowPeak \
    hg38 \
    ${ANALYSIS_DIR}/TES_motifs \
    -size 200 \
    -mask \
    -p 8

# Compare with TEAD1 peaks
findMotifsGenome.pl \
    ${PEAKS_DIR}/TEAD1_peaks.narrowPeak \
    hg38 \
    ${ANALYSIS_DIR}/TEAD1_motifs \
    -size 200 \
    -mask \
    -p 8

echo "Annotation and motif analysis complete!"