#!/bin/bash
#SBATCH --job-name=9_diffbind
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/9_diffbind.err"
#SBATCH --output="./logs/9_diffbind.out"

module load R/4.2.0

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"

# Create DiffBind R script
cat > ${BASE_DIR}/scripts/diffbind_analysis.R << 'EOF'
library(DiffBind)
library(ggplot2)
library(pheatmap)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva")

# Create sample sheet
samples <- data.frame(
    SampleID = c("TES-1", "TES-2", "TES-3", 
                 "TESmut-1", "TESmut-2", "TESmut-3",
                 "TEAD1-1", "TEAD1-2", "TEAD1-3"),
    Tissue = rep("SNB19", 9),
    Factor = c(rep("TES", 3), rep("TESmut", 3), rep("TEAD1", 3)),
    Condition = c(rep("TES", 3), rep("TESmut", 3), rep("TEAD1", 3)),
    Replicate = rep(1:3, 3),
    bamReads = paste0("results/04_filtered/", 
                     c("TES-1", "TES-2", "TES-3", 
                       "TESmut-1", "TESmut-2", "TESmut-3",
                       "TEAD1-1", "TEAD1-2", "TEAD1-3"),
                     "_filtered.bam"),
    ControlID = c(rep("IggMs", 6), rep("IggRb", 3)),
    bamControl = c(rep("results/04_filtered/IggMs_filtered.bam", 6),
                  rep("results/04_filtered/IggRb_filtered.bam", 3)),
    Peaks = paste0("results/05_peaks/", 
                  c("TES-1", "TES-2", "TES-3", 
                    "TESmut-1", "TESmut-2", "TESmut-3",
                    "TEAD1-1", "TEAD1-2", "TEAD1-3"),
                  "_peaks.narrowPeak"),
    PeakCaller = rep("narrow", 9)
)

write.csv(samples, "config/diffbind_samplesheet.csv", row.names=FALSE)

# Load data
dba_obj <- dba(sampleSheet="config/diffbind_samplesheet.csv")

# Count reads
dba_obj <- dba.count(dba_obj, bParallel=TRUE)

# Normalize
dba_obj <- dba.normalize(dba_obj)

# Establish contrast: TES vs TESmut
dba_obj <- dba.contrast(dba_obj, 
                        group1=dba_obj$masks$TES, 
                        group2=dba_obj$masks$TESmut,
                        name1="TES", name2="TESmut")

# Establish contrast: TES vs TEAD1  
dba_obj <- dba.contrast(dba_obj,
                        group1=dba_obj$masks$TES,
                        group2=dba_obj$masks$TEAD1,
                        name1="TES", name2="TEAD1")

# Perform differential analysis
dba_obj <- dba.analyze(dba_obj)

# Get results
res_tes_vs_tesmut <- dba.report(dba_obj, contrast=1, th=1)
res_tes_vs_tead1 <- dba.report(dba_obj, contrast=2, th=1)

# Save results
write.csv(as.data.frame(res_tes_vs_tesmut), 
          "results/07_analysis/DiffBind_TES_vs_TESmut.csv")
write.csv(as.data.frame(res_tes_vs_tead1), 
          "results/07_analysis/DiffBind_TES_vs_TEAD1.csv")

# Plot MA
pdf("results/07_analysis/DiffBind_MA_plots.pdf", width=12, height=6)
par(mfrow=c(1,2))
dba.plotMA(dba_obj, contrast=1)
dba.plotMA(dba_obj, contrast=2)
dev.off()

# Plot PCA
pdf("results/07_analysis/DiffBind_PCA.pdf", width=8, height=8)
dba.plotPCA(dba_obj, label=DBA_ID)
dev.off()

# Plot correlation heatmap
pdf("results/07_analysis/DiffBind_correlation.pdf", width=10, height=10)
dba.plotHeatmap(dba_obj, correlations=TRUE)
dev.off()

# Venn diagram of overlaps
pdf("results/07_analysis/DiffBind_venn.pdf", width=8, height=8)
dba.plotVenn(dba_obj, mask=dba_obj$masks$Consensus)
dev.off()

print("DiffBind analysis complete!")
EOF

# Run DiffBind analysis
Rscript ${BASE_DIR}/scripts/diffbind_analysis.R

echo "Differential binding analysis complete!"