#!/bin/bash
#SBATCH --job-name=10_report
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/10_report.err"
#SBATCH --output="./logs/10_report.out"

module load multiqc/1.12

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"

# Run MultiQC to aggregate all QC reports
multiqc \
    ${BASE_DIR}/results \
    -o ${BASE_DIR}/results/multiqc \
    -n TES_TEAD1_CutTag_Report

# Create custom summary report
cat > ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'
# Cut&Tag Analysis Summary: TES/TEAD1 Epigenetic Study

## Experimental Design
- **TES**: Epigenetic silencer with TEAD1 binding domain + KRAB + DNMT3A/3L
- **TESmut**: Negative control with mutated TEAD1 binding site
- **TEAD1**: Endogenous TEAD1 for comparison
- **Controls**: IggMs (for TES/TESmut), IggRb (for TEAD1)

## Key Findings

### Peak Statistics
EOF

cat ${BASE_DIR}/results/07_analysis/peak_stats.txt >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md

cat >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'

## Conclusions
1. **TES binding specificity**: Compare overlap between TES and TEAD1 peaks
2. **TESmut validation**: Confirm loss of binding in mutant
3. **Target genes**: Analyze common targets between TES and TEAD1
4. **Epigenetic modifications**: Validate KRAB/DNMT activity at target sites

## Next Steps
1. Validate key targets by ChIP-qPCR
2. RNA-seq to measure transcriptional repression
3. Bisulfite sequencing for DNA methylation analysis
4. ATAC-seq to assess chromatin accessibility changes

## Output Files
- FastQC reports: `results/01_fastqc/`
- Trimmed reads: `results/02_trimmed/`
- Aligned BAMs: `results/03_aligned/`
- Filtered BAMs: `results/04_filtered/`
- Peak files: `results/05_peaks/`
- BigWig tracks: `results/06_bigwig/`
- Analysis results: `results/07_analysis/`
- MultiQC report: `results/multiqc/`
EOF

echo "Analysis pipeline complete! Check results/ANALYSIS_SUMMARY.md for summary."