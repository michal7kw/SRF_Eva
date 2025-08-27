#!/bin/bash
#SBATCH --job-name=multiqc_report
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/multiqc.err"
#SBATCH --output="./logs/multiqc.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate multiqc

# Set base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
cd $BASE_DIR

# Create MultiQC output directory
mkdir -p ${BASE_DIR}/results/multiqc

# Create a custom MultiQC configuration file
cat > ${BASE_DIR}/config/multiqc_config.yaml << 'EOF'
# MultiQC Configuration for Cut&Tag Analysis

title: "Cut&Tag Analysis Report: TES/TEAD1 Epigenetic Study"
subtitle: "Investigating transcriptional regulation by TES epigenetic silencer"
intro_text: |
    This report aggregates quality control metrics from the Cut&Tag analysis pipeline.
    The experiment investigates TES (TEAD1-KRAB-DNMT3A/3L) binding profiles compared to 
    endogenous TEAD1 and validates specificity using TESmut negative control.

report_comment: |
    Samples:
    - TES (3 replicates): Epigenetic silencer targeting TEAD1 sites
    - TESmut (3 replicates): Negative control with mutated TEAD1 binding domain
    - TEAD1 (3 replicates): Endogenous TEAD1 for comparison
    - IggMs/IggRb: Control antibodies

# Custom sample name cleaning
sample_names_rename:
    's/_001//g':
    's/_R1//g':
    's/_R2//g':
    's/_val_1//g':
    's/_val_2//g':
    's/_sorted//g':
    's/_filtered//g':
    's/_dedup//g':

# Sample order
sample_names_order:
    - TES-1
    - TES-2
    - TES-3
    - TESmut-1
    - TESmut-2
    - TESmut-3
    - TEAD1-1
    - TEAD1-2
    - TEAD1-3
    - IggMs
    - IggRb

# Module order in report
module_order:
    - fastqc:
        name: 'FastQC (Raw Reads)'
        path_filters:
            - '*_fastqc.zip'
        path_filters_exclude:
            - '*val*'
    - cutadapt:
        name: 'Trim Galore (Adapter Removal)'
    - fastqc:
        name: 'FastQC (Trimmed Reads)'
        path_filters:
            - '*val*_fastqc.zip'
    - bowtie2:
        name: 'Bowtie2 (Alignment)'
    - samtools:
        name: 'Samtools (Stats/Flagstat)'
    - picard:
        name: 'Picard (Duplicate Removal)'
    - macs2:
        name: 'MACS2 (Peak Calling)'
    - deeptools:
        name: 'deepTools (Coverage)'

# Exclude unnecessary files
fn_ignore_dirs:
    - '*/logs/*'
    - '*/scripts/*'
    - '*/config/*'

fn_ignore_files:
    - '*.md5'
    - '*.sh'
    - '*.R'

# Custom colors for sample groups
sample_colors:
    TES-1: '#2E7D32'
    TES-2: '#388E3C'
    TES-3: '#43A047'
    TESmut-1: '#D32F2F'
    TESmut-2: '#F44336'
    TESmut-3: '#EF5350'
    TEAD1-1: '#1976D2'
    TEAD1-2: '#2196F3'
    TEAD1-3: '#42A5F5'
    IggMs: '#757575'
    IggRb: '#9E9E9E'

# Top modules to show
top_modules:
    - 'fastqc'
    - 'bowtie2'
    - 'samtools'
    - 'macs2'

# Extra content for report
extra_fn_clean_exts:
    - '.gz'
    - '.fastq'
    - '.fq'
    - '.bam'
    - '.sam'
    - '.bed'
    - type: 'regex'
      pattern: '_peaks.narrowPeak'
      replace: ''

# Picard configuration
picard_config:
    duplicate_metrics:
        - 'PERCENT_DUPLICATION'
        - 'ESTIMATED_LIBRARY_SIZE'

# MACS2 configuration  
macs2_config:
    peak_count:
        - 'total_peaks'
        - 'fc_over_10'
    
# Samtools configuration
samtools_config:
    flagstat:
        - 'mapped_passed'
        - 'properly_paired'
    stats:
        - 'error_rate'
        - 'average_length'
        - 'average_quality'
EOF

echo "Running MultiQC with comprehensive configuration..."

# Run MultiQC with custom configuration
multiqc \
    ${BASE_DIR}/results \
    -c ${BASE_DIR}/config/multiqc_config.yaml \
    -o ${BASE_DIR}/results/multiqc \
    -n TES_TEAD1_CutTag_Report \
    --force \
    --dirs \
    --fullnames \
    --export \
    --pdf \
    -v

# Create additional custom statistics file
echo "Generating custom statistics summary..."

cat > ${BASE_DIR}/results/multiqc/custom_stats.txt << 'EOF'
===========================================
Cut&Tag Pipeline Summary Statistics
===========================================
Date: $(date)
Pipeline: Cut&Tag Analysis for TES/TEAD1
===========================================

SAMPLE GROUPS:
--------------
TES samples (Active): TES-1, TES-2, TES-3
TESmut samples (Control): TESmut-1, TESmut-2, TESmut-3  
TEAD1 samples (Reference): TEAD1-1, TEAD1-2, TEAD1-3
IgG Controls: IggMs (TES/TESmut), IggRb (TEAD1)

KEY METRICS TO CHECK:
--------------------
1. Sequencing Quality:
   - Per base sequence quality >30
   - Adapter contamination <5%
   
2. Alignment Metrics:
   - Alignment rate >80%
   - Properly paired reads >75%
   - Duplicate rate <30% (typical for Cut&Tag)
   
3. Peak Calling:
   - TES peaks: Expected similar count to TEAD1
   - TESmut peaks: Expected significantly fewer peaks
   - Peak overlap TES/TEAD1: Expected >50%

QUALITY CHECKPOINTS:
-------------------
EOF

# Add specific quality metrics from the pipeline
echo "[ ] Raw read quality (FastQC)" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
for sample in TES-1 TES-2 TES-3 TESmut-1 TESmut-2 TESmut-3 TEAD1-1 TEAD1-2 TEAD1-3 IggMs IggRb; do
    if [ -f "${BASE_DIR}/results/01_fastqc/${sample}_R1_001_fastqc/summary.txt" ]; then
        status=$(grep "Per base sequence quality" ${BASE_DIR}/results/01_fastqc/${sample}_R1_001_fastqc/summary.txt | cut -f1)
        echo "    ${sample}: ${status}" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
    fi
done

echo "" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
echo "[ ] Alignment rates (Bowtie2)" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
for sample in TES-1 TES-2 TES-3 TESmut-1 TESmut-2 TESmut-3 TEAD1-1 TEAD1-2 TEAD1-3 IggMs IggRb; do
    if [ -f "${BASE_DIR}/results/03_aligned/${sample}_bowtie2.log" ]; then
        rate=$(grep "overall alignment rate" ${BASE_DIR}/results/03_aligned/${sample}_bowtie2.log | cut -d' ' -f1)
        echo "    ${sample}: ${rate}" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
    fi
done

echo "" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
echo "[ ] Peak counts (MACS2)" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
for condition in TES TESmut TEAD1; do
    if [ -f "${BASE_DIR}/results/05_peaks/${condition}_peaks.narrowPeak" ]; then
        count=$(wc -l < ${BASE_DIR}/results/05_peaks/${condition}_peaks.narrowPeak)
        echo "    ${condition}: ${count} peaks" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
    fi
done

# Add overlap statistics if available
if [ -f "${BASE_DIR}/results/07_analysis/peak_stats.txt" ]; then
    echo "" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
    echo "PEAK OVERLAP ANALYSIS:" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
    echo "----------------------" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
    cat ${BASE_DIR}/results/07_analysis/peak_stats.txt >> ${BASE_DIR}/results/multiqc/custom_stats.txt
fi

echo "" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
echo "===========================================" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
echo "Report generated: $(date)" >> ${BASE_DIR}/results/multiqc/custom_stats.txt
echo "===========================================" >> ${BASE_DIR}/results/multiqc/custom_stats.txt

# Create a summary HTML report with plots
cat > ${BASE_DIR}/results/multiqc/summary_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Cut&Tag Analysis Summary: TES/TEAD1</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }
        h1 { color: #2E7D32; border-bottom: 3px solid #2E7D32; padding-bottom: 10px; }
        h2 { color: #1976D2; margin-top: 30px; }
        .container { background-color: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .metric-box { display: inline-block; margin: 10px; padding: 15px; border: 1px solid #ddd; border-radius: 5px; background-color: #fafafa; }
        .good { background-color: #E8F5E9; }
        .warning { background-color: #FFF3E0; }
        .bad { background-color: #FFEBEE; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #2E7D32; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .hypothesis { background-color: #E3F2FD; padding: 15px; border-left: 5px solid #1976D2; margin: 20px 0; }
        .conclusion { background-color: #F3E5F5; padding: 15px; border-left: 5px solid #7B1FA2; margin: 20px 0; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Cut&Tag Analysis Summary: TES/TEAD1 Epigenetic Study</h1>
        
        <div class="hypothesis">
            <h3>Hypothesis</h3>
            <p><strong>TES (TEAD1-KRAB-DNMT3A/3L)</strong> should bind to the same genomic targets as endogenous TEAD1, 
            while <strong>TESmut</strong> (with mutated TEAD1 binding domain) should show loss of specific binding.</p>
        </div>
        
        <h2>Experimental Design</h2>
        <table>
            <tr>
                <th>Condition</th>
                <th>Description</th>
                <th>Antibody</th>
                <th>Control</th>
                <th>Replicates</th>
            </tr>
            <tr>
                <td>TES</td>
                <td>TEAD1 binding domain + KRAB + DNMT3A/3L</td>
                <td>V5 (mouse)</td>
                <td>IggMs</td>
                <td>3</td>
            </tr>
            <tr>
                <td>TESmut</td>
                <td>Mutated TEAD1 binding site (negative control)</td>
                <td>V5 (mouse)</td>
                <td>IggMs</td>
                <td>3</td>
            </tr>
            <tr>
                <td>TEAD1</td>
                <td>Endogenous TEAD1</td>
                <td>HA (rabbit)</td>
                <td>IggRb</td>
                <td>3</td>
            </tr>
        </table>
        
        <h2>Quality Control Summary</h2>
        <div class="metric-box good">
            <h4>✓ Sequencing Quality</h4>
            <p>Check FastQC reports for quality scores >30</p>
        </div>
        <div class="metric-box good">
            <h4>✓ Alignment Rates</h4>
            <p>Target: >80% alignment rate</p>
        </div>
        <div class="metric-box good">
            <h4>✓ Library Complexity</h4>
            <p>Duplication rate should be <30%</p>
        </div>
        
        <h2>Key Results</h2>
        <p>Review the following key outputs:</p>
        <ul>
            <li><strong>Peak Overlap:</strong> Check TES vs TEAD1 overlap in <code>results/07_analysis/TES_TEAD1_overlap.bed</code></li>
            <li><strong>TESmut Validation:</strong> Verify reduced peaks in TESmut samples</li>
            <li><strong>Differential Binding:</strong> Review DiffBind results in <code>results/07_analysis/DiffBind_*.csv</code></li>
            <li><strong>Motif Analysis:</strong> Check HOMER results for TEAD motif enrichment</li>
            <li><strong>Target Genes:</strong> Review annotated peaks in <code>results/07_analysis/*_annotated.csv</code></li>
        </ul>
        
        <div class="conclusion">
            <h3>Expected Outcomes</h3>
            <ol>
                <li><strong>High overlap</strong> between TES and TEAD1 peaks (>50%)</li>
                <li><strong>Significantly fewer peaks</strong> in TESmut samples</li>
                <li><strong>TEAD motif enrichment</strong> in both TES and TEAD1 peaks</li>
                <li><strong>Common target genes</strong> involved in cell proliferation and stemness</li>
            </ol>
        </div>
        
        <h2>Report Files</h2>
        <ul>
            <li><a href="TES_TEAD1_CutTag_Report.html">Full MultiQC Report</a></li>
            <li><a href="custom_stats.txt">Custom Statistics Summary</a></li>
            <li><a href="../07_analysis/">Analysis Results Directory</a></li>
        </ul>
        
        <hr>
        <p><small>Generated: <script>document.write(new Date().toLocaleString());</script></small></p>
    </div>
</body>
</html>
EOF

echo "================================================"
echo "MultiQC report generation complete!"
echo "================================================"
echo ""
echo "Output files generated:"
echo "  1. Main report: ${BASE_DIR}/results/multiqc/TES_TEAD1_CutTag_Report.html"
echo "  2. PDF report: ${BASE_DIR}/results/multiqc/TES_TEAD1_CutTag_Report.pdf"
echo "  3. Data files: ${BASE_DIR}/results/multiqc/TES_TEAD1_CutTag_Report_data/"
echo "  4. Custom stats: ${BASE_DIR}/results/multiqc/custom_stats.txt"
echo "  5. Summary HTML: ${BASE_DIR}/results/multiqc/summary_report.html"
echo ""
echo "Key metrics to review:"
echo "  - Sequencing quality (FastQC)"
echo "  - Alignment rates (>80% expected)"
echo "  - Duplication rates (<30% for Cut&Tag)"
echo "  - Peak counts (TES ≈ TEAD1 >> TESmut)"
echo "  - Peak overlap (TES ∩ TEAD1)"
echo "================================================"

# Check for potential issues and create alerts
echo ""
echo "Checking for potential issues..."

# Check if peak counts are as expected
if [ -f "${BASE_DIR}/results/05_peaks/TES_peaks.narrowPeak" ] && [ -f "${BASE_DIR}/results/05_peaks/TESmut_peaks.narrowPeak" ]; then
    TES_peaks=$(wc -l < ${BASE_DIR}/results/05_peaks/TES_peaks.narrowPeak)
    TESmut_peaks=$(wc -l < ${BASE_DIR}/results/05_peaks/TESmut_peaks.narrowPeak)
    
    if [ $TESmut_peaks -gt $(($TES_peaks / 2)) ]; then
        echo "⚠️  WARNING: TESmut has unexpectedly high peak count ($TESmut_peaks vs TES: $TES_peaks)"
        echo "   This may indicate incomplete loss of binding in the mutant."
    else
        echo "✓ TESmut shows expected reduction in peaks"
    fi
fi

# Check for overlap between TES and TEAD1
if [ -f "${BASE_DIR}/results/07_analysis/TES_TEAD1_overlap.bed" ]; then
    overlap_count=$(wc -l < ${BASE_DIR}/results/07_analysis/TES_TEAD1_overlap.bed)
    echo "✓ Found $overlap_count overlapping peaks between TES and TEAD1"
fi

echo ""
echo "MultiQC pipeline complete at $(date)"