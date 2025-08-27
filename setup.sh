#!/bin/bash
#SBATCH --job-name=setup_env
#SBATCH --account=kubacki.michal
#SBATCH --mem=4GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/setup.err"
#SBATCH --output="./logs/setup.out"

# Create directory structure
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva"
cd $BASE_DIR

mkdir -p scripts logs config
mkdir -p results/{01_fastqc,02_trimmed,03_aligned,04_filtered,05_peaks,06_bigwig,07_analysis}

# Create sample configuration file
cat > config/samples.txt << EOF
TES-1
TES-2
TES-3
TESmut-1
TESmut-2
TESmut-3
TEAD1-1
TEAD1-2
TEAD1-3
IggMs
IggRb
EOF

# Create group configuration
cat > config/groups.txt << EOF
TES:TES-1,TES-2,TES-3:IggMs
TESmut:TESmut-1,TESmut-2,TESmut-3:IggMs
TEAD1:TEAD1-1,TEAD1-2,TEAD1-3:IggRb
EOF

echo "Setup complete!"