#!/bin/bash
#SBATCH --job-name=bigwig_cpm
#SBATCH --ntasks=1
#SBATCH --time=2-00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7Gb
#SBATCH --account=srrao
#SBATCH --mail-user=aray@mcw.edu

# BigWig Generation Script with CPM Normalization
# This script generates bigwig files from BAM files using deeptools

echo "=== BigWig Generation with CPM Normalization ==="
echo "Date: $(date)"
echo

# Load required modules
module load samtools
module load deeptools

# Set paths
BAM_DIR="henikoff_bam"
OUTPUT_DIR="henikoff_bigwig"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
echo "Created output directory: ${OUTPUT_DIR}"
echo

# Define BAM files
BAM_FILES=("H0.06K_K27me3.dedup.bam" "H0.6K_K27me3.dedup.bam" "H100K_IgG_R2.dedup.bam" "H2K_K27me3.dedup.bam" "H6K_K27me3.dedup.bam" "H0.2K_K27me3.dedup.bam" "H100K_IgG_R1.dedup.bam" "H20K_K27me3.dedup.bam" "H60K_K27me3.dedup.bam")

# Function to generate bigwig from BAM
generate_bigwig() {
    local bam_file=$1
    local base_name=$(basename "${bam_file}" .dedup.bam)
    local bam_path="${BAM_DIR}/${bam_file}"
    local bigwig_output="${OUTPUT_DIR}/${base_name}.CPM.bw"
    
    echo "=========================================="
    echo "Processing: ${bam_file}"
    echo "Output: ${bigwig_output}"
    echo "=========================================="
    
    # Check if BAM file exists
    if [[ ! -f "${bam_path}" ]]; then
        echo "ERROR: BAM file not found: ${bam_path}"
        return 1
    fi
    
    # Check if BAM file is indexed
    if [[ ! -f "${bam_path}.bai" ]]; then
        echo "Indexing BAM file..."
        samtools index "${bam_path}"
    fi
    
    # Get read count for logging
    local read_count=$(samtools view -c -F 4 "${bam_path}")
    echo "Mapped reads: ${read_count}"
    
    # Generate bigwig with CPM normalization
    echo "Generating bigwig file..."
    bamCoverage \
        -b "${bam_path}" \
        -o "${bigwig_output}" \
        --binSize 10 \
        --normalizeUsing CPM \
        --extendReads \
        --centerReads \
        --numberOfProcessors ${SLURM_CPUS_PER_TASK} \
        --verbose
    
    if [[ $? -eq 0 ]]; then
        echo "✓ Successfully generated: ${bigwig_output}"
    else
        echo "✗ Failed to generate bigwig for: ${bam_file}"
        return 1
    fi
    
    echo
}

# Process each BAM file
for bam in "${BAM_FILES[@]}"; do
    generate_bigwig "${bam}"
done

# Summary
echo "=========================================="
echo "BigWig Generation Complete!"
echo "=========================================="
echo "Output directory: ${OUTPUT_DIR}"
echo "Generated files:"
ls -lh "${OUTPUT_DIR}"/*.bw 2>/dev/null || echo "No bigwig files found"
echo
echo "Completed at: $(date)"


