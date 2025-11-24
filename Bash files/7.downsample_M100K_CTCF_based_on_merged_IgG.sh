#!/bin/bash

#SBATCH --job-name=M100K_CTCF_downsample # Job name
#SBATCH --ntasks=1 # number of tasks
#SBATCH --time=5-20:00 # time limit (day-hour:min)
#SBATCH --cpus-per-task=15 # number of threads
#SBATCH --mem-per-cpu=7Gb # requested memory
#SBATCH --account=srrao # PI's net ID
#SBATCH --mail-user=aray@mcw.edu # Email address to send to

module load macs2
module load samtools
module load bedtools2
module load R

# Directories and file pattern
BAM_DIR="/scratch/g/srrao/aray/cell_no_std/dedup_bams"
MERGED_BAM_DIR="/scratch/g/srrao/aray/cell_no_std/merged_IgG"
OUT_DIR="downsampled_peak_results_for_M100K_CTCF_down_based_merge_M100K_IgG"
PEAK_DIR_DEFAULT="${OUT_DIR}/macs2_peaks_default"
PEAK_DIR_LARGE="${OUT_DIR}/macs2_peaks_large"
LOG_DIR="${OUT_DIR}/macs2_logs"
SUMMARY="${OUT_DIR}/downsample_summary.csv"

# Create output directories
mkdir -p "$OUT_DIR" "$PEAK_DIR_DEFAULT" "$PEAK_DIR_LARGE" "$LOG_DIR"

# Downsample targets for CTCF
targets=(5000000 2500000 1000000 500000 100000 50000 10000)

# Write header for summary
echo "Sample,replicate,target_reads_for_CTCF,actual_reads_for_CTCF,target_reads_for_IgG,actual_reads_for_IgG,num_peaks_default,num_peaks_large" > "$SUMMARY"

# Get the IgG merged BAM file
IGG_BAM="${MERGED_BAM_DIR}/M100K_IgG.merged.sorted.bam" 
if [ ! -f "$IGG_BAM" ]; then
    echo "ERROR: IgG merged BAM file not found: $IGG_BAM"
    exit 1
fi
echo "Using IgG merged BAM: $IGG_BAM"

# Get total reads in IgG merged BAM (NO DOWNSAMPLING)
igg_total_reads=$(samtools view -c -F 4 "$IGG_BAM")
echo "IgG merged BAM total reads: $igg_total_reads (no downsampling)"

# ===== Analyze original CTCF files with original IgG merged BAM =====
echo "Analyzing original CTCF files with original IgG merged BAM..."

# Loop through all three CTCF replicates (R1, R2, R3)
for ctcf_replicate in "R1" "R2" "R3"; do
    ctcf_bam="${BAM_DIR}/M100K_CTCF_${ctcf_replicate}.dedup.bam" 
    if [ ! -f "$ctcf_bam" ]; then
        echo "WARNING: CTCF file not found: $ctcf_bam, skipping..."
        continue
    fi
    echo "Processing CTCF replicate: $ctcf_replicate (original, no downsampling)"
    # Get CTCF total reads (original)
    ctcf_total_reads=$(samtools view -c -F 4 "$ctcf_bam")
    actual_ctcf_reads=$ctcf_total_reads
    echo "CTCF ${ctcf_replicate} reads (original): ${actual_ctcf_reads}"
    echo "IgG merged reads (original, no downsampling): ${igg_total_reads}"
    # Define peak file names for original files
    peak_file_orig_default="${PEAK_DIR_DEFAULT}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_default_peaks.narrowPeak"
    peak_file_orig_large="${PEAK_DIR_LARGE}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_large_peaks.narrowPeak"
    log_file_orig_default="${LOG_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_default.log"
    log_file_orig_large="${LOG_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_large.log"
    # Call peaks with macs2 - Default scale (q-value 0.01, --nolambda)
    echo "Calling peaks with macs2 default scale for M100K_CTCF_${ctcf_replicate} (IgG merged control, original)"
    macs2 callpeak \
        -t "$ctcf_bam" \
        -c "$IGG_BAM" \
        -n "M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_default" \
        -g mm \
        -f BAMPE \
        --outdir "$PEAK_DIR_DEFAULT" \
        -q 0.01 \
        --keep-dup all \
        --nomodel \
        --shift -75 \
        --extsize 150 \
        --nolambda \
        --verbose 2 \
        2> "$log_file_orig_default"
    # Call peaks with macs2 - Scale to large (q-value 0.00001, --nolambda)
    echo "Calling peaks with macs2 scale to large for M100K_CTCF_${ctcf_replicate} (IgG merged control, original)"
    macs2 callpeak \
        -t "$ctcf_bam" \
        -c "$IGG_BAM" \
        -n "M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_large" \
        -g mm \
        -f BAMPE \
        --outdir "$PEAK_DIR_LARGE" \
        -q 0.00001 \
        --keep-dup all \
        --nomodel \
        --shift -75 \
        --extsize 150 \
        --scale-to large \
        --nolambda \
        --verbose 2 \
        2> "$log_file_orig_large"
    # Wait for files to be written
    sleep 3
    # Count number of peaks for original files
    if [ -f "$peak_file_orig_default" ]; then
        num_peaks_orig_default=$(wc -l < "$peak_file_orig_default")
        echo "Found ${num_peaks_orig_default} peaks in ${peak_file_orig_default}"
    else
        num_peaks_orig_default=0
        echo "WARNING: Peak file not found: ${peak_file_orig_default}"
    fi
    if [ -f "$peak_file_orig_large" ]; then
        num_peaks_orig_large=$(wc -l < "$peak_file_orig_large")
        echo "Found ${num_peaks_orig_large} peaks in ${peak_file_orig_large}"
    else
        num_peaks_orig_large=0
        echo "WARNING: Peak file not found: ${peak_file_orig_large}"
    fi
    # Write original results to summary
    echo "M100K_CTCF_${ctcf_replicate},${ctcf_replicate},original,${ctcf_total_reads},original,${igg_total_reads},${num_peaks_orig_default},${num_peaks_orig_large}" >> "$SUMMARY"
    echo "Original results for CTCF_${ctcf_replicate} - Default: ${num_peaks_orig_default} peaks, Large: ${num_peaks_orig_large} peaks"
done

# ===== END Original Analysis Section =====

# ===== Process downsampled versions of CTCF (IgG stays original) =====
echo "Processing downsampled versions of CTCF (IgG stays original, no downsampling)..."
for target in "${targets[@]}"; do
    echo "Processing target: ${target} reads for CTCF"
    # Loop through all three CTCF replicates (R1, R2, R3)
    for ctcf_replicate in "R1" "R2" "R3"; do
        ctcf_bam="${BAM_DIR}/M100K_CTCF_${ctcf_replicate}.dedup.bam"
        if [ ! -f "$ctcf_bam" ]; then
            echo "WARNING: CTCF file not found: $ctcf_bam, skipping..."
            continue
        fi
        echo "Processing CTCF replicate: $ctcf_replicate downsampled to ${target} reads (IgG original, no downsampling)"
        # Get CTCF total reads (original, before downsampling)
        ctcf_total_reads=$(samtools view -c -F 4 "$ctcf_bam")
        # Calculate downsampling fraction for CTCF
        ctcf_fraction=$(awk "BEGIN {printf \"%.6f\", $target/$ctcf_total_reads}")
        # Check if downsampling is needed (if target >= total reads, skip downsampling)
        if (( $(echo "$target >= $ctcf_total_reads" | bc -l) )); then
            echo "Target ${target} >= total reads ${ctcf_total_reads}, using original file without downsampling"
            out_ctcf_bam="$ctcf_bam"
            actual_ctcf_reads=$ctcf_total_reads
        else
            # Use a unique seed for each target and replicate for reproducibility
            ctcf_seed=$((RANDOM % 10000 + 1))
            out_ctcf_bam="${OUT_DIR}/M100K_CTCF_${ctcf_replicate}.down${target}.bam"
            # Downsample CTCF BAM (samtools -s option uses seed.fraction format)
            echo "Downsampling CTCF ${ctcf_replicate} to ${target} reads (fraction: ${ctcf_fraction})"
            samtools view -s ${ctcf_seed}${ctcf_fraction} -b "$ctcf_bam" > "$out_ctcf_bam"
            # Get actual number of reads in downsampled CTCF BAM
            actual_ctcf_reads=$(samtools view -c -F 4 "$out_ctcf_bam")
            echo "Actual reads in downsampled CTCF file: ${actual_ctcf_reads}"
        fi
        echo "CTCF ${ctcf_replicate} reads (downsampled to target ${target}): ${actual_ctcf_reads}"
        echo "IgG merged reads (original, no downsampling): ${igg_total_reads}"
        # Define peak file names
        peak_file_default="${PEAK_DIR_DEFAULT}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_default_peaks.narrowPeak"
        peak_file_large="${PEAK_DIR_LARGE}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_large_peaks.narrowPeak"
        log_file_default="${LOG_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_default.log"
        log_file_large="${LOG_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_large.log"
        # Call peaks with macs2 - Default scale (q-value 0.01, --nolambda)
        echo "Calling peaks with macs2 default scale for M100K_CTCF_${ctcf_replicate} (downsampled to ${target}) vs IgG_merged (original)"
        macs2 callpeak \
            -t "$out_ctcf_bam" \
            -c "$IGG_BAM" \
            -n "M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_default" \
            -g mm \
            -f BAMPE \
            --outdir "$PEAK_DIR_DEFAULT" \
            -q 0.01 \
            --keep-dup all \
            --nomodel \
            --shift -75 \
            --extsize 150 \
            --nolambda \
            --verbose 2 \
            2> "$log_file_default"
        # Call peaks with macs2 - Scale to large (q-value 0.00001, --nolambda)
        echo "Calling peaks with macs2 scale to large for M100K_CTCF_${ctcf_replicate} (downsampled to ${target}) vs IgG_merged (original)"
        macs2 callpeak \
            -t "$out_ctcf_bam" \
            -c "$IGG_BAM" \
            -n "M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_large" \
            -g mm \
            -f BAMPE \
            --outdir "$PEAK_DIR_LARGE" \
            -q 0.00001 \
            --keep-dup all \
            --nomodel \
            --shift -75 \
            --extsize 150 \
            --scale-to large \
            --nolambda \
            --verbose 2 \
            2> "$log_file_large"
        # Wait for files to be written
        sleep 5
        # Count number of peaks for default scale
        if [ -f "$peak_file_default" ]; then
            num_peaks_default=$(wc -l < "$peak_file_default")
            echo "Found ${num_peaks_default} peaks in ${peak_file_default}"
        else
            num_peaks_default=0
            echo "WARNING: Peak file not found: ${peak_file_default}"
        fi
        # Count number of peaks for large scale
        if [ -f "$peak_file_large" ]; then
            num_peaks_large=$(wc -l < "$peak_file_large")
            echo "Found ${num_peaks_large} peaks in ${peak_file_large}"
        else
            num_peaks_large=0
            echo "WARNING: Peak file not found: ${peak_file_large}"
        fi
        # Write to summary
        echo "M100K_CTCF_${ctcf_replicate},${ctcf_replicate},${target},${actual_ctcf_reads},original,${igg_total_reads},${num_peaks_default},${num_peaks_large}" >> "$SUMMARY"
        echo "CTCF_${ctcf_replicate} (downsampled to ${target}) vs IgG (original): Default: ${num_peaks_default} peaks, Large: ${num_peaks_large} peaks"
        # Clean up downsampled CTCF BAM to save space (if it was created)
        if [ "$out_ctcf_bam" != "$ctcf_bam" ] && [ -f "$out_ctcf_bam" ]; then
            rm -f "$out_ctcf_bam"
        fi
    done
done

echo "All done! Summary written to $SUMMARY"
echo "Default scale peaks (q-value 0.01) in: $PEAK_DIR_DEFAULT"
echo "Large scale peaks (q-value 0.00001) in: $PEAK_DIR_LARGE"
echo "Peak files are preserved for further analysis"
echo "Original IgG merged BAM was used for all analyses (no downsampling)"
echo "CTCF files were downsampled based on target read counts"
echo "MACs2 logs are in: $LOG_DIR"
