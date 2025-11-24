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
OUT_DIR="downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF"
PEAK_DIR_DEFAULT="${OUT_DIR}/macs2_peaks_default"
PEAK_DIR_LARGE="${OUT_DIR}/macs2_peaks_large"
LOG_DIR="${OUT_DIR}/macs2_logs"
SEACR_DIR="${OUT_DIR}/seacr_peaks"
BEDGRAPH_DIR="${OUT_DIR}/bedgraph_files"
SUMMARY="${OUT_DIR}/downsample_summary.csv"

# SEACR path
SEACR_PATH="/scratch/g/srrao/aray/seacr/SEACR_1.3.sh"

# Create output directories
mkdir -p "$OUT_DIR" "$PEAK_DIR_DEFAULT" "$PEAK_DIR_LARGE" "$LOG_DIR" "$SEACR_DIR" "$BEDGRAPH_DIR"

# Downsample targets
targets=(5000000 2500000 1000000 500000 100000 50000 10000)

# Write header for summary with CORRECTED columns
echo "Sample,replicate,target_reads_for_CTCF,actual_reads_for_CTCF,target_reads_for_IgG,actual_reads_for_IgG,num_peaks_default,num_peaks_large,num_peaks_seacr_stringent,num_peaks_seacr_relaxed" > "$SUMMARY"

# Function to convert BAM to bedgraph with normalization
bam_to_bedgraph() {
    local bam_file=$1
    local output_prefix=$2
    echo "Converting ${bam_file} to bedgraph..."
    # Get mapped reads only for more accurate normalization
    total_reads=$(samtools view -c -F 4 ${bam_file})
    scale_factor=$(echo "1000000/${total_reads}" | bc -l)
    echo "Total mapped reads: ${total_reads}"
    echo "Scale factor: ${scale_factor}"
    # Convert BAM to bedgraph with normalization
    bedtools genomecov \
        -ibam ${bam_file} \
        -bg \
        -pc \
        -scale ${scale_factor} \
        | sort -k1,1 -k2,2n \
        > "${BEDGRAPH_DIR}/${output_prefix}.bedgraph"
    echo "Created bedgraph file: ${BEDGRAPH_DIR}/${output_prefix}.bedgraph"
}

# Function to run SEACR
run_seacr() {
    local treatment_bg=$1
    local control_bg=$2
    local prefix=$3
    echo "Running SEACR analysis for ${prefix}"
    # 1. Stringent with IgG control
    bash ${SEACR_PATH} \
        ${treatment_bg} \
        ${control_bg} \
        norm \
        stringent \
        "${SEACR_DIR}/${prefix}_stringent"
    # 2. Relaxed with IgG control
    bash ${SEACR_PATH} \
        ${treatment_bg} \
        ${control_bg} \
        norm \
        relaxed \
        "${SEACR_DIR}/${prefix}_relaxed"
}

# Get the IgG merged BAM file
IGG_BAM="${MERGED_BAM_DIR}/M100K_IgG.merged.sorted.bam"
if [ ! -f "$IGG_BAM" ]; then
    echo "ERROR: IgG merged BAM file not found: $IGG_BAM"
    exit 1
fi
echo "Using IgG merged BAM: $IGG_BAM"

# Get total reads in IgG merged BAM
igg_total_reads=$(samtools view -c -F 4 "$IGG_BAM")
echo "IgG merged BAM total reads: $igg_total_reads"

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
    # Define peak file names for original files (matching macs2 output naming)
    peak_file_orig_default="${PEAK_DIR_DEFAULT}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_default_peaks.narrowPeak"
    peak_file_orig_large="${PEAK_DIR_LARGE}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_large_peaks.narrowPeak"
    peak_file_orig_seacr_stringent="${SEACR_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_seacr_stringent.bed"
    peak_file_orig_seacr_relaxed="${SEACR_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original_seacr_relaxed.bed"
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
        --verbose 2 \
        2> "$log_file_orig_large"
    # Convert BAMs to bedgraph for SEACR
    echo "Converting BAMs to bedgraph for SEACR analysis..."
    bam_to_bedgraph "$ctcf_bam" "M100K_CTCF_${ctcf_replicate}_original"
    bam_to_bedgraph "$IGG_BAM" "M100K_IgG_merged_original"
    # Run SEACR
    run_seacr \
        "${BEDGRAPH_DIR}/M100K_CTCF_${ctcf_replicate}_original.bedgraph" \
        "${BEDGRAPH_DIR}/M100K_IgG_merged_original.bedgraph" \
        "M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_original"
    # Wait for files to be written
    sleep 3
    # Count number of peaks for original files
    if [ -f "$peak_file_orig_default" ]; then
        num_peaks_orig_default=$(wc -l < "$peak_file_orig_default" | awk '{print $1}')
        echo "Found ${num_peaks_orig_default} peaks in ${peak_file_orig_default}"
    else
        num_peaks_orig_default=0
        echo "WARNING: Peak file not found: ${peak_file_orig_default}"
    fi
    if [ -f "$peak_file_orig_large" ]; then
        num_peaks_orig_large=$(wc -l < "$peak_file_orig_large" | awk '{print $1}')
        echo "Found ${num_peaks_orig_large} peaks in ${peak_file_orig_large}"
    else
        num_peaks_orig_large=0
        echo "WARNING: Peak file not found: ${peak_file_orig_large}"
    fi
    if [ -f "$peak_file_orig_seacr_stringent" ]; then
        num_peaks_orig_seacr_stringent=$(wc -l < "$peak_file_orig_seacr_stringent" | awk '{print $1}')
        echo "Found ${num_peaks_orig_seacr_stringent} peaks in ${peak_file_orig_seacr_stringent}"
    else
        num_peaks_orig_seacr_stringent=0
        echo "WARNING: SEACR stringent file not found: ${peak_file_orig_seacr_stringent}"
    fi
    if [ -f "$peak_file_orig_seacr_relaxed" ]; then
        num_peaks_orig_seacr_relaxed=$(wc -l < "$peak_file_orig_seacr_relaxed" | awk '{print $1}')
        echo "Found ${num_peaks_orig_seacr_relaxed} peaks in ${peak_file_orig_seacr_relaxed}"
    else
        num_peaks_orig_seacr_relaxed=0
        echo "WARNING: SEACR relaxed file not found: ${peak_file_orig_seacr_relaxed}"
    fi
    # Write original results to summary
    echo "M100K_CTCF_${ctcf_replicate},${ctcf_replicate},original,${ctcf_total_reads},original,${igg_total_reads},${num_peaks_orig_default},${num_peaks_orig_large},${num_peaks_orig_seacr_stringent},${num_peaks_orig_seacr_relaxed}" >> "$SUMMARY"
    echo "Original results for CTCF_${ctcf_replicate} - Default: ${num_peaks_orig_default} peaks, Large: ${num_peaks_orig_large} peaks, SEACR Stringent: ${num_peaks_orig_seacr_stringent} peaks, SEACR Relaxed: ${num_peaks_orig_seacr_relaxed} peaks"
done

# ===== END Original Analysis Section =====

# ===== Process downsampled versions of ONLY IgG (CTCF stays original) =====
echo "Processing downsampled versions of ONLY IgG (CTCF stays original, no downsampling)..."
for target in "${targets[@]}"; do
    echo "Processing target: ${target} reads for IgG"
    # Loop through all three CTCF replicates (R1, R2, R3)
    for ctcf_replicate in "R1" "R2" "R3"; do
        ctcf_bam="${BAM_DIR}/M100K_CTCF_${ctcf_replicate}.dedup.bam"
        if [ ! -f "$ctcf_bam" ]; then
            echo "WARNING: CTCF file not found: $ctcf_bam, skipping..."
            continue
        fi
        echo "Processing CTCF replicate: $ctcf_replicate with IgG downsampled to ${target} reads"
        # Get CTCF total reads (NO DOWNSAMPLING)
        ctcf_total_reads=$(samtools view -c -F 4 "$ctcf_bam")
        actual_ctcf_reads=$ctcf_total_reads
        echo "CTCF ${ctcf_replicate} reads (no downsampling): ${actual_ctcf_reads}"
        # Calculate downsampling fraction for IgG ONLY
        igg_fraction=$(awk "BEGIN {printf \"%.6f\", $target/$igg_total_reads}")
        # Check if downsampling is needed (if target >= total reads, skip downsampling)
        if (( $(echo "$target >= $igg_total_reads" | bc -l) )); then
            echo "Target ${target} >= total reads ${igg_total_reads}, using original file without downsampling"
            out_igg_bam="$IGG_BAM"
            actual_igg_reads=$igg_total_reads
        else
            # Use a unique seed for each target and replicate for reproducibility
            igg_seed=$((RANDOM % 10000 + 1))
            out_igg_bam="${OUT_DIR}/M100K_IgG_merged.down${target}.bam"
            # Downsample IgG merged BAM ONLY (samtools -s option uses seed.fraction format)
            echo "Downsampling IgG merged to ${target} reads (fraction: ${igg_fraction})"
            samtools view -s ${igg_seed}${igg_fraction} -b "$IGG_BAM" > "$out_igg_bam"
            # Get actual number of reads in downsampled IgG merged BAM
            actual_igg_reads=$(samtools view -c -F 4 "$out_igg_bam")
            echo "Actual reads in downsampled IgG file: ${actual_igg_reads}"
        fi
        echo "IgG merged reads (downsampled to target ${target}): ${actual_igg_reads}"
        # Define peak file names (matching macs2 output naming)
        peak_file_default="${PEAK_DIR_DEFAULT}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_default_peaks.narrowPeak"
        peak_file_large="${PEAK_DIR_LARGE}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_large_peaks.narrowPeak"
        peak_file_seacr_stringent="${SEACR_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_seacr_stringent.bed"
        peak_file_seacr_relaxed="${SEACR_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_seacr_relaxed.bed"
        log_file_default="${LOG_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_default.log"
        log_file_large="${LOG_DIR}/M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}_large.log"
        # Convert BAMs to bedgraph for SEACR (CTCF original, IgG downsampled)
        bam_to_bedgraph "$ctcf_bam" "M100K_CTCF_${ctcf_replicate}_original"
        bam_to_bedgraph "$out_igg_bam" "M100K_IgG_merged_${target}"
        # Call peaks with macs2 - Default scale (q-value 0.01, --nolambda)
        echo "Calling peaks with macs2 default scale for M100K_CTCF_${ctcf_replicate} (downsampled IgG_merged to ${target}) vs original CTCF"
        macs2 callpeak \
            -t "$ctcf_bam" \
            -c "$out_igg_bam" \
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
        echo "Calling peaks with macs2 scale to large for M100K_CTCF_${ctcf_replicate} (downsampled IgG_merged to ${target}) vs original CTCF"
        macs2 callpeak \
            -t "$ctcf_bam" \
            -c "$out_igg_bam" \
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
            --verbose 2 \
            2> "$log_file_large"
        # Run SEACR with original CTCF vs downsampled IgG
        run_seacr \
            "${BEDGRAPH_DIR}/M100K_CTCF_${ctcf_replicate}_original.bedgraph" \
            "${BEDGRAPH_DIR}/M100K_IgG_merged_${target}.bedgraph" \
            "M100K_CTCF_${ctcf_replicate}_based_on_IgG_merged_${target}"
        # Wait for files to be written
        sleep 5
        # Count number of peaks for default scale
        if [ -f "$peak_file_default" ]; then
            num_peaks_default=$(wc -l < "$peak_file_default" | awk '{print $1}')
            echo "Found ${num_peaks_default} peaks in ${peak_file_default}"
        else
            num_peaks_default=0
            echo "WARNING: Peak file not found: ${peak_file_default}"
        fi
        # Count number of peaks for large scale
        if [ -f "$peak_file_large" ]; then
            num_peaks_large=$(wc -l < "$peak_file_large" | awk '{print $1}')
            echo "Found ${num_peaks_large} peaks in ${peak_file_large}"
        else
            num_peaks_large=0
            echo "WARNING: Peak file not found: ${peak_file_large}"
        fi
        # Count number of peaks for SEACR
        if [ -f "$peak_file_seacr_stringent" ]; then
            num_peaks_seacr_stringent=$(wc -l < "$peak_file_seacr_stringent" | awk '{print $1}')
            echo "Found ${num_peaks_seacr_stringent} peaks in ${peak_file_seacr_stringent}"
        else
            num_peaks_seacr_stringent=0
            echo "WARNING: SEACR stringent file not found: ${peak_file_seacr_stringent}"
        fi
        if [ -f "$peak_file_seacr_relaxed" ]; then
            num_peaks_seacr_relaxed=$(wc -l < "$peak_file_seacr_relaxed" | awk '{print $1}')
            echo "Found ${num_peaks_seacr_relaxed} peaks in ${peak_file_seacr_relaxed}"
        else
            num_peaks_seacr_relaxed=0
            echo "WARNING: SEACR relaxed file not found: ${peak_file_seacr_relaxed}"
        fi
        # Write to summary
        echo "M100K_CTCF_${ctcf_replicate},${ctcf_replicate},original,${actual_ctcf_reads},${target},${actual_igg_reads},${num_peaks_default},${num_peaks_large},${num_peaks_seacr_stringent},${num_peaks_seacr_relaxed}" >> "$SUMMARY"
        echo "CTCF_${ctcf_replicate} (original) vs IgG downsampled to ${target}: Default: ${num_peaks_default} peaks, Large: ${num_peaks_large} peaks, SEACR Stringent: ${num_peaks_seacr_stringent} peaks, SEACR Relaxed: ${num_peaks_seacr_relaxed} peaks"
        # Clean up downsampled IgG BAM to save space (if it was created)
        if [ "$out_igg_bam" != "$IGG_BAM" ] && [ -f "$out_igg_bam" ]; then
            rm -f "$out_igg_bam"
        fi
    done
done

echo "All done! Summary written to $SUMMARY"
echo "Default scale peaks (q-value 0.01) in: $PEAK_DIR_DEFAULT"
echo "Large scale peaks (q-value 0.00001) in: $PEAK_DIR_LARGE"
echo "SEACR peaks in: $SEACR_DIR"
echo "Bedgraph files in: $BEDGRAPH_DIR"
echo "Peak files are preserved for further analysis"
echo "Original IgG merged BAM analysis included for baseline comparison"
echo "CTCF files were NOT downsampled - only IgG files were downsampled"
echo "MACs2 logs are in: $LOG_DIR"


