#!/bin/bash
#SBATCH --job-name=fripscore
#SBATCH --ntasks=1
#SBATCH --time=2-20:00
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=7Gb
#SBATCH --account=srrao
#SBATCH --mail-user=aray@mcw.edu

# Load required modules
module load bedtools2
module load R
module load samtools
module load python
module load gcc/9.3.0

# Check if bc is available for floating-point calculations
if ! command -v bc &> /dev/null; then
    echo "ERROR: 'bc' command not found. This is required for floating-point calculations."
    echo "Please install bc or load the appropriate module."
    exit 1
fi

# Set error handling
set -e
set -o pipefail

# Define paths
BAM_DIR="/scratch/g/srrao/aray/cell_no_std/dedup_bams"
PEAKS_LARGE_DIR="/scratch/g/srrao/aray/cell_no_std/downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF/macs2_peaks_large"
PEAKS_DEFAULT_DIR="/scratch/g/srrao/aray/cell_no_std/downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF/macs2_peaks_default"
SEACR_DIR="/scratch/g/srrao/aray/cell_no_std/downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF/seacr_peaks"
OUTPUT_DIR="./frip_results_downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF"
LOG_DIR="./frip_logs_downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Define downsampling levels (extracted from peak file names)
DOWNSAMPLE_LEVELS=("original" "5000000" "2500000" "1000000" "500000" "100000" "50000" "10000")
REPLICATES=("R1" "R2" "R3")

# Initialize output CSV file
OUTPUT_CSV="$OUTPUT_DIR/frip_scores_summary.csv"
echo "Sample_Name,IgG_Downsampled_Reads,Replicate,Total_Reads,Peaks_in_Regions,FRiP_macs_default,FRiP_macs_large,FRiP_seacr_relaxed,FRiP_seacr_stringent,Peak_Caller_Used,Notes" > "$OUTPUT_CSV"

# Function to calculate FRiP score
calculate_frip() {
    local bam_file="$1"
    local peak_file="$2"
    local output_file="$3"
    
    # Check if files exist
    if [[ ! -f "$bam_file" ]]; then
        echo "ERROR: BAM file not found: $bam_file" >&2
        return 1
    fi
    
    if [[ ! -f "$peak_file" ]]; then
        echo "ERROR: Peak file not found: $peak_file" >&2
        return 1
    fi
    
    # Get total reads in BAM
    local total_reads=$(samtools view -c "$bam_file")
    
    # Count reads in peaks
    local reads_in_peaks=0
    
    # Handle different peak file formats
    if [[ "$peak_file" == *.narrowPeak ]]; then
        # MACS2 narrowPeak format - use -bed flag for BAM input
        reads_in_peaks=$(bedtools intersect -a "$bam_file" -b "$peak_file" -c -bed | awk '{sum += $NF} END {print sum+0}')
    elif [[ "$peak_file" == *.bed ]]; then
        # SEACR BED format - use -bed flag for BAM input
        reads_in_peaks=$(bedtools intersect -a "$bam_file" -b "$peak_file" -c -bed | awk '{sum += $NF} END {print sum+0}')
    else
        echo "ERROR: Unknown peak file format: $peak_file" >&2
        return 1
    fi
    
    # Ensure reads_in_peaks is a valid number
    if [[ ! "$reads_in_peaks" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Invalid reads_in_peaks value: $reads_in_peaks" >&2
        reads_in_peaks=0
    fi
    
    # Calculate FRiP score
    local frip_score=$(echo "scale=6; $reads_in_peaks / $total_reads" | bc -l)
    
    # Write results to output file
    echo "$total_reads,$reads_in_peaks,$frip_score" > "$output_file"
    
    echo "FRiP calculation completed: $peak_file"
    echo "  Total reads: $total_reads"
    echo "  Reads in peaks: $reads_in_peaks"
    echo "  FRiP score: $frip_score"
}

# Function to process a single sample
process_sample() {
    local replicate="$1"
    local downsample_level="$2"
    
    local sample_name="M100K_CTCF_$replicate"
    local bam_file="$BAM_DIR/${sample_name}.dedup.bam"
    
    echo "Processing $sample_name with $downsample_level downsampling..."
    
    # Initialize variables
    local frip_macs_default="NA"
    local frip_macs_large="NA"
    local frip_seacr_relaxed="NA"
    local frip_seacr_stringent="NA"
    local total_reads="NA"
    local peaks_in_regions="NA"
    local peak_caller_used=""
    local notes=""
    
    # Check if BAM file exists
    if [[ ! -f "$bam_file" ]]; then
        echo "WARNING: BAM file not found: $bam_file"
        notes="BAM file missing"
    else
        # Get total reads
        total_reads=$(samtools view -c "$bam_file")
        
        # Process MACS2 default peaks
        local macs_default_peak="$PEAKS_DEFAULT_DIR/${sample_name}_based_on_IgG_merged_${downsample_level}_default_peaks.narrowPeak"
        if [[ -f "$macs_default_peak" ]]; then
            local temp_output="$OUTPUT_DIR/temp_macs_default_${replicate}_${downsample_level}.txt"
            if calculate_frip "$bam_file" "$macs_default_peak" "$temp_output"; then
                local macs_default_data=$(cat "$temp_output")
                frip_macs_default=$(echo "$macs_default_data" | cut -d',' -f3)
                peaks_in_regions=$(echo "$macs_default_data" | cut -d',' -f2)
                rm "$temp_output"
            fi
        fi
        
        # Process MACS2 large peaks
        local macs_large_peak="$PEAKS_LARGE_DIR/${sample_name}_based_on_IgG_merged_${downsample_level}_large_peaks.narrowPeak"
        if [[ -f "$macs_large_peak" ]]; then
            local temp_output="$OUTPUT_DIR/temp_macs_large_${replicate}_${downsample_level}.txt"
            if calculate_frip "$bam_file" "$macs_large_peak" "$temp_output"; then
                local macs_large_data=$(cat "$temp_output")
                frip_macs_large=$(echo "$macs_large_data" | cut -d',' -f3)
                rm "$temp_output"
            fi
        fi
        
        # Process SEACR relaxed peaks
        local seacr_relaxed_peak="$SEACR_DIR/${sample_name}_based_on_IgG_merged_${downsample_level}_relaxed.relaxed.bed"
        if [[ -f "$seacr_relaxed_peak" ]]; then
            local temp_output="$OUTPUT_DIR/temp_seacr_relaxed_${replicate}_${downsample_level}.txt"
            if calculate_frip "$bam_file" "$seacr_relaxed_peak" "$temp_output"; then
                local seacr_relaxed_data=$(cat "$temp_output")
                frip_seacr_relaxed=$(echo "$seacr_relaxed_data" | cut -d',' -f3)
                rm "$temp_output"
            fi
        fi
        
        # Process SEACR stringent peaks
        local seacr_stringent_peak="$SEACR_DIR/${sample_name}_based_on_IgG_merged_${downsample_level}_stringent.stringent.bed"
        if [[ -f "$seacr_stringent_peak" ]]; then
            local temp_output="$OUTPUT_DIR/temp_seacr_stringent_${replicate}_${downsample_level}.txt"
            if calculate_frip "$bam_file" "$seacr_stringent_peak" "$temp_output"; then
                local seacr_stringent_data=$(cat "$temp_output")
                frip_seacr_stringent=$(echo "$seacr_stringent_data" | cut -d',' -f3)
                rm "$temp_output"
            fi
        fi
        
        # Determine which peak caller was used
        local peak_callers=()
        [[ "$frip_macs_default" != "NA" ]] && peak_callers+=("macs_default")
        [[ "$frip_macs_large" != "NA" ]] && peak_callers+=("macs_large")
        [[ "$frip_seacr_relaxed" != "NA" ]] && peak_callers+=("seacr_relaxed")
        [[ "$frip_seacr_stringent" != "NA" ]] && peak_callers+=("seacr_stringent")
        peak_caller_used=$(IFS=";"; echo "${peak_callers[*]}")
        
        if [[ -z "$peak_caller_used" ]]; then
            notes="No peak files found for this downsampling level"
        fi
    fi
    
    # Write to CSV
    echo "$sample_name,$downsample_level,$replicate,$total_reads,$peaks_in_regions,$frip_macs_default,$frip_macs_large,$frip_seacr_relaxed,$frip_seacr_stringent,$peak_caller_used,$notes" >> "$OUTPUT_CSV"
    
    echo "Completed processing $sample_name with $downsample_level downsampling"
}

# Main processing loop
echo "Starting FRiP score calculation..."
echo "Output will be saved to: $OUTPUT_CSV"
echo "Logs will be saved to: $LOG_DIR"

# Process each replicate and downsampling level
for replicate in "${REPLICATES[@]}"; do
    for downsample_level in "${DOWNSAMPLE_LEVELS[@]}"; do
        echo "=========================================="
        echo "Processing $replicate with $downsample_level downsampling"
        echo "=========================================="
        
        # Log start time
        start_time=$(date)
        echo "Started at: $start_time"
        
        # Process the sample
        process_sample "$replicate" "$downsample_level"
        
        # Log completion time
        end_time=$(date)
        echo "Completed at: $end_time"
        echo ""
    done
done

# Generate summary statistics
echo "Generating summary statistics..."
Rscript - << 'EOF'
library(dplyr)
library(ggplot2)

# Read the CSV file
data <- read.csv("frip_results_downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF/frip_scores_summary.csv", stringsAsFactors = FALSE)

# Convert FRiP scores to numeric, handling "NA" values
numeric_cols <- c("FRiP_macs_default", "FRiP_macs_large", "FRiP_seacr_relaxed", "FRiP_seacr_stringent")
for(col in numeric_cols) {
  data[[col]] <- as.numeric(data[[col]])
}

# Summary statistics
summary_stats <- data %>%
  group_by(Replicate, IgG_Downsampled_Reads) %>%
  summarise(
    Mean_FRiP_macs_default = mean(FRiP_macs_default, na.rm = TRUE),
    Mean_FRiP_macs_large = mean(FRiP_macs_large, na.rm = TRUE),
    Mean_FRiP_seacr_relaxed = mean(FRiP_seacr_relaxed, na.rm = TRUE),
    Mean_FRiP_seacr_stringent = mean(FRiP_seacr_stringent, na.rm = TRUE),
    .groups = 'drop'
  )

# Write summary statistics
write.csv(summary_stats, "frip_results_downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF/frip_summary_statistics.csv", row.names = FALSE)

# Create FRiP score comparison plot
pdf("frip_results_downsampled_peak_results_for_merge_M100K_IgG_no_downsample_M100K_CTCF/frip_scores_comparison.pdf", width = 12, height = 8)

# Reshape data for plotting
plot_data <- data %>%
  select(Sample_Name, IgG_Downsampled_Reads, Replicate, 
         FRiP_macs_default, FRiP_macs_large, FRiP_seacr_relaxed, FRiP_seacr_stringent) %>%
  tidyr::pivot_longer(cols = starts_with("FRiP_"), 
                      names_to = "Peak_Caller", 
                      values_to = "FRiP_Score") %>%
  filter(!is.na(FRiP_Score))

# Create plot
p <- ggplot(plot_data, aes(x = IgG_Downsampled_Reads, y = FRiP_Score, color = Peak_Caller, group = Peak_Caller)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~Replicate, scales = "free_y") +
  theme_minimal() +
  labs(title = "FRiP Scores by Peak Caller and Downsampling Level",
       x = "IgG Downsampled Reads",
       y = "FRiP Score",
       color = "Peak Caller") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()

cat("Summary statistics and plots generated successfully!\n")
EOF

echo "FRiP score calculation completed!"
echo "Results saved to: $OUTPUT_CSV"
echo "Summary statistics saved to: $OUTPUT_DIR/frip_summary_statistics.csv"
echo "Comparison plot saved to: $OUTPUT_DIR/frip_scores_comparison.pdf"
echo "All temporary files cleaned up."

# Final cleanup
rm -f "$OUTPUT_DIR"/temp_*.txt

echo "Script execution completed successfully!" 
