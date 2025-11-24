#!/bin/bash

#SBATCH --job-name=seacr_peaks_merged_IgG
#SBATCH --ntasks=1
#SBATCH --time=2-00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7Gb
#SBATCH --account=srrao
#SBATCH --mail-user=aray@mcw.edu

# SEACR Peak Calling Script for Multi-species CUT&Tag with Merged IgG Controls
# This script processes mouse and rat samples with different cell numbers
# Uses existing bedgraph files (treatment and merged IgG controls)

echo "=== SEACR Peak Calling Script (with merged IgG) ==="
echo "Date: $(date)"
echo

# Load required modules
module load R
module load bedtools2

# Set directory paths
TREATMENT_BEDGRAPH_DIR="/scratch/g/srrao/aray/cell_no_std/bedgraph_files"
CONTROL_BEDGRAPH_DIR="/scratch/g/srrao/aray/cell_no_std/bedgraph_files_merged_IgG"
SEACR_PATH="/scratch/g/srrao/aray/seacr/SEACR_1.3.sh"

# Create output directories
for dir in seacr_peaks_merged_IgG logs_seacr_merged_IgG comparison_seacr_merged_IgG; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory: $dir"
    fi
done

# Function to run SEACR with multiple parameters
run_seacr() {
    local treatment_bg=$1
    local control_bg=$2
    local prefix=$3

    echo "Running SEACR analysis for ${prefix}"

    # 1. Stringent with IgG control
    echo "  Running SEACR stringent..."
    bash ${SEACR_PATH} \
        ${treatment_bg} \
        ${control_bg} \
        norm \
        stringent \
        seacr_peaks_merged_IgG/${prefix}_stringent 2> logs_seacr_merged_IgG/${prefix}_stringent.log

    # 2. Relaxed with IgG control
    echo "  Running SEACR relaxed..."
    bash ${SEACR_PATH} \
        ${treatment_bg} \
        ${control_bg} \
        norm \
        relaxed \
        seacr_peaks_merged_IgG/${prefix}_relaxed 2> logs_seacr_merged_IgG/${prefix}_relaxed.log

    echo "Completed SEACR analysis for ${prefix}"
    echo
}

# Initialize comparison file
echo -e "Sample\tStringent_Peaks\tRelaxed_Peaks" > comparison_seacr_merged_IgG/seacr_peak_summary.txt

# Define parameters
cell_numbers=("100K" "50K" "10K" "1K")
species_prefixes=("M" "R")
marks=("CTCF" "K27me3")
replicates=("R1" "R2" "R3")

# Process each combination
for cell_num in "${cell_numbers[@]}"; do
    for species in "${species_prefixes[@]}"; do
        # Construct control bedgraph file path
        control_file="${species}${cell_num}_IgG.bedgraph"
        control_path="${CONTROL_BEDGRAPH_DIR}/${control_file}"

        # Check if control file exists
        if [[ ! -e "$control_path" ]]; then
            echo "Warning: Control file not found: $control_path"
            continue
        fi

        for mark in "${marks[@]}"; do
            for replicate in "${replicates[@]}"; do
                # Construct treatment bedgraph file path
                treatment_file="${species}${cell_num}_${mark}_${replicate}.bedgraph"
                treatment_path="${TREATMENT_BEDGRAPH_DIR}/${treatment_file}"

                # Sample name for output
                sample_name="${species}${cell_num}_${mark}_${replicate}"

                # Check if treatment file exists
                if [[ -e "$treatment_path" ]]; then
                    echo "=========================================="
                    echo "Processing: $treatment_file"
                    echo "Control (IgG): $control_file"
                    echo "Sample: $sample_name"
                    echo "=========================================="

                    # Run SEACR
                    run_seacr "${treatment_path}" "${control_path}" "${sample_name}"

                    # Count peaks for each method
                    stringent_peaks=$(wc -l < "seacr_peaks_merged_IgG/${sample_name}_stringent.stringent.bed" 2>/dev/null || echo 0)
                    relaxed_peaks=$(wc -l < "seacr_peaks_merged_IgG/${sample_name}_relaxed.relaxed.bed" 2>/dev/null || echo 0)

                    echo -e "${sample_name}\t${stringent_peaks}\t${relaxed_peaks}" >> comparison_seacr_merged_IgG/seacr_peak_summary.txt

                    echo "Completed processing for $sample_name"
                    echo ""
                else
                    echo "Warning: Treatment file does not exist: $treatment_path"
                    echo ""
                fi
            done
        done
    done
done

# Generate summary report
{
    echo "SEACR Peak Calling Summary (with merged IgG)"
    echo "============================================"
    echo ""
    echo "Analysis Date: $(date)"
    echo ""
    echo "Treatment Bedgraph Directory: ${TREATMENT_BEDGRAPH_DIR}"
    echo "Control Bedgraph Directory: ${CONTROL_BEDGRAPH_DIR}"
    echo "Species: Mouse and Rat"
    echo "Cell Numbers: ${cell_numbers[@]}"
    echo "Marks: CTCF, K27me3"
    echo "Replicates: ${replicates[@]}"
    echo ""
    echo "Peak Count Summary:"
    echo "------------------"
    cat comparison_seacr_merged_IgG/seacr_peak_summary.txt
    echo ""
    echo "Analysis completed at: $(date)"
} > comparison_seacr_merged_IgG/seacr_summary_report.txt

echo "=========================================="
echo "Analysis complete!"
echo "=========================================="
echo "Results saved in:"
echo "  - SEACR peaks: seacr_peaks_merged_IgG/"
echo "  - Logs: logs_seacr_merged_IgG/"
echo "  - Summary: comparison_seacr_merged_IgG/seacr_summary_report.txt"
echo "=========================================="


