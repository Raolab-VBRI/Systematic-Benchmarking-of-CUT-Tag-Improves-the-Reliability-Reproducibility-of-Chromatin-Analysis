#!/bin/bash
#SBATCH --job-name=macs2_peaks_merged_IgG
#SBATCH --ntasks=1
#SBATCH --time=2-00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=7Gb
#SBATCH --account=srrao
#SBATCH --mail-user=aray@mcw.edu

# MACS2 Peak Calling Script for CUT&Tag with Merged IgG Controls
# This script processes mouse and rat samples with different cell numbers
# Uses --nolambda flag (common for CUT&Tag analysis)

echo "=== MACS2 Peak Calling Script (with --nolambda and merged IgG) ==="
echo "Date: $(date)"
echo

# Load required modules
module load samtools
module load macs2
module load bedtools2

# Set directory paths
TREATMENT_BAM_DIR="/scratch/g/srrao/aray/cell_no_std/dedup_bams"
CONTROL_BAM_DIR="/scratch/g/srrao/aray/cell_no_std/merged_IgG"

# Create output directories
for dir in macs2_peaks_larger_merged_IgG macs2_peaks_default_merged_IgG logs_merged_IgG comparison_merged_IgG; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory: $dir"
    fi
done

# Function for peak calling with different scaling
call_peaks_with_scaling() {
    local treatment="$1"
    local control="$2"
    local name="$3"
    local peak_type="$4"
    local scaling="$5"
    local output_dir="$6"
    local genome="$7"

    echo "Processing ${name} with ${scaling} scaling, genome: ${genome}"

    # Create log file path
    local log_file="logs_merged_IgG/${name}_${scaling}_macs2.log"

    # Set parameters based on peak type
    if [ "${peak_type}" == "broad" ]; then
        # For broad peaks (K27me3)
        if [ "${scaling}" == "default" ]; then
            macs2 callpeak \
                -t "${treatment}" \
                -c "${control}" \
                -n "${name}" \
                --broad \
                --broad-cutoff 0.01 \
                -g "${genome}" \
                -f BAMPE \
                --outdir "${output_dir}" \
                --keep-dup all \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                --nolambda \
                --verbose 2 \
                2> "${log_file}"
        else
            macs2 callpeak \
                -t "${treatment}" \
                -c "${control}" \
                -n "${name}" \
                --broad \
                --broad-cutoff 0.00001 \
                -g "${genome}" \
                -f BAMPE \
                --outdir "${output_dir}" \
                --keep-dup all \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                --scale-to large \
                --nolambda \
                --verbose 2 \
                2> "${log_file}"
        fi
    else
        # For narrow peaks (CTCF)
        if [ "${scaling}" == "default" ]; then
            macs2 callpeak \
                -t "${treatment}" \
                -c "${control}" \
                -n "${name}" \
                -g "${genome}" \
                -f BAMPE \
                --outdir "${output_dir}" \
                -q 0.01 \
                --keep-dup all \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                --nolambda \
                --verbose 2 \
                2> "${log_file}"
        else
            macs2 callpeak \
                -t "${treatment}" \
                -c "${control}" \
                -n "${name}" \
                -g "${genome}" \
                -f BAMPE \
                --outdir "${output_dir}" \
                -q 0.00001 \
                --keep-dup all \
                --nomodel \
                --shift -75 \
                --extsize 150 \
                --scale-to large \
                --nolambda \
                --verbose 2 \
                2> "${log_file}"
        fi
    fi

    # Check if peak file was created
    local peak_suffix="${peak_type}Peak"
    local peak_file="${output_dir}/${name}_peaks.${peak_suffix}"

    if [ -f "${peak_file}" ]; then
        echo "Successfully called peaks for ${name} with ${scaling} scaling"
        echo "Peak file: ${peak_file}"
        echo "Number of peaks: $(wc -l < "${peak_file}")"
    else
        echo "Warning: No peaks file generated for ${name} with ${scaling} scaling"
        echo "Check log file: ${log_file}"
    fi
}

# Initialize comparison file
echo -e "Sample\tScaling\tPeak_Type\tTotal_Peaks\tMean_Peak_Length" > comparison_merged_IgG/peak_statistics.txt

# Define parameters
cell_numbers=("100K" "50K" "10K" "1K")
species_prefixes=("M" "R")
marks=("CTCF" "K27me3")
replicates=("R1" "R2" "R3")

# Process each combination
for cell_num in "${cell_numbers[@]}"; do
    for species in "${species_prefixes[@]}"; do
        # Construct control file path
        control_file="${species}${cell_num}_IgG.merged.sorted.bam"
        control_path="${CONTROL_BAM_DIR}/${control_file}"

        # Determine genome based on species
        if [[ "$species" == "M" ]]; then
            genome="mm"
        elif [[ "$species" == "R" ]]; then
            genome="2.6e9"
        else
            echo "Error: Unknown species $species. Skipping..."
            continue
        fi

        # Check if control file exists
        if [[ ! -e "$control_path" ]]; then
            echo "Warning: Control file not found: $control_path"
            continue
        fi

        for mark in "${marks[@]}"; do
            # Determine peak type based on mark
            if [[ "${mark}" == "K27me3" ]]; then
                peak_type="broad"
            else
                peak_type="narrow"
            fi

            for replicate in "${replicates[@]}"; do
                # Construct treatment file path
                treatment_file="${species}${cell_num}_${mark}_${replicate}.dedup.bam"
                treatment_path="${TREATMENT_BAM_DIR}/${treatment_file}"

                # Sample name for output (without .dedup.bam extension)
                sample_name="${species}${cell_num}_${mark}_${replicate}"

                # Check if treatment file exists
                if [[ -e "$treatment_path" ]]; then
                    echo "=========================================="
                    echo "Processing: $treatment_file"
                    echo "Control (IgG): $control_file"
                    echo "Genome: $genome"
                    echo "Sample: $sample_name"
                    echo "Peak Type: $peak_type"
                    echo "=========================================="

                    # Call peaks with different scaling methods
                    call_peaks_with_scaling "${treatment_path}" "${control_path}" "${sample_name}" "${peak_type}" "large" "macs2_peaks_larger_merged_IgG" "${genome}"
                    call_peaks_with_scaling "${treatment_path}" "${control_path}" "${sample_name}" "${peak_type}" "default" "macs2_peaks_default_merged_IgG" "${genome}"

                    # Calculate peak statistics for each scaling method
                    for scaling in "large" "default"; do
                        if [[ "${scaling}" == "large" ]]; then
                            output_dir="macs2_peaks_larger_merged_IgG"
                        else
                            output_dir="macs2_peaks_default_merged_IgG"
                        fi

                        peak_suffix="${peak_type}Peak"
                        peak_file="${output_dir}/${sample_name}_peaks.${peak_suffix}"

                        if [ -f "${peak_file}" ]; then
                            total_peaks=$(wc -l < "${peak_file}")
                            mean_length=$(awk '{sum+=($3-$2)} END {if(NR>0) printf "%.2f", sum/NR; else print "0"}' "${peak_file}")

                            echo -e "${sample_name}\t${scaling}\t${peak_type}\t${total_peaks}\t${mean_length}" >> comparison_merged_IgG/peak_statistics.txt
                        else
                            echo "Warning: No peak file found at ${peak_file}"
                        fi
                    done

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
    echo "MACS2 Peak Calling Summary (with --nolambda and merged IgG)"
    echo "=========================================================="
    echo ""
    echo "Analysis Date: $(date)"
    echo ""
    echo "Treatment BAM Directory: ${TREATMENT_BAM_DIR}"
    echo "Control BAM Directory: ${CONTROL_BAM_DIR}"
    echo "Species: Mouse (mm) and Rat (2.6e9)"
    echo "Cell Numbers: ${cell_numbers[@]}"
    echo "Marks: CTCF (narrow), K27me3 (broad)"
    echo "Replicates: ${replicates[@]}"
    echo ""
    echo "Note: Using --nolambda flag (common for CUT&Tag analysis)"
    echo ""
    echo "Peak Statistics by Scaling Method:"
    echo "--------------------------------"
    cat comparison_merged_IgG/peak_statistics.txt
    echo ""
    echo "Analysis completed at: $(date)"
} > comparison_merged_IgG/summary_report.txt

echo "=========================================="
echo "Analysis complete!"
echo "=========================================="
echo "Results saved in:"
echo "  - Default peaks: macs2_peaks_default_merged_IgG/"
echo "  - Larger peaks: macs2_peaks_larger_merged_IgG/"
echo "  - Logs: logs_merged_IgG/"
echo "  - Summary: comparison_merged_IgG/summary_report.txt"
echo "=========================================="


