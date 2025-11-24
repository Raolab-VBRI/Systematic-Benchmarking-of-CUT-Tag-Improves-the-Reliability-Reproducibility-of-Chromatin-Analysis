#!/bin/bash
#SBATCH --job-name=computeMatrix_plotprofile_TSS_enhanced
#SBATCH --ntasks=1
#SBATCH --time=4-00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8Gb
#SBATCH --account=srrao
#SBATCH --mail-user=aray@mcw.edu

module load deeptools
module load samtools

BAM_DIR="/scratch/g/srrao/aray/cell_no_std/old_bigwig_files"
OUTPUT_DIR="/scratch/g/srrao/aray/cell_no_std/plotprofile_TSS_enhanced"

mkdir -p "${OUTPUT_DIR}/matrices" "${OUTPUT_DIR}/profiles" "${OUTPUT_DIR}/logs"

CELL_NUMBERS=("100K" "50K" "10K" "1K")
MARKS=("CTCF" "K27me3")

# Color schemes for different cell numbers
# Blue gradient for CTCF, Red gradient for K27me3
CTCF_COLORS=("#1f77b4" "#ff7f0e" "#2ca02c" "#d62728")  # Blue, Orange, Green, Red
K27ME3_COLORS=("#1f77b4" "#ff7f0e" "#2ca02c" "#d62728")  

run_plotprofile_analysis() {
    local species=$1
    local mark=$2

    echo "Processing ${species} ${mark}..."

    if [[ "$species" == "mouse" ]]; then
        ref_gtf="mm10.refGene.gtf"
    else
        ref_gtf="rn7.refGene.gtf"
    fi

    BIGWIG_FILES=()
    SAMPLE_LABELS=()
    for cell_num in "${CELL_NUMBERS[@]}"; do
        # Handle different naming conventions
        if [[ "$species" == "mouse" ]]; then
            bigwig_file="${BAM_DIR}/M${cell_num}_${mark}.merged.sorted.bw"
        else
            bigwig_file="${BAM_DIR}/R${cell_num}_${mark}.merged.sorted.bw"
        fi
        
        if [[ -f "$bigwig_file" ]]; then
            BIGWIG_FILES+=("$bigwig_file")
            SAMPLE_LABELS+=("${cell_num}")
            echo "  Found: $bigwig_file"
        else
            echo "  Missing: $bigwig_file"
        fi
    done

    if [[ ${#BIGWIG_FILES[@]} -eq 0 ]]; then
        echo "No bigWig files found for ${species} ${mark}"
        return 1
    fi

    echo "  Processing ${#BIGWIG_FILES[@]} bigWig files for ${species} ${mark}"

    # Compute matrix using reference-point mode (TSS-centered)
    computeMatrix reference-point \
        -S "${BIGWIG_FILES[@]}" \
        -R "$ref_gtf" \
        --referencePoint TSS \
        --beforeRegionStartLength 2000 \
        --afterRegionStartLength 2000 \
        --skipZeros \
        --numberOfProcessors 8 \
        -o "${OUTPUT_DIR}/matrices/${species}_${mark}_TSS_matrix.gz" \
        2> "${OUTPUT_DIR}/logs/${species}_${mark}_computematrix.log"

    if [[ $? -ne 0 ]]; then
        echo "Error: computeMatrix failed for ${species} ${mark}"
        return 1
    fi

    # Select colors based on mark
    if [[ "$mark" == "CTCF" ]]; then
        COLORS_ARRAY=("${CTCF_COLORS[@]}")
    else
        COLORS_ARRAY=("${K27ME3_COLORS[@]}")
    fi

    # Generate plotProfile with enhanced styling
    plotProfile \
        -m "${OUTPUT_DIR}/matrices/${species}_${mark}_TSS_matrix.gz" \
        -out "${OUTPUT_DIR}/profiles/${species}_${mark}_TSS_profile.pdf" \
        --perGroup \
        --colors "${COLORS_ARRAY[@]}" \
        --plotTitle "${species^} ${mark} - TSS Profile" \
        --xAxisLabel "Distance from TSS (bp)" \
        --yAxisLabel "Signal Intensity" \
        --samplesLabel "${SAMPLE_LABELS[@]}" \
        --regionsLabel "Genes" \
        --dpi 300 \
        --plotWidth 12 \
        --plotHeight 8 \
        --legendLocation "upper-right" \
        --averageType "mean" \
        --plotType "lines" \
        --lineWidth 2 \
        2> "${OUTPUT_DIR}/logs/${species}_${mark}_plotprofile.log"

    if [[ $? -ne 0 ]]; then
        echo "Error: plotProfile failed for ${species} ${mark}"
        return 1
    fi

    # Also generate PNG version
    plotProfile \
        -m "${OUTPUT_DIR}/matrices/${species}_${mark}_TSS_matrix.gz" \
        -out "${OUTPUT_DIR}/profiles/${species}_${mark}_TSS_profile.png" \
        --perGroup \
        --colors "${COLORS_ARRAY[@]}" \
        --plotTitle "${species^} ${mark} - TSS Profile" \
        --xAxisLabel "Distance from TSS (bp)" \
        --yAxisLabel "Signal Intensity" \
        --samplesLabel "${SAMPLE_LABELS[@]}" \
        --regionsLabel "Genes" \
        --dpi 300 \
        --plotWidth 12 \
        --plotHeight 8 \
        --legendLocation "upper-right" \
        --averageType "mean" \
        --plotType "lines" \
        --lineWidth 2

    # Generate a version with error bars (if you want to show variability)
    plotProfile \
        -m "${OUTPUT_DIR}/matrices/${species}_${mark}_TSS_matrix.gz" \
        -out "${OUTPUT_DIR}/profiles/${species}_${mark}_TSS_profile_with_errorbars.pdf" \
        --perGroup \
        --colors "${COLORS_ARRAY[@]}" \
        --plotTitle "${species^} ${mark} - TSS Profile (with SEM)" \
        --xAxisLabel "Distance from TSS (bp)" \
        --yAxisLabel "Signal Intensity" \
        --samplesLabel "${SAMPLE_LABELS[@]}" \
        --regionsLabel "Genes" \
        --dpi 300 \
        --plotWidth 12 \
        --plotHeight 8 \
        --legendLocation "upper-right" \
        --averageType "mean" \
        --plotType "lines" \
        --lineWidth 2 \
        --errorBars

    echo "  Completed ${species} ${mark}"
    echo "    - Main profile: ${OUTPUT_DIR}/profiles/${species}_${mark}_TSS_profile.pdf"
    echo "    - With error bars: ${OUTPUT_DIR}/profiles/${species}_${mark}_TSS_profile_with_errorbars.pdf"
}

# Function to create summary statistics
create_summary_stats() {
    echo "Creating summary statistics..."
    
    SUMMARY_FILE="${OUTPUT_DIR}/analysis_summary.txt"
    echo "PlotProfile Analysis Summary" > "$SUMMARY_FILE"
    echo "Generated on: $(date)" >> "$SUMMARY_FILE"
    echo "=================================" >> "$SUMMARY_FILE"
    echo >> "$SUMMARY_FILE"
    
    for mark in "${MARKS[@]}"; do
        echo "=== ${mark} Analysis ===" >> "$SUMMARY_FILE"
        for species in "mouse" "rat"; do
            matrix_file="${OUTPUT_DIR}/matrices/${species}_${mark}_TSS_matrix.gz"
            if [[ -f "$matrix_file" ]]; then
                echo "${species^} ${mark}: Matrix generated successfully" >> "$SUMMARY_FILE"
            else
                echo "${species^} ${mark}: Matrix generation failed" >> "$SUMMARY_FILE"
            fi
        done
        echo >> "$SUMMARY_FILE"
    done
    
    echo "Summary saved to: $SUMMARY_FILE"
}

main() {
    echo "Starting enhanced plotProfile analysis with TSS reference point..."
    echo "Output directory: ${OUTPUT_DIR}"
    echo "Cell numbers: ${CELL_NUMBERS[*]}"
    echo "Marks: ${MARKS[*]}"
    echo

    for mark in "${MARKS[@]}"; do
        echo "=== Processing ${mark} ==="
        run_plotprofile_analysis "mouse" "$mark"
        run_plotprofile_analysis "rat" "$mark"
        echo
    done

    create_summary_stats

    echo "=== Analysis Complete ==="
    echo "Generated plots:"
    echo "1. Mouse CTCF TSS Profile (with and without error bars)"
    echo "2. Mouse K27me3 TSS Profile (with and without error bars)" 
    echo "3. Rat CTCF TSS Profile (with and without error bars)"
    echo "4. Rat K27me3 TSS Profile (with and without error bars)"
    echo
    echo "All plots saved in: ${OUTPUT_DIR}/profiles/"
    echo "Matrices saved in: ${OUTPUT_DIR}/matrices/"
    echo "Logs saved in: ${OUTPUT_DIR}/logs/"
    echo "Summary saved in: ${OUTPUT_DIR}/analysis_summary.txt"
}

main "$@"

