#!/bin/bash
#SBATCH --job-name=align_CUT_Tag
#SBATCH --ntasks=1
#SBATCH --time 2-00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=7Gb
#SBATCH --account=srrao
#SBATCH --mail-user=aray@mcw.edu

# Load modules
module load picard
module load samtools
module load deeptools

# Create output directory if it doesn't exist
mkdir -p /scratch/g/srrao/aray/cell_no_std/merged_bam

# Define input and output directories
INPUT_DIR="/scratch/g/srrao/aray/cell_no_std/dedup_bams"
OUTPUT_DIR="/scratch/g/srrao/aray/cell_no_std/merged_bam"

# Merge BAM files
echo "Merging BAM files..."

# IgG
java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/M100K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/M100K_IgG_R2.dedup.bam \
    I=${INPUT_DIR}/M100K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/M100K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/M50K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/M50K_IgG_R2.dedup.bam \
    I=${INPUT_DIR}/M50K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/M50K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true


java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/M10K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/M10K_IgG_R2.dedup.bam \
    I=${INPUT_DIR}/M10K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/M10K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/M1K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/M1K_IgG_R2.dedup.bam \
    I=${INPUT_DIR}/M1K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/M1K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true


# Rat IgG
java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/R100K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/R100K_IgG_R2.dedup.bam \
    I=${INPUT_DIR}/R100K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/R100K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/R50K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/R50K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/R50K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/R10K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/R10K_IgG_R2.dedup.bam \
    I=${INPUT_DIR}/R10K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/R10K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

java -jar $PICARD MergeSamFiles \
    I=${INPUT_DIR}/R1K_IgG_R1.dedup.bam \
    I=${INPUT_DIR}/R1K_IgG_R3.dedup.bam \
    O=${OUTPUT_DIR}/R1K_IgG.merged.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true


# Sort and index merged BAM files
echo "Sorting and indexing merged BAM files..."
for bam in ${OUTPUT_DIR}/*.merged.bam; do
    if [ -f "$bam" ]; then
        echo "Processing $bam"
        samtools sort "$bam" -o "${bam%.bam}.sorted.bam"
        samtools index "${bam%.bam}.sorted.bam"
    else
        echo "Warning: $bam not found"
    fi
done

# Convert to bigwig
echo "Converting to bigwig..."
for bam in ${OUTPUT_DIR}/*.merged.sorted.bam; do
    if [ -f "$bam" ]; then
        echo "Processing $bam"
        bamCoverage -b "$bam" -o "${bam%.bam}.bw" --normalizeUsing CPM
    else
        echo "Warning: $bam not found"
    fi
done

echo "Analysis complete!"

