#!/bin/bash
#SBATCH --job-name=align_H_samples # Job name
#SBATCH --ntasks=4 # number of parallel tasks
#SBATCH --time 5-20:00 # time limit (day-hour:min)
#SBATCH --cpus-per-task=8 # number of threads per task
#SBATCH --mem-per-cpu=7Gb # requested memory
#SBATCH --account=srrao # PI's net ID
#SBATCH --mail-user=aray@mcw.edu # Email address to send to

# Load required modules
module load bowtie/2.5.0
module load samtools
module load parallel

# Set paths
ref="/scratch/g/srrao/aray/References/GRh38_NCBI_references/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/genome"

# Create output directories
#mkdir -p aligned_SH_bams
#mkdir -p alignment_SH_logs

# Function to process one sample
process_sample() {
    local sample=$1
    local fastq1=$2
    local fastq2=$3
    
    echo "Processing sample: $sample"
    
    # Alignment with bowtie2
    echo "Aligning ${sample}..."
    bowtie2 \
        --very-sensitive-local \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -I 10 \
        -X 700 \
        --dovetail \
        --no-unal \
        -p ${SLURM_CPUS_PER_TASK} \
        -x "$ref" \
        -1 "$fastq1" \
        -2 "$fastq2" \
        -S aligned_SH_bams/"${sample}".sam \
        2> alignment_SH_logs/"${sample}"_bowtie2.log
    
    # Convert SAM to BAM and sort
    echo "Converting and sorting ${sample}..."
    samtools view -bS aligned_SH_bams/"${sample}".sam | \
        samtools sort -@ ${SLURM_CPUS_PER_TASK} -o aligned_SH_bams/"${sample}".sorted.bam -
    
    # Index the sorted BAM
    samtools index aligned_SH_bams/"${sample}".sorted.bam
    
    # Generate alignment statistics
    samtools flagstat aligned_SH_bams/"${sample}".sorted.bam > alignment_SH_logs/"${sample}"_flagstat.txt
    
    # Calculate fragment size distribution
    echo "Calculating fragment size distribution for ${sample}..."
    samtools view -f 2 aligned_SH_bams/"${sample}".sorted.bam | \
        awk '{print sqrt(($9)^2)}' > alignment_SH_logs/"${sample}"_fragment_sizes.txt
    
    # Clean up SAM file to save space
    rm -f aligned_SH_bams/"${sample}".sam
    
    echo "Completed processing for ${sample}"
}

# Export function for parallel processing
export -f process_sample
export ref
export SLURM_CPUS_PER_TASK

# Get list of H-prefixed FASTQ files
echo "Scanning for H-prefixed FASTQ files..."

# Find all _1.fastq files and extract sample names
FASTQ_FILES_1=($(ls /scratch/g/srrao/aray/cell_no_std/henikoff/fastq/H*_1.fastq 2>/dev/null))
if [ ${#FASTQ_FILES_1[@]} -eq 0 ]; then
    echo "No H*_1.fastq files found. Checking for compressed files..."
    FASTQ_FILES_1=($(ls /scratch/g/srrao/aray/cell_no_std/henikoff/fastq/H*_1.fastq.gz 2>/dev/null))
fi

if [ ${#FASTQ_FILES_1[@]} -eq 0 ]; then
    echo "ERROR: No H-prefixed FASTQ files found!"
    echo "Expected format: H0.06K_K27me3_1.fastq, H20K_K27me3_2.fastq, etc."
    exit 1
fi

echo "Found ${#FASTQ_FILES_1[@]} samples:"
for file in "${FASTQ_FILES_1[@]}"; do
    echo "  $file"
done

# Process samples in parallel
echo "Starting parallel alignment of ${#FASTQ_FILES_1[@]} samples..."

for fastq1 in "${FASTQ_FILES_1[@]}"; do
    # Extract sample name (remove _1.fastq or _1.fastq.gz)
    sample=$(basename "$fastq1" | sed 's/_1\.fastq.*$//')
    
    # Find corresponding _2.fastq file
    fastq_dir=$(dirname "$fastq1")
    if [[ "$fastq1" == *".fastq" ]]; then
        fastq2="${fastq_dir}/${sample}_2.fastq"
    else
        fastq2="${fastq_dir}/${sample}_2.fastq.gz"
    fi
    
    # Check if paired file exists
    if [ ! -f "$fastq2" ]; then
        echo "WARNING: Paired file $fastq2 not found for $fastq1. Skipping..."
        continue
    fi
    
    echo "Processing sample: $sample"
    echo "  Fastq1: $fastq1"
    echo "  Fastq2: $fastq2"
    
    # Process alignment
    process_sample "$sample" "$fastq1" "$fastq2"
done

# Wait for all parallel jobs to complete
wait

# Generate summary report
echo "Generating alignment summary report..."
echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapping_Rate(%)" > alignment_summary.txt

for bam in aligned_SH_bams/*.sorted.bam; do
    if [ -f "$bam" ]; then
        sample=$(basename "$bam" .sorted.bam)
        
        # Extract statistics from flagstat
        total_reads=$(grep "in total" alignment_SH_logs/"${sample}"_flagstat.txt | awk '{print $1}')
        mapped_reads=$(grep "mapped" alignment_SH_logs/"${sample}"_flagstat.txt | head -1 | awk '{print $1}')
        
        # Calculate mapping rate
        if [ "$total_reads" -gt 0 ]; then
            mapping_rate=$(awk -v mapped="$mapped_reads" -v total="$total_reads" 'BEGIN {printf "%.2f", (mapped/total)*100}')
        else
            mapping_rate="ERROR"
        fi
        
        echo -e "${sample}\t${total_reads}\t${mapped_reads}\t${mapping_rate}%" >> alignment_summary.txt
    fi
done

echo "Alignment complete! Check alignment_summary.txt for results."
echo "Aligned BAM files are in: aligned_SH_bams/"
echo "Log files are in: alignment_SH_logs/"
