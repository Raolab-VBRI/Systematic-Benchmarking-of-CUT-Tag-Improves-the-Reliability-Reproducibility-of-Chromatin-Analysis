#!/bin/bash
#SBATCH --job-name=dedup # Job name
#SBATCH --ntasks=4 # number of parallel tasks
#SBATCH --time=12:00:00 # time limit (hour:min:sec)
#SBATCH --cpus-per-task=8 # number of threads per task
#SBATCH --mem-per-cpu=7Gb # requested memory
#SBATCH --account=srrao # PI's net ID
#SBATCH --mail-user=aray@mcw.edu # Email address to send to

# Load required modules
module load samtools
module load parallel
module load picard

picardCMD="java -jar $PICARD"

# Create output directories
mkdir -p dedup_bams
mkdir -p dedup_metrics

# Function to process deduplication for one sample
process_dedup() {
    local i=$1
    
    echo "Removing duplicates for ${i}..."
    
    # Remove duplicates with Picard
    $picardCMD MarkDuplicates \
        I=aligned_bams/"${i}".sorted.bam \
        O=dedup_bams/"${i}".dedup.bam \
        M=dedup_metrics/"${i}".dedup.metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT

    # Index the deduplicated BAM
    samtools index dedup_bams/"${i}".dedup.bam

    # Generate post-deduplication statistics
    samtools flagstat dedup_bams/"${i}".dedup.bam > alignment_logs/"${i}"_dedup_flagstat.txt

    # Calculate fragment size distribution
    echo "Calculating fragment size distribution for ${i}..."
    samtools view -f 2 dedup_bams/"${i}".dedup.bam | \
        awk '{print sqrt(($9)^2)}' > alignment_logs/"${i}"_fragment_sizes.txt

    echo "Completed deduplication for ${i}"
}

# Export function and variables for parallel
export -f process_dedup
export picardCMD
export SLURM_CPUS_PER_TASK

# Get list of samples from sorted BAM files
SAMPLES=($(ls aligned_bams/*.sorted.bam | sed 's/.*\///' | sed 's/\.sorted\.bam//' | uniq))

# Process samples in parallel
echo "Starting parallel deduplication of ${#SAMPLES[@]} samples..."
parallel --jobs ${SLURM_NTASKS} \
    --joblog parallel_dedup.log \
    process_dedup ::: "${SAMPLES[@]}"

# Wait for all parallel jobs to complete
wait

# Generate summary report with CORRECTED paired-end parsing
echo "Generating summary report..."
echo -e "Sample\tTotal_Reads\tMapped_Reads\tDuplicated_Reads\tDuplication_Rate(%)" > dedup_summary.txt

for i in dedup_metrics/*.dedup.metrics.txt
do
    sample=$(basename $i .dedup.metrics.txt)
    
    echo "Processing metrics for sample: $sample"
    
    # Get the actual data line (not the header)
    metrics_line=$(grep -v "^#" "$i" | grep -v "^$" | grep -v "^LIBRARY" | head -n 1)
    
    if [ -n "$metrics_line" ]; then
        echo "Metrics line found: $metrics_line"
        
        # CORRECTED: Extract the right columns for paired-end data
        # Column 2: UNPAIRED_READS_EXAMINED (single-end reads)
        # Column 3: READ_PAIRS_EXAMINED (paired-end read pairs)
        # Column 6: UNPAIRED_READ_DUPLICATES (single-end duplicates)
        # Column 7: READ_PAIR_DUPLICATES (paired-end duplicates)
        
        unpaired_reads=$(echo "$metrics_line" | cut -f 2)
        read_pairs=$(echo "$metrics_line" | cut -f 3)
        unpaired_duplicates=$(echo "$metrics_line" | cut -f 6)
        read_pair_duplicates=$(echo "$metrics_line" | cut -f 7)
        
        # Calculate total reads: (unpaired_reads + read_pairs * 2)
        total_reads=$((unpaired_reads + read_pairs * 2))
        
        # Calculate total duplicates: (unpaired_duplicates + read_pair_duplicates * 2)
        total_duplicates=$((unpaired_duplicates + read_pair_duplicates * 2))
        
        echo "Unpaired reads: $unpaired_reads"
        echo "Read pairs: $read_pairs"
        echo "Total reads: $total_reads"
        echo "Unpaired duplicates: $unpaired_duplicates"
        echo "Read pair duplicates: $read_pair_duplicates"
        echo "Total duplicates: $total_duplicates"
        
        # Calculate duplication rate as percentage
        if [ "$total_reads" -gt 0 ] && [ "$total_reads" -ge "$total_duplicates" ]; then
            duplication_rate=$(awk -v dup="$total_duplicates" -v total="$total_reads" 'BEGIN {printf "%.2f", (dup/total)*100}')
        else
            echo "WARNING: Invalid read counts for $sample - Total: $total_reads, Duplicated: $total_duplicates"
            duplication_rate="ERROR"
        fi
        
        # Get mapped reads from deduplicated BAM
        mapped_reads=$(samtools view -c -F 4 dedup_bams/"${sample}".dedup.bam)
        
        echo "Mapped reads: $mapped_reads"
        echo "Duplication rate: $duplication_rate%"
        echo "---"
        
        echo -e "${sample}\t${total_reads}\t${mapped_reads}\t${total_duplicates}\t${duplication_rate}%" >> dedup_summary.txt
    else
        echo "ERROR: Could not parse metrics for sample ${sample}"
        echo -e "${sample}\tERROR\tERROR\tERROR\tERROR" >> dedup_summary.txt
    fi
done

echo "Deduplication complete! Check dedup_summary.txt for results."
