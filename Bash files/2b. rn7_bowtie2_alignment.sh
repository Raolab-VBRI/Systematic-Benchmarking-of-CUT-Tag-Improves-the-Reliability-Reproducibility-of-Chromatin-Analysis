#!/bin/bash
#SBATCH --job-name=align_dedup # Job name
#SBATCH --ntasks=4 # number of parallel tasks
#SBATCH --time 5-20:00 # time limit (day-hour:min)
#SBATCH --cpus-per-task=8 # number of threads per task
#SBATCH --mem-per-cpu=7Gb # requested memory
#SBATCH --account=srrao # PI's net ID
#SBATCH --mail-user=aray@mcw.edu # Email address to send to

# Load required modules
module load bowtie/2.5.0
module load samtools
module load picard
module load parallel

picardCMD="java -jar $PICARD"

# Set paths
ref="/scratch/g/srrao/aray/AR_CUT_TAG/Rat/rn7"

# Create output directories
#mkdir -p aligned_bams
#mkdir -p alignment_logs


# Function to process one sample
process_sample() {
    local i=$1
    
    echo "Processing sample: $i"
    
    # Alignment with bowtie2
    echo "Aligning ${i}..."
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
        -x ${ref} \
        -1 Trim_cutadapt_"${i}".1.fastq \
        -2 Trim_cutadapt_"${i}".2.fastq \
        2> alignment_logs/"${i}"_bowtie2.log \
        | samtools view -bS - \
        | samtools sort -@ ${SLURM_CPUS_PER_TASK} -o aligned_bams/"${i}".sorted.bam 

    # Index the sorted BAM
    samtools index aligned_bams/"${i}".sorted.bam 

    #Generate alignment statistics
    samtools flagstat aligned_bams/"${i}".sorted.bam > alignment_logs/"${i}"_flagstat.txt

}

# Export function and variables for parallel
export -f process_sample
export ref
export SLURM_CPUS_PER_TASK

# Get list of samples
SAMPLES=($(ls Trim_cutadapt_R*.1.fastq | sed 's/Trim_cutadapt_//' | sed 's/\.1\.fastq//' | uniq))

# Process samples in parallel
echo "Starting parallel processing of ${#SAMPLES[@]} samples..."
parallel --jobs ${SLURM_NTASKS} \
    --joblog parallel_job.log \
    process_sample ::: "${SAMPLES[@]}"

# Wait for all parallel jobs to complete
wait


echo "All processing complete!"

# Print completion time and parallel job log
echo "Parallel job log summary:"
cat parallel_job.log
