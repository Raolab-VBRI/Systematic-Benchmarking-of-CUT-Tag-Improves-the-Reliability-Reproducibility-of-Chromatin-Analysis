#!/bin/bash
#SBATCH --job-name=align_CUT_Tag # Job name
#SBATCH --ntasks=1 # number of tasks
#SBATCH --time 5-20:00 # time limit (day-hour:min)
#SBATCH --cpus-per-task=15 # number of threads
#SBATCH --mem-per-cpu=7Gb # requested memory
#SBATCH --account=srrao # PI's net ID
#SBATCH --mail-user=aray@mcw.edu # Email adress to send to

## adaptor trimming ##

module load cutadapt

for i in $(ls *_1.fq.gz | sed 's/\_1.fq\.gz//' |uniq)
do
    cutadapt \
        -a CTGTCTCTTATACACATCT \
        -A CTGTCTCTTATACACATCT  \
        -m 20 \
        -q 20 \
	-j 8 \
        -o Trim_cutadapt_"$i".1.fastq \
        -p Trim_cutadapt_"$i".2.fastq \
        "$i"_1.fq.gz \
        "$i"_2.fq.gz
done
