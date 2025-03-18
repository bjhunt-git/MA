#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=split # Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b.hunt2@exeter.ac.uk # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=8 # Number of CPU cores per task
#SBATCH --mem=50gb # Job memory request
#SBATCH --output=split.log # Standard output and error log
#SBATCH --account=c.bass


pwd; hostname; date
export OMP_NUM_THREADS=8
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
module load  Workspace slurm/18.08.4 shared samtools/1.10 TrimGalore/0.6.5 FastQC/0.11.8 HISAT2/2.1.0
DATA=/nobackup/beegfs/workspace/bh471/data/projects/ns_annotation/rnaseq
GENOME=/nobackup/beegfs/workspace/bh471/data/static_data/ns_genome
INDEX=$GENOME/NS_HISAT_INDEX
mkdir $DATA/trimmed
mkdir $DATA/aligned

trim_galore --paired -o $DATA/trimmed/ --basename NStrim -j 8 --illumina --fastqc $DATA/raw_data/NS_R1.fastq.gz $DATA/raw_data/NS_R2.fastq.gz
hisat2 -p 8 --rna-strandness RF -x ${INDEX} -1 ${DATA}/trimmed/NStrim_val_1.fq.gz -2 ${DATA}/trimmed/NStrim_val_2.fq.gz  | samtools sort > $DATA/aligned/NStrim.sorted.bam

samtools view -F 16 -o $DATA/aligned/NStrim.forward.bam $DATA/aligned/NStrim.sorted.bam
samtools view -f 16 -o $DATA/aligned/NStrim.reverse.bam $DATA/aligned/NStrim.sorted.bam