#!/bin/bash -l
#SBATCH --job-name=wgbs
#SBATCH --partition=hmq
#SBATCH --mail-type=FAIL,END 
#SBATCH --mail-user=b.hunt2@exeter.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100gb
#SBATCH --output=wgbs.log
#SBATCH --account=c.bass
##
pwd; hostname; date
export OMP_NUM_THREADS=16
BISMARK=$HOME/software/bismark_v0.22.1/

module load bowtie/2.3.4.3 samtools TrimGalore

$BISMARK/bismark_genome_preparation ns_scaffold

#per sample
while read -r j; do
	trim_galore --paired -o trim/$j --basename $j -j 2 --illumina --2colour 20 --clip_R1 10 --clip_R2 35 --fastqc $j.R1.fastq.gz $j.R2.fastq.gz
	$BISMARK/bismark ns_scaffold -1 $j_val_1.fq.gz -2 $j_val_2.fq.gz -o /aligned/$j
	$BISMARK/deduplicate_bismark -p --bam /aligned/$j_val_1_bismark_bt2_pe.bam --output_dir /aligned/$j
	$BISMARK/bismark_methylation_extractor -p --no_overlap --comprehensive --bedgraph --report --cytosine_report --genome_folder ns_scaffold --buffer_size 50% /aligned/$j/$j_val_1_bismark_bt2_pe.deduplicated.bam -o /extract/$j_CpG
	$BISMARK/bismark_methylation_extractor -p --no_overlap --comprehensive --merge_non_CpG --gzip --bedgraph --CX --report --cytosine_report --genome_folder ns_scaffold --buffer_size 50% /aligned/$j/${j}_val_1_bismark_bt2_pe.deduplicated.bam -o /extract/${j}"_all
	$BISMARK/coverage2cytosine -o /extract/$j_CpG/$j_mergeCpG --merge_CpGs --genome_folder ns_scaffold /extract/$j_CpG/$j_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
done < "${SAMPLES}"

