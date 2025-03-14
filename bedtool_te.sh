#!/bin/bash -l
#SBATCH --job-name=bedtools_te    			# Job name
#SBATCH --partition=fxq
#SBATCH --mail-type=END,FAIL         			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b.hunt2@exeter.ac.uk   		# Where to send mail	
#SBATCH --ntasks=1                   			# Run a single task	
#SBATCH --cpus-per-task=4            			# Number of CPU cores per task
#SBATCH --mem=10gb                    			# Job memory request
#SBATCH --output=bedtools_te.log     			# Standard output and error log
pwd; hostname; date

export OMP_NUM_THREADS=4

module load ks575/BedTools

GFF=/nobackup/beegfs/workspace/bh471/data/static_data/ns_genome/te_annotation_jg
DATA=/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/wgbs

bedtools intersect -loj -wa -wb -a $DATA/ancestral.CpG.txt -b $GFF/te.newID.gff3 > $DATA/te.CpG.annotated.txt
