#!/bin/bash -l
#SBATCH --job-name=bedtools_feature    			# Job name
#SBATCH --partition=fxq
#SBATCH --mail-type=END,FAIL         			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b.hunt2@exeter.ac.uk   		# Where to send mail	
#SBATCH --ntasks=1                   			# Run a single task	
#SBATCH --cpus-per-task=4            			# Number of CPU cores per task
#SBATCH --mem=10gb                    			# Job memory request
#SBATCH --output=bedtools_feature.log     			# Standard output and error log
pwd; hostname; date

export OMP_NUM_THREADS=4

module load ks575/BedTools

GFF=/nobackup/beegfs/workspace/bh471/results/static_data/ns_annotation/ns_jg_softmask_bams
DATA=/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/wgbs

#create the stripped down gff file
awk '$3 == "exon" || $3 == "intron"' $GFF/ns.jg.annotate.edit.gff3 > $GFF/ns.methstats.gff3
grep ".t1" $GFF/ns.methstats.gff3 > $GFF/ns.methstats.t1.gff3
awk '$3 == "promoter"' $GFF/ns.jg.annotate.edit.gff3 >> $GFF/ns.methstats.t1.gff3
sort -k1,1 -k4,4n -o $GFF/ns.methstats.t1.gff3 $GFF/ns.methstats.t1.gff3
sed -i 's/NA;NA;//g' $GFF/ns.methstats.t1.gff3

bedtools intersect -loj -wa -wb -a $DATA/ancestral.CpG.txt -b $GFF/ns.methstats.t1.gff3 > $DATA/ancestral.CpG.annotated.txt
