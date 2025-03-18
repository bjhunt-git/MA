#!/bin/bash -l
#SBATCH --job-name=imi_stringtie
#SBATCH --mail-type=FAIL,END 
#SBATCH --mail-user=b.hunt2@exeter.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --output=imi_stringtie.log
#SBATCH --account=c.bass
##
pwd; hostname; date
export OMP_NUM_THREADS=8

BASE=/nobackup/beegfs/workspace/bh471/data/projects/ma_imi
DATA=$BASE/rnaseq/aligned
OUT=$BASE/rnaseq/stringtie_new_annotation
SAMPLES=/nobackup/beegfs/workspace/bh471/scripts/ma_imi/rnaseq/rnaseq.txt
INDEX=/nobackup/beegfs/workspace/bh471/data/static_data/ns_genome/NS_HISAT_INDEX
GTF=/nobackup/beegfs/workspace/bh471/results/static_data/ns_annotation/ns_jg_softmask_bams/ns.gtf

if [ ! -d "${OUT}" ]; then                                                                                 
    mkdir -p "${OUT}" 
fi

while read -r j; do
    mkdir $OUT/$j
    stringtie -p 8 -e -B -A $OUT/${j}/${j}.abund.txt -G $GTF -l $j -o $OUT/${j}/${j}.gtf -l $j $DATA/${j}/${j}.sorted.bam
done < "${SAMPLES}"

