#!/bin/bash -l
#SBATCH --job-name=ns_jg_softmask_bams # Job name
#SBATCH --partition=hmq
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b.hunt2@exeter.ac.uk # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=16 # Number of CPU cores per task
#SBATCH --mem=100gb # Job memory request
#SBATCH --output=ns_jg_softmask_bams.log # Standard output and error log
#SBATCH --account=c.bass


pwd; hostname; date
export OMP_NUM_THREADS=32
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
module load  Workspace slurm/18.08.4 shared Conda/Python3 GeneMark-ES braker ActivePerl
DATA=/nobackup/beegfs/workspace/bh471/data/projects/ns_annotation/rnaseq
GENOME=/nobackup/beegfs/workspace/bh471/data/static_data/ns_genome
OUTPUT=/nobackup/beegfs/workspace/bh471/results/static_data/ns_annotation
cd $OUTPUT

wd=ns_jg_softmask_bams
export AUGUSTUS_SCRIPTS_PATH=/cm/shared/admin-apps/braker/2.1.6/bin
export PROTHINT_PATH=/cm/shared/admin-apps/GeneMark-ES/4.68_lic/ProtHint/bin
export PERL_PATH=/cm/shared/admin-apps/braker/2.1.6/bin
if [ -d $wd ]; then
    rm -r $wd
fi
#--useexisting
#UTR addition is apparently unstable https://github.com/Gaius-Augustus/BRAKER/issues/419

braker.pl --cores 16 --gff3 --species=myzus_persicae_jg_soft --workingdir=$wd --softmasking --genome=$GENOME/ns.fasta --bam=$DATA/aligned/NStrim.forward.bam,$DATA/aligned/NStrim.reverse.bam --stranded=+,- 
