#!/bin/bash -l
#SBATCH --job-name=bedtools_union    			# Job name
#SBATCH --partition=fxq
#SBATCH --mail-type=END,FAIL         			# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b.hunt2@exeter.ac.uk   		# Where to send mail	
#SBATCH --ntasks=1                   			# Run a single task	
#SBATCH --cpus-per-task=8            			# Number of CPU cores per task
#SBATCH --mem=50gb                    			# Job memory request
#SBATCH --output=bedtools_union.log     			# Standard output and error log
pwd; hostname; date

export OMP_NUM_THREADS=8

module load ks575/BedTools
WD=/nobackup/beegfs/workspace/bh471/data/projects/ma_imi/wgbs/extract/bedgraph/

cd $WD

gunzip *.bedGraph.gz
for i in *bedGraph; do
  echo "sorting" $i 
  bedtools sort -i $i > $i.sorted.bedGraph
done

bedtools unionbedg -i D01_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D02_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D03_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D04_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D05_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D06_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D07_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D08_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D09_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D10_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D11_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D12_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D13_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D14_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D15_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D16_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D17_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D18_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D19_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D20_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D21_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D22_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D23_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D24_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D25_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D26_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D27_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D28_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D29_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D30_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D31_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D32_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D33_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D34_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D35_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D36_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D37_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D38_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D39_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D40_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D41_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D42_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D43_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D44_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D45_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D46_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D47_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D48_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D49_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D50_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D51_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D52_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D53_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D54_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D55_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D56_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D57_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D58_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D59_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph D60_val_1_bismark_bt2_pe.deduplicated.bedGraph.sorted.bedGraph -filler NA -header -names D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 D14 D15 D16 D17 D18 D19 D20 D21 D22 D23 D24 D25 D26 D27 D28 D29 D30 D31 D32 D33 D34 D35 D36 D37 D38 D39 D40 D41 D42 D43 D44 D45 D46 D47 D48 D49 D50 D51 D52 D53 D54 D55 D56 D57 D58 D59 D60 > union.bedGraph

grep -v "NA" union.bedGraph > union.na.omit.bedGraph