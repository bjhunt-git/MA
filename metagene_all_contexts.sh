#!/bin/bash

#open the environment in the command line first
#source ~/miniconda3/bin/activate
#conda activate deeptools

sort -k1,1 -k2,2n -k3,3n ancestral.bedGraph > ancestral.sorted.bedGraph
grep -w -E "HiC_scaffold_1|HiC_scaffold_2|HiC_scaffold_3|HiC_scaffold_4|HiC_scaffold_5|HiC_scaffold_6" cg.bedGraph > cg_six.bedGraph
./bedGraphToBigWig cg_six.bedGraph chrom_sizes.txt cg.bw
grep -w -E "HiC_scaffold_1|HiC_scaffold_2|HiC_scaffold_3|HiC_scaffold_4|HiC_scaffold_5|HiC_scaffold_6" chg.bedGraph > chg_six.bedGraph
./bedGraphToBigWig chg_six.bedGraph chrom_sizes.txt chg.bw
grep -w -E "HiC_scaffold_1|HiC_scaffold_2|HiC_scaffold_3|HiC_scaffold_4|HiC_scaffold_5|HiC_scaffold_6" chh.bedGraph > chh_six.bedGraph
./bedGraphToBigWig chh_six.bedGraph chrom_sizes.txt chh.bw

computeMatrix scale-regions -S cg.bw chg.bw chh.bw -R NS_six.gtf -b 1000 -a 1000 --metagene -bs 50 --regionBodyLength 1000 -o all_contexts_matrix.gz
plotProfile -m all_contexts_matrix.gz --perGroup -o all_contexts.png --outFileNameData all_contexts.tab --colors blue red green --plotHeight 8 --plotWidth 8 --legendLocation best --samplesLabel CG CHG CHH

