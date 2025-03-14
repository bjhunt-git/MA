#!/bin/bash -l
#SBATCH --job-name=wgs
#SBATCH --partition=hmq
#SBATCH --mail-type=FAIL,END 
#SBATCH --mail-user=b.hunt2@exeter.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100gb
#SBATCH --output=wgs.log
#SBATCH --account=c.bass
##
pwd; hostname; date
export OMP_NUM_THREADS=16
module load samtools bwa/0.7.17 Picard/2.18.26 Conda/Python2/2.7.15 GATK/4.1.0.0 TrimGalore FastQC  

samtools faidx ns.fasta
picard CreateSequenceDictionary R=ns.fasta O=ns.dict
bwa index -p ns ns.fasta

#per sample
while read -r j; do
	trim_galore --paired -o trim/$j/ --basename ${j} -j 4 --illumina --2colour 20 $j.R1.fastq.gz $j.R2.fastq.gz
	SM=$j 
	PL="ILLUMINA"
	RGID=$i
	PU="unit1"	
	bwa mem -M -t 8 -R \"@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$SM\tPU:$PU\" ns trim/${j}/${j}_val_1.fq.gz trim/${j}/${j}_val_2.fq.gz | samtools view -b -S -F 4 -o aligned/$j/$j.bam 
	samtools sort aligned/$j/$j.bam -o aligned/$j/$j.sorted.bam -@ 8
	picard MarkDuplicates I=aligned/$j/$j.sorted.bam O=aligned/$j/$j.sorted.dup.bam M=aligned/$j/$j.marked_dup_metrics.txt TMP_DIR=$TMP
	samtools index aligned/$j/$j.sorted.dup.bam -@ 7
	gatk HaplotypeCaller -R ns.fasta -I aligned/"${j}"/"${j}".sorted.dup.bam --native-pair-hmm-threads 8 --tmp-dir tmp --output-mode EMIT_ALL_CONFIDENT_SITES --emit-ref-confidence GVCF -O vcf/"${j}".g.vcf -bamout aligned/"${j}"/"${j}".bamout.bam
done < "${SAMPLES}"

gatk CombineGVCFs -R ns.fasta -V vcf/A10237_1/A10237_1.g.vcf -V vcf/A10238_1/A10238_1.g.vcf -V vcf/A10239_1/A10239_1.g.vcf -V vcf/A10240_1/A10240_1.g.vcf -V vcf/A10241_1/A10241_1.g.vcf -V vcf/A10242_1/A10242_1.g.vcf -V vcf/A10243_1/A10243_1.g.vcf -V vcf/A10244_1/A10244_1.g.vcf -V vcf/A10245_1/A10245_1.g.vcf -V vcf/A10246_1/A10246_1.g.vcf -V vcf/A10247_1/A10247_1.g.vcf -V vcf/A10248_1/A10248_1.g.vcf -V vcf/A10249_1/A10249_1.g.vcf -V vcf/A10250_1/A10250_1.g.vcf -V vcf/A10251_1/A10251_1.g.vcf -V vcf/A10252_1/A10252_1.g.vcf -V vcf/A10253_1/A10253_1.g.vcf -V vcf/A10254_1/A10254_1.g.vcf -O vcf/all.g.vcf

gatk GenotypeGVCFs -R ns.fasta -V vcf/all.g.vcf --include-non-variant-sites -O vcf/genotyped.vcf

#exposed lines
OUT=/vcf/exposed
mkdir -p $OUT
gatk SelectVariants -R ns.fasta -V vcf/genotyped.vcf -O $OUT/exposed.vcf \
  -sn A10245_1 \
  -sn A10246_1 \
  -sn A10247_1 \
  -sn A10248_1 \
  -sn A10249_1 \
  -sn A10250_1 \
  -sn A10254_1 \
  --selectExpressions 'vc.getGenotype("A10254_1").isHomRef()'

gatk SelectVariants -R ns.fasta -V $OUT/exposed.vcf -O $OUT/snps.vcf \
  --select-type-to-include SNP

#filter at INFO level and further filter on FORMAT level for genotype quality and read depth
gatk VariantFiltration -R ns.fasta -V $OUT/snps.vcf -O $OUT/filtered.snps.vcf \
	  --filter-name "QD5" \
	  --filter-expression "QD < 5.0" \
	  --filter-name "FS60" \
	  --filter-expression "FS > 60.0" \
    --filter-name "SOR3" \
	  --filter-expression "SOR > 3.0" \
	  --filter-name "MQ50" \
	  --filter-expression "MQ < 55.0" \
	  --filter-name "MQRS-12.5" \
	  --filter-expression "MQRankSum < -12.5" \
  	--filter-name "RPRS-8" \
	  --filter-expression "ReadPosRankSum < -8.0" 
     
#filter at format level
gatk VariantFiltration -R ns.fasta -V $OUT/filtered.snps.vcf -O $OUT/gt.filtered.snps.vcf \
	  --genotype-filter-name "lowGQ" \
    --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "lowDP" \
    --genotype-filter-expression "DP < 10" 
    
#convert FORMAT filters to no call because downstream tools don't parse the FORMAT level filter
gatk SelectVariants -R ns.fasta -V $OUT/gt.filtered.snps.vcf --set-filtered-gt-to-nocall -O $OUT/nocall.gt.filtered.snps.vcf

#remove sites that have a no call in one sample or more
gatk SelectVariants -R ns.fasta -V $OUT/nocall.gt.filtered.snps.vcf --max-nocall-number 0 --exclude-filtered TRUE -O $OUT/snps.final.vcf
bcftools stats $OUT/snps.final.vcf > $OUT/snps.final.vcf.stats.txt 
#do the same for all variants
gatk VariantFiltration -R ns.fasta -V $OUT/exposed.vcf -O $OUT/exposed.gt.filtered.vcf \
    --genotype-filter-name "lowDP" \
    --genotype-filter-expression "DP < 10" 
    
gatk SelectVariants -R ns.fasta -V $OUT/exposed.gt.filtered.vcf --set-filtered-gt-to-nocall -O $OUT/exposed.nocall.gt.filtered.vcf

gatk SelectVariants -R ns.fasta -V $OUT/exposed.nocall.gt.filtered.vcf --max-nocall-number 0 --exclude-filtered TRUE -O $OUT/exposed.final.vcf

bcftools stats $OUT/unexposed.final.vcf > $OUT/unexposed.final.vcf.stats.txt

cd $OUT
mkdir unique_vcfs
while read -r i; do
  bcftools view --private --samples $i snps.final.vcf > unique_vcfs/"${i}".unique.vcf
  bcftools query -f '%CHROM\t%POS[\t%GT\t]\n' unique_vcfs/"${i}".unique.vcf > unique_vcfs/"${i}"_unique.genotypes.txt  
done < ${SAMPLES}

#unexposed lines
OUT=vcf/unexposed
mkdir -p $OUT
gatk SelectVariants -R ns.fasta -V vcf/genotyped.vcf -O $OUT/unexposed.vcf \
  -sn A10237_1 \
  -sn A10238_1 \
  -sn A10239_1 \
  -sn A10240_1 \
  -sn A10241_1 \
  -sn A10242_1 \
  -sn A10243_1 \
  -sn A10244_1 \
  -sn A10254_1 \
  --selectExpressions 'vc.getGenotype("A10254_1").isHomRef()'

gatk SelectVariants -R ns.fasta -V $OUT/unexposed.vcf -O $OUT/snps.vcf \
  --select-type-to-include SNP

#filter at INFO level and further filter on FORMAT level for genotype quality and read depth
gatk VariantFiltration -R ns.fasta -V $OUT/snps.vcf -O $OUT/filtered.snps.vcf \
	  --filter-name "QD5" \
	  --filter-expression "QD < 5.0" \
	  --filter-name "FS60" \
	  --filter-expression "FS > 60.0" \
    --filter-name "SOR3" \
	  --filter-expression "SOR > 3.0" \
	  --filter-name "MQ50" \
	  --filter-expression "MQ < 55.0" \
	  --filter-name "MQRS-12.5" \
	  --filter-expression "MQRankSum < -12.5" \
  	--filter-name "RPRS-8" \
	  --filter-expression "ReadPosRankSum < -8.0" 
     
#filter at format level
gatk VariantFiltration -R ns.fasta -V $OUT/filtered.snps.vcf -O $OUT/gt.filtered.snps.vcf \
	  --genotype-filter-name "lowGQ" \
    --genotype-filter-expression "GQ < 20" \
    --genotype-filter-name "lowDP" \
    --genotype-filter-expression "DP < 10" 
    
#convert FORMAT filters to no call because downstream tools don't parse the FORMAT level filter
gatk SelectVariants -R ns.fasta -V $OUT/gt.filtered.snps.vcf --set-filtered-gt-to-nocall -O $OUT/nocall.gt.filtered.snps.vcf

#remove sites that have a no call in one sample or more
gatk SelectVariants -R ns.fasta -V $OUT/nocall.gt.filtered.snps.vcf --max-nocall-number 0 --exclude-filtered TRUE -O $OUT/snps.final.vcf
bcftools stats $OUT/snps.final.vcf > $OUT/snps.final.vcf.stats.txt 
#do the same for all variants
gatk VariantFiltration -R ns.fasta -V $OUT/unexposed.vcf -O $OUT/unexposed.gt.filtered.vcf \
    --genotype-filter-name "lowDP" \
    --genotype-filter-expression "DP < 10" 
    
gatk SelectVariants -R ns.fasta -V $OUT/unexposed.gt.filtered.vcf --set-filtered-gt-to-nocall -O $OUT/unexposed.nocall.gt.filtered.vcf

gatk SelectVariants -R ns.fasta -V $OUT/unexposed.nocall.gt.filtered.vcf --max-nocall-number 0 --exclude-filtered TRUE -O $OUT/unexposed.final.vcf

bcftools stats $OUT/unexposed.final.vcf > $OUT/unexposed.final.vcf.stats.txt

cd $OUT
mkdir unique_vcfs
while read -r i; do
  bcftools view --private --samples $i snps.final.vcf > unique_vcfs/"${i}".unique.vcf
  bcftools query -f '%CHROM\t%POS[\t%GT\t]\n' unique_vcfs/"${i}".unique.vcf > unique_vcfs/"${i}"_unique.genotypes.txt  
done < ${SAMPLES}
