#!/bin/bash

## Directories + programs + files
ref='/mnt/ps/home/CORP/amanda.mitchell/genomes/human/hg19/hg19.fa'
genome='/mnt/ps/home/CORP/amanda.mitchell/genomes/human/hg19/'
fastqtovcf='/mnt/ps/home/CORP/amanda.mitchell/data/data/RXRX/scripts/fastqtovcf'
samtools='/mnt/ps/home/CORP/amanda.mitchell/tools/samtools-1.3.1/samtools'
PICARD='/mnt/ps/home/CORP/amanda.mitchell/miniconda3/pkgs/picard-2.18.23-0/bin/picard'
GATK='/mnt/ps/home/CORP/amanda.mitchell/miniconda3/pkgs/gatk-3.6-1/bin/gatk'
fastqc='/mnt/ps/home/CORP/amanda.mitchell/miniconda3/pkgs/fastqc-0.11.8-0/bin/fastqc'
logs='/mnt/bh1/data/bulkrnaseq/wgs/220607/logs'
fastq='/mnt/bh1/data/bulkrnaseq/wgs/220607/files'
results='/mnt/bh1/data/bulkrnaseq/wgs/220607/results'
vcf='/mnt/bh1/data/bulkrnaseq/wgs/220607/results/vcf'
plinkdir='/mnt/bh1/data/bulkrnaseq/wgs/220607/results/plink'
refp='/mnt/ps/home/CORP/amanda.mitchell/genomes/human/hg19/hg19.dict'
AMdir='/mnt/bh1/data/bulkrnaseq/wgs/220607'


for i in A1_CSFP220025719-1a_HY7H7DSX3 B1_CSFP220025720-1a_HY2NCDSX3 B1_CSFP220025720-1a_HY7LFDSX3 C1_CSFP220025721-1a_HY7H7DSX3;
do 
  bwa mem -t 32 -M ${ref} ${fastq}/${i}_L2_1.fq  ${fastq}/${i}_L2_2.fq > ${results}/${i}_paired.sam
done  


for i in C1_CSFP220025721-1a_HY7H7DSX3 B1_CSFP220025720-1a_HY7LFDSX3 B1_CSFP220025720-1a_HY2NCDSX3 ;
do 
  ${samtools} view -h ${results}/${i}_paired.sam -o ${results}/${i}_paired.bam #30 minutes, 95 files 
  ${samtools} sort ${results}/${i}_paired.bam -o ${results}/${i}_paireds.bam #10 minutes 95 files
  $PICARD MarkDuplicates I=${results}/${i}_paireds.bam O=${results}/${i}_dedup.bam M=${logs}/${i}_pmarked_dup_metrics.txt #34 minutes
  ${samtools} addreplacerg -r '@RG\tID:samplename\tSM:samplename' -o ${results}/${i}_sdh.bam ${results}/${i}_dedup.bam #10 minutes
  ${samtools} view -h ${results}/${i}_sdh.bam -o ${results}/${i}_final.bam #10 minutes
  $PICARD BuildBamIndex I=${results}/${i}_final.bam 
  $samtools flagstat ${results}/${i}_final.bam > ${logs}/${i}_bam_logstat.txt #1 minute
  mv ${results}/${i}.bam ${results}/${i}_final.bam
  mv ${results}/${i}.bai ${results}${i}_final.bai
  rm ${results}/${i}*paired.bam ${results}/${i}*paireds.bam ${results}/${i}_ded* ${results}/${i}_sdh*
done



for i in B1_CSFP220025720-1a_HY2NCDSX3;
do 
gatk HaplotypeCaller -R ${ref} -I ${results}/${i}_final.bam -O ${vcf}/${i}_raw_variants.vcf.gz -ERC GVCF #426 minutes
gunzip ${vcf}/${i}_raw_variants.vcf.gz
gatk SelectVariants -R ${ref} -V ${vcf}/${i}_raw_variants.vcf -select-type SNP -O ${vcf}/${i}_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${vcf}/${i}_raw_variants.vcf -select-type INDEL -O ${vcf}/${i}_raw_indels.vcf

# Filter SNPs
gatk VariantFiltration -R ${ref} -V ${vcf}/${i}_raw_snps.vcf -O ${vcf}/${i}_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

# Filter Indels
gatk VariantFiltration -R ${ref} -V ${vcf}/${i}_raw_indels.vcf -O ${vcf}/${i}_filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

# Exclude filtered variants
gatk SelectVariants --exclude-filtered -V ${vcf}/${i}_filtered_snps.vcf -O ${vcf}/${i}_bqsr_snps.vcf
gatk SelectVariants --exclude-filtered -V ${vcf}/${i}_filtered_indels.vcf -O ${vcf}/${i}_bqsr_indels.vcf
$PICARD AddOrReplaceReadGroups I=${results}/${i}_final3.bam O=${results}/${i}_final4.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
gatk BaseRecalibrator -I ${results}/${i}_final4.bam -R ${ref} --known-sites ${vcf}/${i}_bqsr_snps.vcf --known-sites ${vcf}/${i}_bqsr_indels.vcf -O ${vcf}/${i}_ecal_data.table
gatk ApplyBQSR -R ${ref} -I ${results}/${i}_final4.bam -bqsr ${vcf}/${i}_ecal_data.table -O ${vcf}/${i}_recal_reads.bam \
gatk BaseRecalibrator-R ${ref} -I ${vcf}/${i}_recal_reads.bam --known-sites ${vcf}/${i}_bqsr_snps.vcf --known-sites ${vcf}/${i}_bqsr_indels.vcf -O ${vcf}/${i}_post_recal_data.table
gatk AnalyzeCovariates -before ${vcf}/${i}_ecal_data.table -after ${vcf}/${i}_post_recal_data.table -plots ${vcf}/${i}_recalibration_plots.pdf
gatk HaplotypeCaller -R ${ref} -I ${vcf}/${i}_recal_reads.bam -O ${vcf}/${i}_raw_variants_recal.vcf #292 minutes

# Extract SNPs and Indels
gatk SelectVariants -R ${ref} -V ${vcf}/${i}_raw_variants_recal.vcf -select-type SNP -O ${vcf}/${i}_raw_snps_recal.vcf
gatk SelectVariants -R ${ref} -V ${vcf}/${i}_raw_variants.vcf -select-type INDEL -O ${vcf}/${i}_raw_indels_recal.vcf

# Filter SNPs
gatk VariantFiltration -R ${ref} -V ${vcf}/${i}_raw_snps_recal.vcf -O ${vcf}/${i}_filtered_snps_final.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        #-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        #-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

# Filter Indels
gatk VariantFiltration -R ${ref} -V ${vcf}/${i}_raw_indels_recal.vcf -O ${vcf}/${i}_filtered_indels_final.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

# Annotate SNPs and predict effects
snpEff -Xmx8g -v GRCh37.75 ${vcf}/${i}_filtered_snps_final.vcf > ${vcf}/${i}_filtered_snps_final.ann.vcf #hg19
snpEff -Xmx8g -v GRCh37.75 ${vcf}/${i}_filtered_indels_final.vcf > ${vcf}/${i}_filtered_indel_final.ann.vcf #hg19

done


echo "Finished."
