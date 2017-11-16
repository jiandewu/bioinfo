#!/bin/bash

PRODIR=/media/jwu/data2/Ayesha_Project/
GE_REF=/media/jwu/data2/Ayesha_Project/hg38ref/hg38.fa
GATK=/home/jwu/Tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
PICARD=/home/jwu/bin/picard.jar
SNPEFF=/home/jwu//bin/snpEff/snpEff.jar
RESULT=pipe

find ${PRODIR} -name "*.fastq.gz" -printf "%h\n" | uniq | parallel mkdir -p {}/${RESULT}

echo "step 1.  time: $(date)" > ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "*.fastq.gz" -printf "%h/_%f\n" | grep -v _R2_ | \
parallel --colsep '_' "bwa mem -M -R '@RG\tID:{7}.{8}.{10}\tLB:{7}.{8}\tPL:ILLUMINA\tPM:HISEQ\tSM:{8}\tPU:{7}.{10}.{8}' ${GE_REF} {1}_{2}_{3}_{4}_{5}_{6}{7}_{8}_{9}_{10}_R1_{12} {1}_{2}_{3}_{4}_{5}_{6}{7}_{8}_{9}_{10}_R2_{12} > {1}_{2}_{3}_{4}_{5}_{6}${RESULT}/{7}_{8}_{9}_{10}_{12}.sam"

echo "step 2.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "*.fastq.gz.sam" | parallel "java -jar ${PICARD} SortSam INPUT={} OUTPUT={.}.bam SORT_ORDER=coordinate"
find ${PRODIR} -type f -name "*.fastq.gz.sam" -exec rm -f {} \;

find ${PRODIR} -name "*.fastq.gz" -printf "%h\n" | uniq | parallel "samtools merge {}/sorted.bam {}/'*'.bam"
find ${PRODIR} -type f -name "*.fastq.gz.bam" -exec rm -f {} \;

echo "step 3.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "sorted.bam" | parallel "java -jar ${PICARD} CollectAlignmentSummaryMetrics R=${GE_REF} I={} O={//}/alignment_metrics.txt"
find ${PRODIR} -name "sorted.bam" | parallel "java -jar ${PICARD} CollectInsertSizeMetrics INPUT={} OUTPUT={//}/insert_metrics.txt HISTOGRAM_FILE={//}/insert_size_histogram.pdf"
find ${PRODIR} -name "sorted.bam" | parallel "samtools depth -a {} '>' {//}/depth_out.txt"

echo "step 4.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "sorted.bam" | parallel "java -jar ${PICARD} MarkDuplicates INPUT={} OUTPUT={//}/dedup.bam METRICS_FILE={//}/metrics.txt"
find ${PRODIR} -type f -name "sorted.bam" -exec rm -f {} \;

echo "step 5.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "dedup.bam" | parallel "java -jar ${PICARD} BuildBamIndex INPUT={}"

echo "step 6.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "dedup.bam" | parallel "java -jar ${GATK} -T RealignerTargetCreator -R ${GE_REF} -I {} -o {//}/realignment_targets.list"

echo "step 7.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "dedup.bam" | parallel "java -jar ${GATK} -T IndelRealigner -R ${GE_REF} -I {} -targetIntervals {//}/realignment_targets.list -o {//}/realigned.bam"
find ${PRODIR} -type f -name "dedup.bam" -exec rm -f {} \;

echo "step 8.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants.vcf"

echo "step 9.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "raw_variants.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps.vcf"
find ${PRODIR} -name "raw_variants.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels.vcf"

echo "step 10.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "raw_snps.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' {//}/filtered_snps.vcf"

echo "step 11.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "raw_indels.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels.vcf"

echo "step 12.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -o {//}/recal_data.table"

echo "step 13.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -BQSR {//}/recal_data.table -o {//}/post_recal_data.table"

# Step 14	Analyze Covariates
# gatk -T AnalyzeCovariates -R ../new/hg38.fa -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf

echo "step 15.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T PrintReads -R ${GE_REF} -I {} -BQSR {//}/recal_data.table -o {//}/recal.bam"

echo "step 16.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "recal.bam" | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants_recal.vcf"

echo "step 17.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "raw_variants_recal.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps_recal.vcf"
find ${PRODIR} -name "raw_variants_recal.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels_recal.vcf"

echo "step 18.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "raw_snps_recal.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o {//}/filtered_snps_final.vcf"

echo "step 19.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "raw_indels_recal.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels_final.vcf"

echo "step 20.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "filtered_snps_final.vcf" | parallel "java -jar ${SNPEFF} -v 'GRCh38.p7.RefSeq' {} > {//}/filtered_snps_final.ann.vcf"

echo "step 21.  time: $(date)" >> ${PRODIR}/${RESULT}.log
find ${PRODIR} -name "recal.bam" | parallel "bedtools genomecov -bga -ibam {} > {//}/genomecov.bedgraph"
