#!/bin/bash

# please no space in the directory path
PRJDIR=/media/jwu/data3/Liu_project
GE_REF=/media/jwu/data2/Ayesha_Project/hg38ref/hg38.fa
GATK=/home/jwu/Tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
PICARD=/home/jwu/bin/picard.jar
SNPEFF=/home/jwu//bin/snpEff/snpEff.jar
RESULT=pipe

# find ${PRJDIR} -name "*.clean.fq.gz" -printf "%h\n" | uniq | parallel mkdir -p {}/${RESULT}

echo "step 1.  time: $(date)" > ${PRJDIR}/${RESULT}.log
# /media/jwu/data3/Liu_project/F13TSFUSAT0104_HUMjxxX/result/1-N/clean_data/131129_I649_FCC332NACXX_L1_RHUMjxxXAAAAAAA-16_1.clean.fq.gz
# @RG\tID:flowcell.sample.lane\tLB:flowcell.sample\tPL:ILLUMINA\tPM:HISEQ\tSM:sample\tPU:flowcell.lane.sample
find ${PRJDIR} -name "*.clean.fq.gz" | grep -v _2.clean.fq.gz | \
parallel --colsep '_' "bwa mem -M -R '@RG\tID:{6}.{8}.{7}\tLB:{6}.{8}\tPL:ILLUMINA\tPM:HISEQ\tSM:{8}\tPU:{6}.{7}.{8}' ${GE_REF} {1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}_1.clean.fq.gz {1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}_2.clean.fq.gz > {1}_{2}_{3}_{4}_{5}_{6}_{7}_{8}.clean.fq.gz.sam"

echo "step 2.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "*.clean.fq.gz.sam" | parallel "java -jar ${PICARD} SortSam INPUT={} OUTPUT={.}.bam SORT_ORDER=coordinate"
find ${PRJDIR} -type f -name "*.clean.fq.gz.sam" -exec rm -f {} \;

find ${PRJDIR} -name "*.clean.fq.gz.bam" -printf "%h\n" | uniq | parallel "samtools merge {}/sorted.bam {}/*.bam"
find ${PRJDIR} -type f -name "*.clean.fq.gz.bam" -exec rm -f {} \;

echo "step 3.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "sorted.bam" | parallel "java -jar ${PICARD} CollectAlignmentSummaryMetrics R=${GE_REF} I={} O={//}/alignment_metrics.txt"
find ${PRJDIR} -name "sorted.bam" | parallel "java -jar ${PICARD} CollectInsertSizeMetrics INPUT={} OUTPUT={//}/insert_metrics.txt HISTOGRAM_FILE={//}/insert_size_histogram.pdf"
find ${PRJDIR} -name "sorted.bam" | parallel "samtools depth -a {} > {//}/depth_out.txt"

echo "step 4.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "sorted.bam" | parallel "java -jar ${PICARD} MarkDuplicates INPUT={} OUTPUT={//}/dedup.bam METRICS_FILE={//}/metrics.txt"
find ${PRJDIR} -type f -name "sorted.bam" -exec rm -f {} \;

echo "step 5.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "dedup.bam" | parallel "java -jar ${PICARD} BuildBamIndex INPUT={}"

echo "step 6.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "dedup.bam" | parallel "java -jar ${GATK} -T RealignerTargetCreator -R ${GE_REF} -I {} -o {//}/realignment_targets.list"

echo "step 7.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "dedup.bam" | parallel "java -jar ${GATK} -T IndelRealigner -R ${GE_REF} -I {} -targetIntervals {//}/realignment_targets.list -o {//}/realigned.bam"
find ${PRJDIR} -type f -name "dedup.bam" -exec rm -f {} \;

echo "step 8.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants.vcf"

echo "step 9.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_variants.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps.vcf"
find ${PRJDIR} -name "raw_variants.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels.vcf"

echo "step 10.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_snps.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' {//}/filtered_snps.vcf"

echo "step 11.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_indels.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels.vcf"

echo "step 12.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -o {//}/recal_data.table"

echo "step 13.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -BQSR {//}/recal_data.table -o {//}/post_recal_data.table"

# Step 14	Analyze Covariates
# gatk -T AnalyzeCovariates -R ../new/hg38.fa -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf

echo "step 15.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | parallel "java -jar ${GATK} -T PrintReads -R ${GE_REF} -I {} -BQSR {//}/recal_data.table -o {//}/recal.bam"

echo "step 16.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "recal.bam" | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants_recal.vcf"

echo "step 17.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_variants_recal.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps_recal.vcf"
find ${PRJDIR} -name "raw_variants_recal.vcf" | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels_recal.vcf"

echo "step 18.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_snps_recal.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o {//}/filtered_snps_final.vcf"

echo "step 19.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_indels_recal.vcf" | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels_final.vcf"

echo "step 20.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "filtered_snps_final.vcf" | parallel "java -jar ${SNPEFF} -v 'GRCh38.p7.RefSeq' {} > {//}/filtered_snps_final.ann.vcf"

echo "step 21.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "recal.bam" | parallel "bedtools genomecov -bga -ibam {} > {//}/genomecov.bedgraph"
