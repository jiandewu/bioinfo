#!/bin/bash

GATK=/home/jwu/Tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
PICARD=/home/jwu/bin/picard.jar
SNPEFF=/home/jwu/bin/snpEff/snpEff.jar

echo "step 2.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "*.script.sam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} SortSam INPUT={} OUTPUT={.}.bam SORT_ORDER=coordinate"
# find ${PRJDIR} -type f -name "*.script.sam" -exec rm -f {} \;

find ${PRJDIR} -name "*.script.bam" -printf "%h\n" | uniq | grep /${RESULT}/ | parallel "samtools merge {}/sorted.bam {}/*.script.bam"
# find ${PRJDIR} -type f -name "*.script.bam" -exec rm -f {} \;

echo "step 3.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} CollectAlignmentSummaryMetrics R=${GE_REF} I={} O={//}/alignment_metrics.txt"
find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} CollectInsertSizeMetrics INPUT={} OUTPUT={//}/insert_metrics.txt HISTOGRAM_FILE={//}/insert_size_histogram.pdf"
find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "samtools depth -a {} > {//}/depth_out.txt"
find ${PRJDIR} -name "depth_out.txt" | parallel "gzip {}"

echo "step 4.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} MarkDuplicates INPUT={} OUTPUT={//}/dedup.bam METRICS_FILE={//}/metrics.txt"
# find ${PRJDIR} -type f -name "sorted.bam" -exec rm -f {} \;

echo "step 5.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "dedup.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} BuildBamIndex INPUT={}"

echo "step 6.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "dedup.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T RealignerTargetCreator ${NEEDfixMisencodedQuals} -R ${GE_REF} -I {} -o {//}/realignment_targets.list"

echo "step 7.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "dedup.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T IndelRealigner ${NEEDfixMisencodedQuals} -R ${GE_REF} -I {} -targetIntervals {//}/realignment_targets.list -o {//}/realigned.bam"
# find ${PRJDIR} -type f -name "dedup.ba*" -exec rm -f {} \;

echo "step 8.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants.vcf"

echo "step 9.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_variants.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps.vcf"
find ${PRJDIR} -name "raw_variants.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels.vcf"

echo "step 10.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_snps.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o {//}/filtered_snps.vcf"

echo "step 11.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_indels.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels.vcf"

echo "step 12.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -o {//}/recal_data.table"

echo "step 13.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -BQSR {//}/recal_data.table -o {//}/post_recal_data.table"

echo "step 14.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "post_recal_data.table" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T AnalyzeCovariates -R ${GE_REF} -before {//}/recal_data.table -after {} -plots {//}/recalibration_plots.pdf"

echo "step 15.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T PrintReads -R ${GE_REF} -I {} -BQSR {//}/recal_data.table -o {//}/recal.bam"

echo "step 16.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "recal.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants_recal.vcf"

echo "step 17.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_variants_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps_recal.vcf"
find ${PRJDIR} -name "raw_variants_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels_recal.vcf"

echo "step 18.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_snps_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o {//}/filtered_snps_final.vcf"

echo "step 19.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "raw_indels_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels_final.vcf"

AWK_CMD='{if($0 ~ /\#/) print; else if($7 == "PASS") print}'
find ${PRJDIR} -name "filtered_*_final.vcf" | grep /${RESULT}/ | parallel "awk -F '\t' '${AWK_CMD}' {} > {.}_PASS.vcf"

echo "step 20.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "filtered_snps_final.vcf" | grep /${RESULT}/ | parallel "(cd {//} && java -jar ${SNPEFF} -v 'GRCh38.p7.RefSeq' {} > {//}/filtered_snps_final.ann.vcf)"

echo "step 21.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "recal.bam" | grep /${RESULT}/ | parallel "bedtools genomecov -bga -ibam {} > {//}/genomecov.bedgraph"