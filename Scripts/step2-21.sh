#!/bin/bash

GATK=/home/jwu/Tools/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
PICARD=/home/jwu/bin/picard.jar
SNPEFF=/home/jwu/bin/snpEff/snpEff.jar

case $STARTSTEP in
[0-2])
	echo "step 2.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "*.script.sam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} SortSam INPUT={} OUTPUT={.}.bam SORT_ORDER=coordinate"
	# find ${PRJDIR} -type f -name "*.script.sam" -exec rm -f {} \;

	find ${PRJDIR} -name "*.script.bam" -printf "%h\n" | uniq | grep /${RESULT} | parallel "samtools merge {}/sorted.bam {}/*.script.bam"
	# find ${PRJDIR} -type f -name "*.script.bam" -exec rm -f {} \;
	;&
[0-3])
	echo "step 3.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} CollectAlignmentSummaryMetrics R=${GE_REF} I={} O={//}/alignment_metrics.txt"
	find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} CollectInsertSizeMetrics INPUT={} OUTPUT={//}/insert_metrics.txt HISTOGRAM_FILE={//}/insert_size_histogram.pdf"
	find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "samtools depth -a {} > {//}/depth_out.txt"
	find ${PRJDIR} -name "depth_out.txt" | parallel "gzip {}"
	;&
[0-4])
	echo "step 4.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} MarkDuplicates INPUT={} OUTPUT={//}/dedup.bam METRICS_FILE={//}/metrics.txt"
	# find ${PRJDIR} -type f -name "sorted.bam" -exec rm -f {} \;
	;&
[0-5])
	echo "step 5.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "dedup.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} BuildBamIndex INPUT={}"
	;&
[0-6])
	echo "step 6.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "dedup.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T RealignerTargetCreator ${NEEDfixMisencodedQuals} -R ${GE_REF} -I {} -o {//}/realignment_targets.list"
	;&
[0-7])
	echo "step 7.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "dedup.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T IndelRealigner ${NEEDfixMisencodedQuals} -R ${GE_REF} -I {} -targetIntervals {//}/realignment_targets.list -o {//}/realigned.bam"
	# find ${PRJDIR} -type f -name "dedup.ba*" -exec rm -f {} \;
	;&
[0-8])
	echo "step 8.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants.vcf"
	;&
[0-9])
	echo "step 9.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "raw_variants.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps.vcf"
	find ${PRJDIR} -name "raw_variants.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels.vcf"
	;&
[0-10])
	echo "step 10.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "raw_snps.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o {//}/filtered_snps.vcf"
	;&
[0-11])
	echo "step 11.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "raw_indels.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels.vcf"
	;&
[0-12])
	echo "step 12.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -o {//}/recal_data.table"
	;&
[0-13])
	echo "step 13.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T BaseRecalibrator -R ${GE_REF} -I {} -knownSites {//}/filtered_snps.vcf -knownSites {//}/filtered_indels.vcf -BQSR {//}/recal_data.table -o {//}/post_recal_data.table"
	;&
[0-14])
	echo "step 14.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "post_recal_data.table" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T AnalyzeCovariates -R ${GE_REF} -before {//}/recal_data.table -after {} -plots {//}/recalibration_plots.pdf"
	;&
[0-15])
	echo "step 15.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "realigned.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T PrintReads -R ${GE_REF} -I {} -BQSR {//}/recal_data.table -o {//}/recal.bam"
	;&
[0-16])
	echo "step 16.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "recal.bam" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T HaplotypeCaller -R ${GE_REF} -I {} -o {//}/raw_variants_recal.vcf"
	;&
[0-17])
	echo "step 17.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "raw_variants_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType SNP -o {//}/raw_snps_recal.vcf"
	find ${PRJDIR} -name "raw_variants_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T SelectVariants -R ${GE_REF} -V {} -selectType INDEL -o {//}/raw_indels_recal.vcf"
	;&
[0-18])
	echo "step 18.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "raw_snps_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName 'basic_snp_filter' -o {//}/filtered_snps_final.vcf"
	;&
[0-19])
	echo "step 19.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "raw_indels_recal.vcf" | grep /${RESULT}/ | parallel "java -jar ${GATK} -T VariantFiltration -R ${GE_REF} -V {} --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName 'basic_indel_filter' -o {//}/filtered_indels_final.vcf"

	AWK_CMD='{if($0 ~ /\#/) print; else if($6 > 30 && $7 == "PASS") print}'
	find ${PRJDIR} -name "filtered_*_final.vcf" | grep /${RESULT}/ | parallel "awk -F '\t' '${AWK_CMD}' {} > {.}_PASS_Qgt30.vcf"
	;&
[0-20])
	echo "step 20.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "filtered_snps_final.vcf" | grep /${RESULT}/ | parallel "(cd {//} && java -jar ${SNPEFF} -v 'GRCh38.p7.RefSeq' {} > {//}/filtered_snps_final.ann.vcf)"
	;&
[0-21])
	echo "step 21.  time: $(date)" >>${LOGFILE}
	find ${PRJDIR} -name "recal.bam" | grep /${RESULT}/ | parallel "bedtools genomecov -bga -ibam {} > {//}/genomecov.bedgraph"
	;;
esac
