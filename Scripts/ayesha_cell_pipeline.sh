#!/bin/bash
# ref: http://www.cell.com/cell/fulltext/S0092-8674(17)30486-5

# no space in the directory path
PRJDIR=/media/jwu/data2/Ayesha_Project
GE_REF=/media/jwu/data2/RefGene/hg38ref/hg38.fa
RESULT=pipe2
PICARD=/home/jwu/bin/picard.jar

find ${PRJDIR} -name "*.fastq.gz" -printf "%h\n" | uniq | parallel mkdir -p {}/${RESULT}

echo "step 1.  time: $(date)" > ${PRJDIR}/${RESULT}.log
# /media/jwu/data2/Ayesha_Project/family2/101401_DVDS_v1.4.7.4_ACE3_mother/FASTQ/HLWTJADXX_101401_TAGGATGA_L002_R1_001.fastq.gz
# @RG\tID:flowcell.sample.lane\tLB:flowcell.sample\tPL:ILLUMINA\tPM:HISEQ\tSM:sample\tPU:flowcell.lane.sample
find ${PRJDIR} -name "*.fastq.gz" -printf "%h/_%f\n" | grep -v _R2_ | \
parallel --colsep '_' "bwa mem -M -R '@RG\tID:{7}.{8}.{10}\tLB:{7}.{8}\tPL:ILLUMINA\tPM:HISEQ\tSM:{8}\tPU:{7}.{10}.{8}' ${GE_REF} {1}_{2}_{3}_{4}_{5}_{6}{7}_{8}_{9}_{10}_R1_{12} {1}_{2}_{3}_{4}_{5}_{6}{7}_{8}_{9}_{10}_R2_{12} > {1}_{2}_{3}_{4}_{5}_{6}${RESULT}/{7}_{8}_{9}_{10}_{12}.script.sam"


echo "step 2.  time: $(date)" >> ${PRJDIR}/${RESULT}.log
find ${PRJDIR} -name "*.script.sam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} SortSam INPUT={} OUTPUT={.}.bam SORT_ORDER=coordinate"
# find ${PRJDIR} -type f -name "*.script.sam" -exec rm -f {} \;

find ${PRJDIR} -name "*.script.bam" -printf "%h\n" | uniq | grep /${RESULT}/ | parallel "samtools merge {}/sorted.bam {}/*.script.bam"
# find ${PRJDIR} -type f -name "*.script.bam" -exec rm -f {} \;

find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "samtools index {}"

### MarkDuplicates
find ${PRJDIR} -name "sorted.bam" | grep /${RESULT}/ | parallel "java -jar ${PICARD} MarkDuplicates I={} O={//}/redup.bam METRICS_FILE={//}/metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
find ${PRJDIR} -name "redup.bam" | grep /${RESULT}/ | parallel "samtools index {}"

### variant calling
find ${PRJDIR} -name "redup.bam" | grep /${RESULT}/ | parallel "samtools  mpileup -DSug -C 50 -Q 20 -q 40 -f ${GE_REF} {} | bcftools view -cvg - > {//}/variant.vcf"

### filte indel 
find ${PRJDIR} -name "variant.vcf" | grep /${RESULT}/ | parallel "grep -v '##'' {} | grep -v 'INDEL' > {//}/snp.vcf"