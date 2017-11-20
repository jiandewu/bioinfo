#!/bin/bash

# no space in the directory path
PRJDIR=/media/jwu/data2/Ayesha_Project
GE_REF=/media/jwu/data2/RefGene/hg38ref/hg38.fa
RESULT=pipe

find ${PRJDIR} -name "*.fastq.gz" -printf "%h\n" | uniq | parallel mkdir -p {}/${RESULT}

echo "step 1.  time: $(date)" > ${PRJDIR}/${RESULT}.log
# /media/jwu/data2/Ayesha_Project/family2/101401_DVDS_v1.4.7.4_ACE3_mother/FASTQ/HLWTJADXX_101401_TAGGATGA_L002_R1_001.fastq.gz
# @RG\tID:flowcell.sample.lane\tLB:flowcell.sample\tPL:ILLUMINA\tPM:HISEQ\tSM:sample\tPU:flowcell.lane.sample
find ${PRJDIR} -name "*.fastq.gz" -printf "%h/_%f\n" | grep -v _R2_ | \
parallel --colsep '_' "bwa mem -M -R '@RG\tID:{7}.{8}.{10}\tLB:{7}.{8}\tPL:ILLUMINA\tPM:HISEQ\tSM:{8}\tPU:{7}.{10}.{8}' ${GE_REF} {1}_{2}_{3}_{4}_{5}_{6}{7}_{8}_{9}_{10}_R1_{12} {1}_{2}_{3}_{4}_{5}_{6}{7}_{8}_{9}_{10}_R2_{12} > {1}_{2}_{3}_{4}_{5}_{6}${RESULT}/{7}_{8}_{9}_{10}_{12}.script.sam"

source ${GITSRC}/bioinfo/Scripts/step2-21.sh