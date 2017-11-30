#!/bin/bash
# ref: https://gencore.bio.nyu.edu/variant-calling-pipeline/

# no space in the directory path
PRJDIR=/media/jwu/data3/Liu_project
GE_REF=/media/jwu/data2/RefGene/hg38ref/hg38.fa
RESULT=pipe
STARTSTEP=0

case $STARTSTEP in
0)
	find ${PRJDIR} -name "*.clean.fq.gz" -printf "%h\n" | uniq | parallel mkdir -p {}/${RESULT}
	;&
[0-1])
	echo "step 1.  time: $(date)" >${PRJDIR}/${RESULT}.log
	# /media/jwu/data3/Liu_project/F13TSFUSAT0104_HUMjxxX/result/1-N/clean_data/_131129_I649_FCC332NACXX_L1_RHUMjxxXAAAAAAA-16_1.clean.fq.gz
	# @RG\tID:flowcell.sample.lane\tLB:RHUMjxxXAAAAAAA-16\tPL:ILLUMINA\tPM:HISEQ\tSM:sample\tPU:flowcell.lane.sample
	find ${PRJDIR} -name "*.clean.fq.gz" -printf "%h/_%f\n" | grep -v _2.clean.fq.gz |
		parallel --colsep '_' "bwa mem -M -R '@RG\tID:{7}.{5}.{8}\tLB:{9}\tPL:ILLUMINA\tPM:HISEQ\tSM:{5}\tPU:{5}.{6}.{7}.{8}.{9}' ${GE_REF} {1}_{2}_{3}_{4}{5}_{6}_{7}_{8}_{9}_1.clean.fq.gz {1}_{2}_{3}_{4}{5}_{6}_{7}_{8}_{9}_2.clean.fq.gz > {1}_{2}_{3}_{4}${RESULT}/{5}_{6}_{7}_{8}_{9}.clean.fq.gz.script.sam"
	;;
esac

NEEDfixMisencodedQuals=-fixMisencodedQuals

source ${GITSRC}/bioinfo/Scripts/step2-21.sh
