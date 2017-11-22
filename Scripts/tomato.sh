#!/bin/bash

PRJDIR=/media/jwu/data2/tomato
GE_REF=/media/jwu/data2/RefGene/tomato/tomato.SL2.50.fa
PICARD=/home/jwu/bin/picard.jar

# bwa index ${GE_REF}

cd ${PRJDIR}
# wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR527/SRR5274882/SRR5274882.sra
# wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR527/SRR5274881/SRR5274881.sra
# wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR527/SRR5274880/SRR5274880.sra

# fastq-dump -F --skip-technical --split-files  --defline-qual '+' --defline-seq '@$ac-$si/$ri' --split-3 SRR5274880.sra
# fastq-dump -F --skip-technical --split-files  --defline-qual '+' --defline-seq '@$ac-$si/$ri' --split-3 SRR5274881.sra
# fastq-dump -F --skip-technical --split-files  --defline-qual '+' --defline-seq '@$ac-$si/$ri' --split-3 SRR5274882.sra

# ### mapping    
# bwa mem -M -t 20 -R "@RG\tID:ejmt\tSM:ejmtpool\tPL:Illumina" ${GE_REF} SRR5274882_1.fastq SRR5274882_2.fastq > ejmt.sam
# ### convert sort index remove_duplicate index
# samtools view -bS -o ejmt.bam ejmt.sam
samtools sort -m 2G -@ 20 ejmt.bam -o ejmt.sorted.bam
samtools index ejmt.sorted.bam
### MarkDuplicates
java -Xmx20G -jar ${PICARD} MarkDuplicates I=ejmt.sorted.bam O=ejmt.sorted.redup.bam \
METRICS_FILE=ejmt.sort.metrics REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

samtools index ejmt.sorted.redup.bam
## vaiant calling
samtools  mpileup -DSug -C 50 -Q 20 -q 40 -f ${GE_REF} ejmt.sorted.redup.bam | bcftools view -cvg - > ejmt.vcf



# ### mapping    
# bwa mem -M -t 20 -R "@RG\tID:jmt\tSM:jmtpool\tPL:Illumina" ${GE_REF} SRR5274881_1.fastq SRR5274881_2.fastq > jmt.sam
# ### convert sort index remove_duplicate index
# samtools view -bS -o jmt.bam jmt.sam
samtools sort -m 2G -@ 20 jmt.bam -o jmt.sorted.bam
samtools index jmt.sorted.bam
### MarkDuplicates
java -Xmx20G -jar ${PICARD} MarkDuplicates I=jmt.sorted.bam O=jmt.sorted.redup.bam \
METRICS_FILE=jmt.sort.metrics REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

samtools index jmt.sorted.redup.bam
## vaiant calling
samtools  mpileup -DSug -C 50 -Q 20 -q 40 -f ${GE_REF} jmt.sorted.redup.bam | bcftools view -cvg - > jmt.vcf



# ### mapping    
# bwa mem -M -t 20 -R "@RG\tID:wt\tSM:wtpool\tPL:Illumina" ${GE_REF} SRR5274880_1.fastq SRR5274880_2.fastq > wt.sam
# ### convert sort index remove_duplicate index
# samtools view -bS -o wt.bam wt.sam
samtools sort -m 2G -@ 20 wt.bam -o wt.sorted.bam
samtools index wt.sorted.bam
### MarkDuplicates
java -Xmx20G -jar ${PICARD} MarkDuplicates I=wt.sorted.bam O=wt.sorted.redup.bam \
METRICS_FILE=wt.sort.metrics REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

samtools index wt.sorted.redup.bam
## vaiant calling
samtools  mpileup -DSug -C 50 -Q 20 -q 40 -f ${GE_REF} wt.sorted.redup.bam | bcftools view -cvg - > wt.vcf
