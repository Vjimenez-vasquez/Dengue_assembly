#!/bin/bash
# -*- ENCODING: UTF-8 -*-
bwa index GCF_000871845.1_2.fasta ;
for r1 in *fastq 
do
prefix=$(basename $r1 _1.fastq)
r2=${prefix}_2.fastq
bwa mem -t 10 GCF_000871845.1_2.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -bS -T GCF_000871845.1_2.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index ${prefix}.bam ;
rm ${prefix}_uno.bam ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix}.fasta -q 25 -t 0.55 ;
samtools fastq -F 4 ${prefix}.bam -1 ${prefix}_f_aligned.fastq -2 ${prefix}_r_aligned.fastq ;
done ;
rm *fastq.fa *qual.txt ;
cat *.fa > genomes_obtained_by_mapping.fasta ;
aliview genomes_obtained_by_mapping.fasta ;
mkdir aligned_for_denovo ;
mv *aligned.fastq aligned_for_denovo ;
cd aligned_for_denovo ;
for r1 in *fastq
do
prefix=$(basename $r1 _f_aligned.fastq)
r2=${prefix}_r_aligned.fastq
fastq_pair $r1 $r2 ; 
rm *.single.fq ;
done ;
for r1 in *fq
do
prefix=$(basename $r1 _f_aligned.fastq.paired.fq)
r2=${prefix}_r_aligned.fastq.paired.fq
spades --pe1-1 $r1 --pe1-2 $r2 -o ${prefix} ; 
done 
exit
