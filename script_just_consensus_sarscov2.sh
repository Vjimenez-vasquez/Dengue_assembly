#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#2# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz

#5# obtener el genoma consenso con las siguientes caracteristicas#
#score: 25, frecuencia de nucleotido predominante: 60%#
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 10 ;
done ; 

#6# juntar los genomas generados en un solo archivo y visualizarlos en aliview#
rm *qual.txt *fastq.gz.fa ;
mkdir genomes_mapping_m5 ;
mv *.fa genomes_mapping_m5 ;
cp reference.fasta genomes_mapping_m5 ;
cd genomes_mapping_m5 ;
cat *.fa > genomes.fasta ;
mafft --thread 15 --addfragments genomes.fasta --adjustdirection --auto --inputorder reference.fasta > genome_alin.fasta ;
aliview genome_alin.fasta ;
exit
