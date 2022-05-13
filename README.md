![Captura desde 2022-05-07 13-05-18](https://user-images.githubusercontent.com/89874227/167266465-43cf8374-3d37-40b7-b385-f7002f162144.png)
# Dengue_assembly
Collection of codes for Dengue virus genome assembly

## Input
```r
1. ubuntu 22.04
2. reference genome
3. fastq files (forward and reverse)
```

## Requiremnts
```r
sudo apt install sra-toolkit
sudo apt install spades
sudo apt install fastqc
sudo apt install ivar
sudo apt install samtools 
sudo apt install fastq-pair
sudo apt install iva
aliview
```
## Just in case generate a naked script.sh file
```r
cd $HOME && touch script.sh && chmod +x script.sh ;
cd $HOME && echo '#!/bin/bash' > script.sh && echo '# -*- ENCODING: UTF-8 -*-' >> script.sh ;
```

## Details
```r
#1# indexar el genoma de referencia#
bwa index GCF_000871845.1_2.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fastq 
do
prefix=$(basename $r1 _1.fastq)
r2=${prefix}_2.fastq

#3# instrucciones para generar el archivo .bam#
bwa mem -t 10 GCF_000871845.1_2.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -bS -T GCF_000871845.1_2.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index ${prefix}.bam ;

#4# remover los archivos intermediarios#
rm ${prefix}_uno.bam ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;

#5# obtener el genoma consenso con las siguientes caracteristicas#
#score: 25, frecuencia de nucleotido predominante: 55%#
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix}.fasta -q 25 -t 0.55 -m 20 ;

#6# juntar los genomas generados en un solo archivo y visualizarlos en aliview#
mkdir genomes_mapping ;
mv *.fa genomes_mapping ;
rm genomes_mapping/*fastq.fa ;

#7# extraer del archivo .bam, todos los reads alineados con la referencia y generar dos fastq (f y r)#
samtools fastq -F 4 ${prefix}.bam -1 ${prefix}_f_aligned.fastq -2 ${prefix}_r_aligned.fastq ;
done ;

#8# remover archivos intermedios#
rm *fastq.fa *qual.txt ;

#9# mover estos nuevos fastq a una carpeta nueva para prerar el ensamblaje de novo#
mkdir aligned_for_denovo ;
mv *aligned.fastq aligned_for_denovo ;
cd aligned_for_denovo ;

#10# generar archivos fastq solo con reads pareados#
for r1 in *fastq
do
prefix=$(basename $r1 _f_aligned.fastq)
r2=${prefix}_r_aligned.fastq
fastq_pair $r1 $r2 ; 
rm *.single.fq ;
done ;

#11# realizar en ensamblaje de novo con spades#
for r1 in *fq
do
prefix=$(basename $r1 _f_aligned.fastq.paired.fq)
r2=${prefix}_r_aligned.fastq.paired.fq
spades --pe1-1 $r1 --pe1-2 $r2 -o ${prefix}_spades ;
mv ${prefix}_spades/scaffolds.fasta ${prefix}_spades/${prefix}_spades_scaffolds.fasta ;
mv ${prefix}_spades/${prefix}_spades_scaffolds.fasta . ;
done ;
rmdir *fq_spades ;
mv *scaffolds.fasta .. ;

#12# realizar el ensamblaje de novo con iva#
for r1 in *fq
do
prefix=$(basename $r1 _f_aligned.fastq.paired.fq)
r2=${prefix}_r_aligned.fastq.paired.fq
iva -t 4 -f $r1 -r $r2 ${prefix}_iva ;
mv ${prefix}_iva/contigs.fasta ${prefix}_iva/${prefix}_iva_contigs.fasta ;
mv ${prefix}_iva/${prefix}_iva_contigs.fasta . ;
done ;
mv *iva_contigs.fasta .. ;

#13# mover los ensamblados a carpetas especificas#
cd .. ;
mkdir iva_genomes ;
mv *iva_contigs.fasta iva_genomes ;
mkdir spades_genomes ;
mv *scaffolds.fasta spades_genomes ;
rm aligned_for_denovo/*aligned.fastq ;
```r

##14# modificar los respectivos headers y alinear con la referencia# 
```r
#para los contigs de IVA# 
cd iva_genomes ;
R 
library(seqinr) ;
r <- dir() ;
head <- gsub("_L001.*","",r) ; 
a <- 0 ;
for (i in 1:length(head)){ ;
a <- read.fasta(r[i]) ;
names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
write.fasta(a,names(a), file.out=paste0(head[i],"_1_iva_contigs.fasta")) ;
} ;
q("no") ;

cat *_1_iva_contigs.fasta > allcontigsiva.fasta ;
aliview allcontigsiva.fasta ;
mafft --addfragments allcontigsiva.fasta --adjustdirection --auto --inputorder GCF_000862125.1_1.fasta > iva_alineado.fasta ;
aliview iva_alineado.fasta ;

# para los contigs de SPADES#
cd spades_genomes ;
R 
library(seqinr) ;
r <- dir() ;
head <- gsub("_L001.*","",r) ;
a <- 0 ;
for (i in 1:length(head)){ ;
a <- read.fasta(r[i]) ;
names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
write.fasta(a,names(a), file.out=paste0(head[i],"_1_spades_scaffolds.fasta")) ;
} ;
q("no") ;

cat *_1_spades_scaffolds.fasta > allscaffolds.fasta ;
aliview allscaffolds.fasta ;
mafft --addfragments allscaffolds.fasta --adjustdirection --auto --inputorder GCF_000862125.1_1.fasta > spades_alineado.fasta ;
aliview spades_alineado.fasta ;
```

```
## Output
```r
1. Genomes obtained by mapping by reference (IVAR software)
2. Aligned reads in fastq forward and reverse different files
3. Genome scaffolds obtained by de-novo assembly (SPADES and IVA)
```
