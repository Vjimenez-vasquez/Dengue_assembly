# Dengue_assembly
Collection of codes for Dengue virus genome assembly

## Usage
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
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix}.fasta -q 25 -t 0.55 ;

#6# extraer del archivo .bam, todos los reads alineados con la referencia y generar dos fastq (f y r)#
samtools fastq -F 4 ${prefix}.bam -1 ${prefix}_f_aligned.fastq -2 ${prefix}_r_aligned.fastq ;
done ;

#7# remover archivos intermedios#
rm *fastq.fa *qual.txt ;

#8# juntar los genomas generados en un solo archivo y visualizarlos en aliview#
cat *.fa > genomes_obtained_by_mapping.fasta ;
aliview genomes_obtained_by_mapping.fasta ;

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

#11# realizar en ensamblade de novo con spades#
for r1 in *fq
do
prefix=$(basename $r1 _f_aligned.fastq.paired.fq)
r2=${prefix}_r_aligned.fastq.paired.fq
spades --pe1-1 $r1 --pe1-2 $r2 -o ${prefix} ; 
done 

#12# disfrutar#

```
## Output
```r
1. Genomes obtained by mapping by reference (IVAR software)
2. Aligned reads in fastq forward and reverse different files
3. Genome scaffolds obtained by de-novo assembly (SPADES)
```
