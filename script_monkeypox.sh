#1#alineado con el genoma de referencia#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz

#2# instrucciones para generar el archivo .bam#
bwa mem -t 15 NC_063383.1.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 15 -bS -T NC_063383.1.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 15 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 15 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 15 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 15 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 15 ${prefix}.bam ;

#3#remocion de intermedios#
rm *_uno.bam *_uno.sam *_unoa.bam *_dosa.bam *_tresa.bam *_cuatroa.bam ;

#4#obtencion del genoma consenso#
samtools mpileup -aa -A -d 0 -Q 0 ${prefix}.bam | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 20 ; 

#5#obtencion de archivos bam con solo mapeados#
samtools view -b -F 4 ${prefix}.bam > ${prefix}.mapped.bam ; 
samtools fastq -1 ${prefix}_f.fq -2 ${prefix}_r.fq -0 /dev/null -s /dev/null -n ${prefix}.mapped.bam ; 
rm *.fastq.gz.fa *.fastq.gz.mapped.bam *.fastq.gz.qual.txt *.fastq.gz_f.fq *.fastq.gz_r.fq ;
samtools index -@ 15 ${prefix}.mapped.bam ; 
done ;

#6#mover los archivos#
mkdir mapped ; 
mv *.mapped.bam *.fa *.fq *.mapped.bam.bai mapped/ ;
cat mapped/*.fa > mapped/genomes.fasta ; 
aliview mapped/genomes.fasta ;
mafft --thread 15 --addfragments mapped/genomes.fasta --adjustdirection --auto --inputorder NC_063383.1.fasta > mapped/genomes_alin.fasta ;
aliview mapped/genome_alin.fasta ;
