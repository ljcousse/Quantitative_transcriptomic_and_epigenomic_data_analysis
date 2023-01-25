# !/bin/sh 
# Run on midas
# wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
# STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /home/louisc/Hoofdstuk_Wim/99bpOverhang --genomeFastaFiles /data/igenomes/Homo_Sapiens_GRCh37/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /home/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
# Get indexed genome from midas
# scp -rp -P 2222 louisc@midas.ugent.be:/home/louisc/Hoofdstuk_Wim/99bpOverhang /data/igenomes 


# # fastq to sam
cd /data/louisc/Hoofdstuk_Wim/
# pigz -d *

dos2unix Acc_List.txt

amount=6

for iteration in 1 #2 3 4 5 6 7 8
do 
	echo "$iteration"
	if [ "$iteration" -eq 8 ]; then
		start=$(($amount*($iteration-1)+1))
		end=44
	else
		start=$(($amount*($iteration-1)+1))
		end=$(($amount*($iteration)))
		
	fi
	
	echo $start
	echo $end
	
	# Download fastq files
	echo $(head -$(($amount*$iteration)) Acc_List.txt | tail -$amount) > temp
	for accession in $(cat temp)
	do
		fastq-dump -I --split-files $accession &
	done
	wait
	
	# Fastqc
	echo $(ls | grep .fastq$)
	for fastq in $(ls | grep .fastq$)
	do 
		fastqc $fastq 
	done
	
	# Map with STAR (untrimmed)
	echo $(cat temp)
	for sample in $(cat temp)
	do
		/data/louisc/STAR-2.5.2b/source/STAR --runThreadN 8 --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/$sample\_1.fastq /data/louisc/Hoofdstuk_Wim/$sample\_2.fastq --outFileNamePrefix $sample --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	done
	
	# Trim with Trimmomatic
	echo $(cat temp)
	for sample in $(cat temp)
	do
		java -jar /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 $sample\_1.fastq $sample\_2.fastq $sample\_1_trim_paired.fastq $sample\_1_trim_unpaired.fastq $sample\_2_trim_paired.fastq $sample\_2_trim_unpaired.fastq SLIDINGWINDOW:4:15 ILLUMINACLIP:Adapters.fa:2:30:10
	done
	
	# Map with STAR (trimmed)
	echo $(cat temp)
	for sample in $(cat temp)
	do
		/data/louisc/STAR-2.5.2b/source/STAR --runThreadN 8 --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/$sample\_1_trim_paired.fastq /data/louisc/Hoofdstuk_Wim/$sample\_2_trim_paired.fastq --outFileNamePrefix $sample\_trim_bis --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	done
	
	# java -jar /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads $amount -phred33 SRR5962198_1.fastq SRR5962198_2.fastq SRR5962198_1_trim_paired.fastq SRR5962198_1_trim_unpaired.fastq SRR5962198_2_trim_paired.fastq SRR5962198_2_trim_unpaired.fastq SLIDINGWINDOW:4:15
	# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN $amount --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/SRR5962198_1.fastq /data/louisc/Hoofdstuk_Wim/SRR5962198_2.fastq --outFileNamePrefix SRR5962198 --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 50
	# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN $amount --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/SRR5962198_1_trim_paired.fastq /data/louisc/Hoofdstuk_Wim/SRR5962198_2_trim_paired.fastq --outFileNamePrefix SRR5962198bis --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 50
	
	# java -jar /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads $amount -phred33 SRR5962198_1.fastq SRR5962198_2.fastq SRR5962198_1_trim_paired.fastq SRR5962198_1_trim_unpaired.fastq SRR5962198_2_trim_paired.fastq SRR5962198_2_trim_unpaired.fastq SLIDINGWINDOW:4:20
	# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN $amount --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/SRR5962198_1_trim_paired.fastq /data/louisc/Hoofdstuk_Wim/SRR5962198_2_trim_paired.fastq --outFileNamePrefix SRR5962198bisbis --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	
	# htseq-count -s reverse -f sam SRR5962198Aligned.out.sam /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > test\_htseqobj.txt
	
	# Data Summary: htseq-count
	# Gene level
	echo $(ls | grep Aligned.out.sam$ | grep -v trim)
	for file in $(ls | grep Aligned.out.sam$ | grep -v trim) 
	do
		htseq-count -s reverse -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj.txt &
	done
	wait
	# Exon level
	echo $(ls | grep Aligned.out.sam$ | grep -v trim)
	for file in $(ls | grep Aligned.out.sam$ | grep -v trim) 
	do
		htseq-count -s reverse -i exon_id -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj_exon.txt &
	done
	wait
	
	# Gene level (trimmed)
	echo $(ls | grep Aligned.out.sam$ | grep trim)
	for file in $(ls | grep Aligned.out.sam$ | grep trim) 
	do
		htseq-count -s reverse -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj.txt &
	done
	wait
	# Exon level (trimmed)
	echo $(ls | grep Aligned.out.sam$ | grep trim)
	for file in $(ls | grep Aligned.out.sam$ | grep trim) 
	do
		htseq-count -s reverse -i exon_id -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj_exon.txt &
	done
	wait
	
	rm *.fastq
	rm *Aligned.out.sam
	rm -r *fastqc
	rm -r *STARgenome
	rm *Log.out
	rm *Log.progress.out
	rm *SJ.out.tab
	rm temp
done


cut -d$'\t' -f1 $(ls | grep htseqobj.txt | grep -v trim | head -1) > counts.txt
echo $(ls | grep htseqobj.txt | grep -v trim)
for file in $(ls | grep htseqobj.txt | grep -v trim)
do
	pr -mTs $'\t' <(cat counts.txt) <(cut -d$'\t' -f2 $file) > counts_temp.txt
	cat counts_temp.txt > counts.txt
	rm counts_temp.txt
done

cut -d$'\t' -f1 $(ls | grep sam_htseqobj_exon.txt | grep -v trim | head -1) > counts_exon.txt
echo $(ls | grep htseqobj_exon.txt | grep -v trim)
for file in $(ls | grep htseqobj_exon.txt | grep -v trim)
do
	pr -mTs $'\t' <(cat counts_exon.txt) <(cut -d$'\t' -f2 $file) > counts_exon_temp.txt
	cat counts_exon_temp.txt > counts_exon.txt
	rm counts_exon_temp.txt
done

cut -d$'\t' -f1 $(ls | grep htseqobj.txt | grep trim | head -1) > counts_trim.txt
echo $(ls | grep htseqobj.txt | grep trim)
for file in $(ls | grep htseqobj.txt | grep trim)
do
	pr -mTs $'\t' <(cat counts_trim.txt) <(cut -d$'\t' -f2 $file) > counts_trim_temp.txt
	cat counts_trim_temp.txt > counts_trim.txt
	rm counts_trim_temp.txt
done

cut -d$'\t' -f1 $(ls | grep htseqobj_exon.txt | grep trim | head -1) > counts_exon_trim.txt
echo $(ls | grep htseqobj_exon.txt | grep trim)
for file in $(ls | grep htseqobj_exon.txt | grep trim)
do
	pr -mTs $'\t' <(cat counts_exon_trim.txt) <(cut -d$'\t' -f2 $file) > counts_exon_trim_temp.txt
	cat counts_exon_trim_temp.txt > counts_exon_trim.txt
	rm counts_exon_trim_temp.txt
done




# dos2unix Acc_List_sub.txt

# N=7
# (
# for accession in $(head - Acc_List_sub.txt)
# do
	# ((i=i%N)); ((i++==0)) && wait
	# fastq-dump -I --split-files $accession &
# done
# )


# echo $(ls | cut -d_ -f 1 | uniq)
# for sample in $(ls | grep "SRR" | cut -d_ -f 1 | uniq)
# do
	# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN 6 --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/$sample\_1.fastq /data/louisc/Hoofdstuk_Wim/$sample\_2.fastq --outFileNamePrefix $sample --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
# done


# N=7
# (
# for file in $(ls | grep Aligned.out.sam$) 
# do
	# ((i=i%N)); ((i++==0)) && wait
	# htseq-count -s no -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj.txt &
# done
# )

# sam to bam 
# N=7
# (
# for files in $(ls | grep Aligned.out.sam)
# do 
	# ((i=i%N)); ((i++==0)) && wait
	# samtools view -bS -o $files.bam $files &
# done
# )

# merge bams
# for prefix in $(ls | cut -d- -f 3 | cut -d_ -f 1 | sort -u)
# do 
	# samtools merge $prefix\_merged $(ls -rt -d -1 $PWD/{*,.*}/* | grep .bam | grep $prefix\_)
# done

# make counts
# for file in $(ls | grep bam) 
# do
	# samtools sort $file $file\_sorted
	
# done

# N=7
# (
# for file in $(ls | grep sorted.bam) 
# do
	# ((i=i%N)); ((i++==0)) && wait
	# htseq-count -s yes -f bam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh38.75.gtf > $file\_htseqobj.txt &
# done
# )




# N=7
# (
# for file in $(ls | grep sorted.bam) 
# do
	# ((i=i%N)); ((i++==0)) && wait
	# htseq-count -s no -f bam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh38.75.gtf > $file\_htseqobj.txt &
# done
# )


# pr -mTs $'\t' <(cut -d$'\t' -f2 001_htseqobj.txt) <(cut -d$'\t' -f2 002_htseqobj.txt) > htseqobj.txt
# cut -d$'\t' -f1 001_htseqobj.txt > Maastricht_annot.txt
