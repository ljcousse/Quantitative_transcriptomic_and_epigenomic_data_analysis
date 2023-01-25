# !/bin/sh 

# Run this part on midas, rest is on porthos
# wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
# STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /home/louisc/Hoofdstuk_Wim/99bpOverhang --genomeFastaFiles /data/igenomes/Homo_Sapiens_GRCh37/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /home/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
# Get indexed genome from midas
# scp -rp -P 2222 louisc@midas.ugent.be:/home/louisc/Hoofdstuk_Wim/99bpOverhang /data/igenomes 


# # fastq to sam
cd /data/louisc/Hoofdstuk_Wim/
# pigz -d *

dos2unix Acc_List.txt

amount=6

for iteration in 1 2 3 4 5 6 7 8
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
		fastq-dump -I --split-files $accession 
	done
	
	# Fastqc
	# echo $(ls | grep .fastq$)
	# for fastq in $(ls | grep .fastq$)
	# do 
		# fastqc $fastq 
	# done
	
	# Map with STAR (untrimmed)
	# echo $(cat temp)
	# for sample in $(cat temp)
	# do
		# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN 8 --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/$sample\_1.fastq /data/louisc/Hoofdstuk_Wim/$sample\_2.fastq --outFileNamePrefix $sample --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	# done
	
	# Trim with Trimmomatic
	echo $(cat temp)
	for sample in $(cat temp)
	do
		java -jar /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 $sample\_1.fastq $sample\_2.fastq $sample\_1_trim_paired.fastq $sample\_1_trim_unpaired.fastq $sample\_2_trim_paired.fastq $sample\_2_trim_unpaired.fastq SLIDINGWINDOW:4:15 ILLUMINACLIP:/data/louisc/Hoofdstuk_Wim/Adapters.fa:2:30:10
	done
	
	# Map versus rRNA
	echo $(cat temp)
	for sample in $(cat temp)
	do
		/data/louisc/STAR-2.5.2b/source/STAR --runThreadN 8 --genomeDir /data/louisc/Hoofdstuk_Wim/rRNA/ --readFilesIn /data/louisc/Hoofdstuk_Wim/$sample\_1_trim_paired.fastq /data/louisc/Hoofdstuk_Wim/$sample\_2_trim_paired.fastq --outFileNamePrefix $sample\_rRNA --sjdbOverhang 99 --outReadsUnmapped Fastx
	done
	
	# Map unmapped
	echo $(cat temp)
	for sample in $(cat temp)
	do
		/data/louisc/STAR-2.5.2b/source/STAR --runThreadN 8 --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/$sample\_rRNAUnmapped.out.mate1 /data/louisc/Hoofdstuk_Wim/$sample\_rRNAUnmapped.out.mate2 --outFileNamePrefix $sample\_trim_bis_norRNA --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	done
	
	# java -jar /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads $amount -phred33 SRR5962198_1.fastq SRR5962198_2.fastq SRR5962198_1_trim_paired.fastq SRR5962198_1_trim_unpaired.fastq SRR5962198_2_trim_paired.fastq SRR5962198_2_trim_unpaired.fastq SLIDINGWINDOW:4:15
	# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN $amount --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/SRR5962198_1.fastq /data/louisc/Hoofdstuk_Wim/SRR5962198_2.fastq --outFileNamePrefix SRR5962198 --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 50
	# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN $amount --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/SRR5962198_1_trim_paired.fastq /data/louisc/Hoofdstuk_Wim/SRR5962198_2_trim_paired.fastq --outFileNamePrefix SRR5962198bis --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 50
	
	# java -jar /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads $amount -phred33 SRR5962198_1.fastq SRR5962198_2.fastq SRR5962198_1_trim_paired.fastq SRR5962198_1_trim_unpaired.fastq SRR5962198_2_trim_paired.fastq SRR5962198_2_trim_unpaired.fastq SLIDINGWINDOW:4:20
	# /data/louisc/STAR-2.5.2b/source/STAR --runThreadN $amount --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/SRR5962198_1_trim_paired.fastq /data/louisc/Hoofdstuk_Wim/SRR5962198_2_trim_paired.fastq --outFileNamePrefix SRR5962198bisbis --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	
	# htseq-count -s reverse -f sam SRR5962198Aligned.out.sam /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > test\_htseqobj.txt
	
	# Data Summary: htseq-count
	# Gene level
	# echo $(ls | grep Aligned.out.sam$ | grep -v trim)
	# for file in $(ls | grep Aligned.out.sam$ | grep -v trim) 
	# do
		# htseq-count -s reverse -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj.txt &
	# done
	# wait
	# Exon level
	# echo $(ls | grep Aligned.out.sam$ | grep -v trim)
	# for file in $(ls | grep Aligned.out.sam$ | grep -v trim) 
	# do
		# htseq-count -s reverse -i exon_id -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj_exon.txt &
	# done
	# wait
	
	# Gene level (trimmed)
	echo $(ls | grep Aligned.out.sam$ | grep trim_bis_norRNA)
	for file in $(ls | grep Aligned.out.sam$ | grep trim_bis_norRNA) 
	do
		htseq-count -s reverse -f sam $file /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf > $file\_norRNA_htseqobj.txt &
	done
	wait
	# Exon level (trimmed)
	
	
	rm *.fastq
	rm *Aligned.out.sam
	rm -r *fastqc
	rm -r *STARgenome
	rm *Log.out
	rm *Log.progress.out
	rm *SJ.out.tab
	rm temp
done
