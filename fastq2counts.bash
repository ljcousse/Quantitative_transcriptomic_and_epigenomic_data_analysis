# !/bin/sh 

#########################
##                     ##
##  fastq2counts.bash  ##
##                     ##
#########################

# This bash script performs preprocessing, trimming, alignment and summarization into count files. 
# Most of it is written as a parallelized to decrease run time. Note that you require the following software installed and added to PATH:
# samtools
# dos2unix
# STAR (here installed at /data/louisc/STAR-2.5.2b/source/STAR)
# wget
# Trimmomatic (here installed at /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/)
# htseq-count

# To run this code we used a 64 core server with 5GB RAM/core (= 320GB RAM). When running this code as such you need at least 12 cores and 30-40GB RAM. 
# Tune down parallelization if you do not have these resources. Indexing the reference genome is the most memory intensive step, this can be circumvented 
# by downloading an readily available pre-made index.

# Note: don't forget to download the Acc_List.txt file from the github folder and put it in the directory you create for this analysis.
# Optional code, e.g. to obtain untrimmed results will be included, however will remain "commented out" since this code is not required for the final analysis.



#############
# Resources #
#############

# Create output directory and change directory 
####################
mkdir /data/louisc/Quantitative_transcriptomic_and_epigenomic_data_analysis
cd /data/louisc/Quantitative_transcriptomic_and_epigenomic_data_analysis/

# Download the reference genome
####################
wget https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz

# Download genome annotation
####################
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

# Unzip resources
####################
pigz -d *.gz

# Create indexed genome for alignment
####################
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /home/louisc/Quantitative_transcriptomic_and_epigenomic_data_analysis/99bpOverhang --genomeFastaFiles /home/louisc/Quantitative_transcriptomic_and_epigenomic_data_analysis/Homo_sapiens.GRCh37.75.dna.toplevel.fa --sjdbGTFfile /home/louisc/Quantitative_transcriptomic_and_epigenomic_data_analysis/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99



############
# Analysis #
############

# Convert accession list to a unix-readable file
####################
dos2unix Acc_List.txt
# the Acc_list.txt file contains all run ids present in the data

# Define parallelization
####################
amount=6
# amount of parallel cores

# amount of iterations is equal to the amount of samples divided by the amount of parallel cores rounded up to the next integer
for iteration in 1 #2 3 4 5 6 7 8
do 
	echo "$iteration"
	if [ "$iteration" -eq 8 ]; then
		start=$(($amount*($iteration-1)+1))
		end=44 # adapt to the total amount of samples
	else
		start=$(($amount*($iteration-1)+1))
		end=$(($amount*($iteration)))
		
	fi
	
	echo $start
	echo $end
	
	# Download fastq files
	####################
	echo $(head -$(($amount*$iteration)) Acc_List.txt | tail -$amount) > temp
	for accession in $(cat temp)
	do
		fastq-dump -I --split-files $accession &
	done
	wait
	
	# Fastqc for raw data
	####################
	echo $(ls | grep .fastq$)
	for fastq in $(ls | grep .fastq$)
	do 
		fastqc $fastq 
	done
	
	# Alignment with STAR (untrimmed)
	####################
	echo $(cat temp)
	for sample in $(cat temp)
	do
		/data/louisc/STAR-2.5.2b/source/STAR --runThreadN 8 --genomeDir "99bpOverhang/" --readFilesIn $sample\_1.fastq $sample\_2.fastq --outFileNamePrefix $sample --sjdbGTFfile Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	done
	
	# Trimming with Trimmomatic
	####################
	# Performs trimming for quality scores (sliding window) and adapter sequences based on Adapters.fa file simultaneously
	echo $(cat temp)
	for sample in $(cat temp)
	do
		java -jar /data/louisc/Hoofdstuk_Wim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 $sample\_1.fastq $sample\_2.fastq $sample\_1_trim_paired.fastq $sample\_1_trim_unpaired.fastq $sample\_2_trim_paired.fastq $sample\_2_trim_unpaired.fastq SLIDINGWINDOW:4:15 ILLUMINACLIP:Adapters.fa:2:30:10
	done
	
	# Alignment with STAR (trimmed)
	####################
	echo $(cat temp)
	for sample in $(cat temp)
	do
		/data/louisc/STAR-2.5.2b/source/STAR --runThreadN 8 --genomeDir /data/igenomes/99bpOverhang --readFilesIn /data/louisc/Hoofdstuk_Wim/$sample\_1_trim_paired.fastq /data/louisc/Hoofdstuk_Wim/$sample\_2_trim_paired.fastq --outFileNamePrefix $sample\_trim_bis --sjdbGTFfile /data/louisc/Hoofdstuk_Wim/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 99
	done
	
	# Data Summary: htseq-count
	####################
	# Gene level (untrimmed)
	# echo $(ls | grep Aligned.out.sam$ | grep -v trim)
	# for file in $(ls | grep Aligned.out.sam$ | grep -v trim) 
	# do
		# htseq-count -s reverse -f sam $file Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj.txt &
	# done
	# wait
	
	# Exon level (untrimmed)
	# echo $(ls | grep Aligned.out.sam$ | grep -v trim)
	# for file in $(ls | grep Aligned.out.sam$ | grep -v trim) 
	# do
		# htseq-count -s reverse -i exon_id -f sam $file Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj_exon.txt &
	# done
	# wait
	
	# Gene level (trimmed)
	echo $(ls | grep Aligned.out.sam$ | grep trim)
	for file in $(ls | grep Aligned.out.sam$ | grep trim) 
	do
		htseq-count -s reverse -f sam $file Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj.txt &
	done
	wait
	
	# Exon level (trimmed)
	# echo $(ls | grep Aligned.out.sam$ | grep trim)
	# for file in $(ls | grep Aligned.out.sam$ | grep trim) 
	# do
		# htseq-count -s reverse -i exon_id -f sam $file Homo_sapiens.GRCh37.75.gtf > $file\_htseqobj_exon.txt &
	# done
	# wait
	
	# Removal of temporary files (to limit disk space usage)
	####################
	rm *.fastq
	rm *Aligned.out.sam
	rm -r *fastqc
	rm -r *STARgenome
	rm *Log.out
	rm *Log.progress.out
	rm *SJ.out.tab
	rm temp
done

# Concatenate sepearate htseq-count output files
####################

# Gene level (untrimmed)
# cut -d$'\t' -f1 $(ls | grep htseqobj.txt | grep -v trim | head -1) > counts.txt
# echo $(ls | grep htseqobj.txt | grep -v trim)
# for file in $(ls | grep htseqobj.txt | grep -v trim)
# do
	# pr -mTs $'\t' <(cat counts.txt) <(cut -d$'\t' -f2 $file) > counts_temp.txt
	# cat counts_temp.txt > counts.txt
	# rm counts_temp.txt
# done

# Exon level (untrimmed)
# cut -d$'\t' -f1 $(ls | grep sam_htseqobj_exon.txt | grep -v trim | head -1) > counts_exon.txt
# echo $(ls | grep htseqobj_exon.txt | grep -v trim)
# for file in $(ls | grep htseqobj_exon.txt | grep -v trim)
# do
	# pr -mTs $'\t' <(cat counts_exon.txt) <(cut -d$'\t' -f2 $file) > counts_exon_temp.txt
	# cat counts_exon_temp.txt > counts_exon.txt
	# rm counts_exon_temp.txt
# done

# Gene level (trimmed)
cut -d$'\t' -f1 $(ls | grep htseqobj.txt | grep trim | head -1) > counts_trim.txt
echo $(ls | grep htseqobj.txt | grep trim)
for file in $(ls | grep htseqobj.txt | grep trim)
do
	pr -mTs $'\t' <(cat counts_trim.txt) <(cut -d$'\t' -f2 $file) > counts_trim_temp.txt
	cat counts_trim_temp.txt > counts_trim.txt
	rm counts_trim_temp.txt
done

# Exon level (trimmed)
# cut -d$'\t' -f1 $(ls | grep htseqobj_exon.txt | grep trim | head -1) > counts_exon_trim.txt
# echo $(ls | grep htseqobj_exon.txt | grep trim)
# for file in $(ls | grep htseqobj_exon.txt | grep trim)
# do
	# pr -mTs $'\t' <(cat counts_exon_trim.txt) <(cut -d$'\t' -f2 $file) > counts_exon_trim_temp.txt
	# cat counts_exon_trim_temp.txt > counts_exon_trim.txt
	# rm counts_exon_trim_temp.txt
# done
