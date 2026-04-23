#!/bin/bash


#The NGS Pipeline - From raw data to alignement and variant calls


# --- Step 1: Download data ---
# Downloaing the raw fastq data, annotation files and reference files
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# Rename files due to "odd" extension (fastq files usually end in .gz so this will prevent any issues downstream as you use more tools)
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

# --- Step 2: Organise project and files ---
# Note: You should be in your home directory, use the pwd command to check where you are, and the ls command to ensure your files are also located in the home directory.
# Using -p creates the whole path at once
mkdir -p ~/ngs_assignment/dnaseq/data/untrimmed_fastq 
mkdir -p ~/ngs_assignment/dnaseq/data/trimmed_fastq 
mkdir -p ~/ngs_assignment/dnaseq/meta 
mkdir -p ~/ngs_assignment/dnaseq/results 
mkdir -p ~/ngs_assignment/dnaseq/logs 

# Define Directory Variables (Shortcuts) for the rest of the script
# This makes the script easier to read and maintain
RAW_DATA=~/ngs_assignment/dnaseq/data/untrimmed_fastq
TRIM_DATA=~/ngs_assignment/dnaseq/data/trimmed_fastq
RES_DIR=~/ngs_assignment/dnaseq/results

# Move files to their organized homes
mv *.gz $RAW_DATA/
mv annotation.bed ~/ngs_assignment/dnaseq/data/
mv hg19.fa.gz ~/ngs_assignment/dnaseq/data/


# --- Step 3: Install pipeline tools ---
# We assume Conda is already installed on the system
# We set up the channels first to ensure we find the bioinformatics tools
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge

# Install the tools
# -y is used to automatically answer 'yes' to prompts for non-interactive execution
#conda install -y samtools bwa freebayes picard trimmomatic fastqc
#sudo apt -y install libvcflib-tools
#sudo apt -y install bedtools


# --- Step 4: Pre-Alignment QC ---
# First we use fastqc to assess the raw data quality.
cd $RAW_DATA 
fastqc -t 4 *.fastq.gz

# Make a new directoy to store the results and move the results into it
mkdir $RES_DIR/fastqc_untrimmed_reads/
mv *fastqc* $RES_DIR/fastqc_untrimmed_reads/

# The fastqc results are currently stored in a zip(compressed) file, to access the summary of the results create a loop to unzip all of them at once.
cd $RES_DIR/fastqc_untrimmed_reads/
for zip in *.zip; do unzip $zip; done

# Concatenate the results and store them in the logs to keep a track that's easy to read later.
cat */summary.txt > ~/ngs_assignment/dnaseq/logs/fastqc_summaries.txt

# Next move onto trimming the data
# Note: ensure the adapter path is correct for the version you have. To check your adapter do: find ~/anaconda3 -name "NexteraPE-PE.fa"
trimmomatic PE \
        -threads 4 \
        -phred33 \
        $RAW_DATA/NGS0001.R1.fastq.gz $RAW_DATA/NGS0001.R2.fastq.gz \
        -baseout $TRIM_DATA/NGS0001_trimmed_R \
        ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.40-hdfd78af_0/share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa:2:30:10 \
        TRAILING:25 MINLEN:50

# Now use fastqc to assess the trimmed data quality.
cd $TRIM_DATA
fastqc -t 4 *_trimmed_R*P

# Make a new directoy to store the results and move the results into it
mkdir $RES_DIR/fastqc_trimmed_reads/
mv *fastqc* $RES_DIR/fastqc_trimmed_reads/

# Use the same loop to unzip all of them at once
cd $RES_DIR/fastqc_trimmed_reads/
for zip in *.zip; do unzip $zip; done

# Concatenate the results and store them in the logs to keep a track that's easy to read later.
cat */summary.txt > ~/ngs_assignment/dnaseq/logs/fastqc_trimmed_summaries.txt

# --- Step 5: Alignment ---
# Building the reference index:
# Create a folder for the reference and its index files and then run bwa to generate the index files
mkdir -p ~/ngs_assignment/dnaseq/data/reference
mv ~/ngs_assignment/dnaseq/data/hg19.fa.gz ~/ngs_assignment/dnaseq/data/reference
bwa index ~/ngs_assignment/dnaseq/data/reference/hg19.fa.gz

# Create a folder for your aligned data to be stored later
# Then run the BWA mem alignment - ensure the read group identifiers are correct
# Note: run zcat ~/ngs_assignment/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz | head -n 1 to check the fastq header to verify the read group info
mkdir -p ~/ngs_assignment/dnaseq/data/aligned_data
ALIGNED_DIR=~/ngs_assignment/dnaseq/data/aligned_data
bwa mem -t 4 -v 1 \
	-R '@RG\tID:11V6WR1.111.D1375ACXX.1\tSM:NGS0001\tPL:ILLUMINA\tLB:NGS0001-lib1\tDT:2017-02-23\tPU:D1375ACXX' \
	-I 250,50 \
	~/ngs_assignment/dnaseq/data/reference/hg19.fa.gz \
	$TRIM_DATA/NGS0001_trimmed_R_1P \
	$TRIM_DATA/NGS0001_trimmed_R_2P \
	> $ALIGNED_DIR/NGS0001.sam

# Convert SAM to BAM
cd $ALIGNED_DIR
samtools view -h -b NGS0001.sam > NGS0001.bam

# Sort the BAM
samtools sort NGS0001.bam -o NGS0001_sorted.bam

# Index the sorted BAM
samtools index NGS0001_sorted.bam

# Mark Duplicates
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

# Re-index the marked file
samtools index NGS0001_sorted_marked.bam

# Filter the marked file then index again
samtools view -b -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

# - Alignment statistics -
# 1. Flagstats
samtools flagstat NGS0001_sorted_filtered.bam > NGS0001_flagstats.txt

# 2. Idxstats
samtools idxstats NGS0001_sorted_filtered.bam > NGS0001_idxstats.txt

# 3. Insert size metrics
picard CollectInsertSizeMetrics \
I=NGS0001_sorted_filtered.bam \
O=insert_size_metrics.txt \
H=insert_size_histogram.pdf

# 4. Depth of coverage
bedtools coverage -a ~/ngs_assignment/dnaseq/data/annotation.bed -b NGS0001_sorted_filtered.bam > NGS0001_bedtools_coverage_results.txt

# --- Step 6: Variant Calling ---
# # 1. Decompress the genome
zcat ~/ngs_assignment/dnaseq/data/reference/hg19.fa.gz > ~/ngs_assignment/dnaseq/data/reference/hg19.fa

# 2. Index the genome
samtools faidx ~/ngs_assignment/dnaseq/data/reference/hg19.fa

# 3. Variant calling
freebayes --bam $ALIGNED_DIR/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_assignment/dnaseq/data/reference/hg19.fa > $RES_DIR/NGS0001.vcf

# 4. Quality filter the variants
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" $RES_DIR/NGS0001.vcf > $RES_DIR/NGS0001_quality_filtered.vcf

# 5. Spatial Filter
# INTERSECT VARIANTS WITH TARGET REGIONS 
# Note: Using 'head -n -1' to remove a truncated/corrupted line at the end 
# of the VCF caused by a disk-space/write error during the previous step. 
# Using 'stdin' to stream data and save disk space.
head -n -1 $RES_DIR/NGS0001_quality_filtered.vcf | bedtools intersect -header -wa -a stdin -b ~/ngs_assignment/dnaseq/data/annotation.bed > $RES_DIR/NGS0001_final_filtered.vcf

# 6. Compress and Index final filtered file
bgzip $RES_DIR/NGS0001_final_filtered.vcf
tabix -p vcf $RES_DIR/NGS0001_final_filtered.vcf.gz

# --- Step 7: Variant Annotation and Prioritisation ---
# Ensure the ANNOVAR source code is uploaded onto your openstack
# Unpack the uploaded file and et up annovar
tar -xvf ~/ngs_assignment/annovar.tar

# Download the databases needed for annotation
cd annovar
# Note: databases are summed to be downloaded during setup phase for this pipeline
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

# Convert VCF to Annovar input format
./convert2annovar.pl -format vcf4 $RES_DIR/NGS0001_final_filtered.vcf.gz > $RES_DIR/NGS0001_final_filtered.avinput

# Run Annovar table function
./table_annovar.pl $RES_DIR/NGS0001_final_filtered.avinput humandb/ -buildver hg19 -out $RES_DIR/NGS0001_annovar_final -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout


# --- snpEFF ---
# Note: I couldnt figure out how to solve the JNI error so I haven't included this in the workflow
# These following commands are what I aimed to use to download snpEFF, the hg19 ref database and execute further annotation
# download snpEFF: wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
# unzip snpEff_latest_core.zip
# download database: java -jar snpEff.jar download -v hg19
# Annotation commands: 
# cd ~/ngs_assignment
# java -jar ~/snpEff/snpEff.jar ann -v -lof -stats ~/ngs_assignment/dnaseq/results/NGS0001_snpeff_stats.html hg19 ~/ngs_assignment/dnaseq/results/NGS0001_final_filtered.vcf.gz > ~/ngs_assignment/dnaseq/results/NGS0001_snpeff_annotated.vcf
# Note: ann to annotate file, -v to show progress, -lof to add loss of function and nonsense mediated decay predictions, -stats to generate stats

# --- Basic Variant Prioritisation ---
# filter for exonic variants not seen in the dbSNP
# "." means there is no exac frequency, therefore the variant is not present in the population database
awk -F',' '$6 ~ /exonic/ && $21 == "."' $RES_DIR/NGS0001_annovar_final.hg19_multianno.csv > $RES_DIR/NGS0001_prioritized_variants.csv

#
#
#
#
#
#
#
#
#
#
