#####################################################
#      Cattle Liver Fluke RNA-seq Time Course       #
#####################################################

# Repository DOI badge: 
# Author: Carolina N. Correia
# Last updated on: 22/08/2017

###########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering
cd !$

# Run FastQC in one file to check if it works:
fastqc -o /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering \
--noextract --nogroup -t 10 \
/home/workspace/acampos/flukeRNAseq/fastq/raw_data/A01_S23_L001_R1_001.fastq.gz

# Transfer compressed folder to personal laptop via SCP
# and check the HTML report:
scp ccorreia@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/A01_S23_L001_R1_001_fastqc.zip .

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data \
-name *.fastq.gz`; \
do echo "fastqc --noextract --nogroup -t 10 \
-o /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering $file" \
>> fastqc.sh; \
done

# Check number of lines in script:
wc -l fastqc.sh

# Split and run all scripts on Stampede:
split -d -l 264 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Analysis complete" >> succesful_fastqc.txt
done

wc -l succesful_fastqc.txt

# Deleted all HTML files:
rm -r *.html


# Collect FastQC stats:
mkdir /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp; \
done

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp \
-name summary.txt`; do more $file >> summary_pre-filtering.txt; \
done

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Transfer compressed folders to personal laptop via SCP
# and check HTML reports:
scp -r \
ccorreia@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/quality_check/pre-filtering/tmp .

# Remove temporary folder and its files:
rm -r tmp

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# Create a working directory for filtered reads:
mkdir /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq
cd !$

# I put file with adapters in this directory from my terminal (not from rodeo)
scp -r scp -r Adapter_sequence.txt acampos@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/

# Run ngsShoRT in one pair of reads to check if it's working:
perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 /home/workspace/acampos/flukeRNAseq/fastq/raw_data/A01_S23_L001_R1_001.fastq.gz \
-pe2 /home/workspace/acampos/flukeRNAseq/fastq/raw_data/A01_S23_L001_R2_001.fastq.gz \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L001 \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Compress files with discarded reads:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ \
-name extracted_*.txt`; do echo "gzip -9 $file" >> discarded_compression.sh; \
done;

# Split and run all scripts on Stampede:
split -d -l 70 discarded_compression.sh discarded_compression.sh.
for script in `ls discarded_compression.sh.*`
do
chmod 755 $script
nohup ./$script &
done

# Gather ngsShoRT reports from sample into one file:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ \
-name final_PE_report.txt`; \
do echo echo \
"\`dirname $file | perl -pe 's/.*(A\d\d\d\d.*\d)/\$1/'\` \
\`grep 'Read Pair Count:' $file\` \
\`grep 'Removed PE Pair\* Count:' $file\` >> \
ngsshort_cattle_A01_S23_L001.txt" >> ngsshort_summary_cattle_A01_S23_L001.sh
done

chmod 755 ngsshort_summary_cattle_A01_S23_L001.sh
./ngsshort_summary_cattle_A01_S23_L001.sh

# Because I have only 1 sample, I do not add header.

# scp gsshort_summary_cattle_A01_S23_L001.txt to share folder

################################################
# FastQC quality check of filtered FASTQ files #
################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock
cd /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock

# Run FastQC in one file to see if it's working well:
fastqc -o /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock --noextract --nogroup -t 10 /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L001/trimmed_A01_S23_L001_R1_001.fastq.gz

# Secure copy to my laptop
scp -r acampos@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock/trimmed_A01_S23_L001_R1_001_fastqc.zip .

# unzip content of file
unzip trimmed_A01_S23_L001_R1_001_fastqc.zip

# with do with R2
fastqc -o /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock --noextract --nogroup -t 10 /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L001/trimmed_A01_S23_L001_R2_001.fastq.gz

# Secure copy to share folder
scp -r acampos@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock/trimmed_A01_S23_L001_R2_001_fastqc.zip .

# Run ngsShoRT in one pair of reads to check if it's working:
perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 /home/workspace/acampos/flukeRNAseq/fastq/raw_data/A01_S23_L002_R1_001.fastq.gz \
-pe2 /home/workspace/acampos/flukeRNAseq/fastq/raw_data/A01_S23_L002_R2_001.fastq.gz \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L002 \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &

# Compress files with discarded reads:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ \
-name extracted_*.txt`; do echo "gzip -9 $file" >> discarded_compression.sh; \
done;

# Split and run all scripts on Stampede:
split -d -l 70 discarded_compression.sh discarded_compression.sh.
for script in `ls discarded_compression.sh.*`
do
chmod 755 $script
nohup ./$script &
done

# Gather ngsShoRT reports from sample into one file:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ \
-name final_PE_report.txt`; \
do echo echo \
"\`dirname $file | perl -pe 's/.*(A\d\d\d\d.*\d)/\$1/'\` \
\`grep 'Read Pair Count:' $file\` \
\`grep 'Removed PE Pair\* Count:' $file\` >> \
ngsshort_cattle_A01_S23_L002.txt" >> ngsshort_summary_cattle_A01_S23_L002.sh
done

chmod 755 ngsshort_summary_cattle_A01_S23_L002.sh
./ngsshort_summary_cattle_A01_S23_L002.sh

# Go to directory for FASTQC
cd /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock

# Run FastQC in one file to see if it's working well:

fastqc -o /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock --noextract --nogroup -t 10 /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L002/trimmed_A01_S23_L002_R1_001.fastq.gz
fastqc -o /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock --noextract --nogroup -t 10 /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L002/trimmed_A01_S23_L002_R2_001.fastq.gz

# Copy file on share folder and unzip
### and checked the HTML report. It worked fine.
