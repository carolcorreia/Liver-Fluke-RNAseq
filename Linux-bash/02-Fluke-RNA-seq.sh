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
-methods 5adpt_lqr -5a_f Adapter_sequence.txt -5a_mp 90 -5a_del 0 \
-5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip &


# file A01_S23_L001 is created. However, the following message pops up in rodeo: "ERROR \
# (extract_5adpt_sequences_to_arrays_and_print_header): the lines listed in the \
# 5' adapter file, Adapter_sequence.txt, do not seem to follow the standard format.\
# Please use five_prime_adapter_seq_TEMPLATE.txt as your template."

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
# Because all this is a one-sample trial, I make a directory to keep all files and directories created separately from the future analysis
mkdir /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/mock_one_sample

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
fastqc -o /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering-mock \
--noextract --nogroup -t 10 \
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L001/trimmed_A01_S23_L001_R1_001.fastq.gz

#I have an error message:
# Failed to process /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L001/trimmed_A01_S23_L001_R1_001.fastq.gz
# java.io.EOFException
#	at java.util.zip.GZIPInputStream.readUByte(GZIPInputStream.java:268)
#	at java.util.zip.GZIPInputStream.readUShort(GZIPInputStream.java:258)
#	at java.util.zip.GZIPInputStream.readHeader(GZIPInputStream.java:164)
#	at java.util.zip.GZIPInputStream.<init>(GZIPInputStream.java:79)
#	at java.util.zip.GZIPInputStream.<init>(GZIPInputStream.java:91)
#	at uk.ac.babraham.FastQC.Utilities.MultiMemberGZIPInputStream.<init>(MultiMemberGZIPInputStream.java:37)
#	at uk.ac.babraham.FastQC.Sequence.FastQFile.<init>(FastQFile.java:80)
#	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:106)
#	at uk.ac.babraham.FastQC.Sequence.SequenceFactory.getSequenceFile(SequenceFactory.java:62)
#	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.processFile(OfflineRunner.java:129)
#	at uk.ac.babraham.FastQC.Analysis.OfflineRunner.<init>(OfflineRunner.java:102)
#	at uk.ac.babraham.FastQC.FastQCApplication.main(FastQCApplication.java:316)


