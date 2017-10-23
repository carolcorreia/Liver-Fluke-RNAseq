#####################################################
#      Cattle Liver Fluke RNA-seq Time Course       #
#####################################################

# Repository DOI badge: 
# Author: Andres Garcia Campos
# Last updated on: 18/10/2017

##################################################################
# Adapter-contamination and quality filtering of raw FASTQ files #
##################################################################

# Required software is ngsShoRT (version 2.2). More information can be found
# here: http://research.bioinformatics.udel.edu/genomics/ngsShoRT/index.html

# I have verified that trimming and FASTQ files worked with 2 of my samples. I proceed to do the samples with all samples
# Move to directory for filtered reads:
cd /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq

# Create bash scripts to perform filtering of each FASTQ file, keeping the
# sequencing lane information:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L001_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.001.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L002_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.002.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L003_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.003.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L004_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.004.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L005_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.005.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L006_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.006.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L007_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.007.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L008_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 5 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.008.sh; done;

# Change permission for execute and run all scripts on Stampede:

chmod 755 filtering.00*.sh

for script in `ls filtering.00*.sh`
do
nohup ./$script > ${script}.nohup &
done

# Scripts are running. There is a message that pops up:

## nohup: ignoring input and redirecting stderr to stdout

# Check that all files were processed:
for file in `ls filtering.00*.sh.nohup`; \
do grep -o 'Done-MAIN' $file | wc -l; done

# Compress files with discarded reads:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ \
-name extracted_*.txt`; do echo "gzip -9 $file" >> discarded_compression.sh; \
done;

# Split and run all scripts on Stampede:

split -d -l 70 discarded_compression.sh discarded_compression.sh.
chmod 755 discarded_compression.sh.*

for script in `ls discarded_compression.sh.*`; \
do
nohup ./$script &
done

# Gather ngsShoRT reports from all samples into one file:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ \
-name final_PE_report.txt`; \
do echo echo \
"\`dirname $file | perl -pe 's/.*(\A\w\d\d_S\d\d_L00\d)/\$1/'\` \
\`grep 'Read Pair Count:' $file\` \
\`grep 'Removed PE Pair\* Count:' $file\` >> \
ngsshort_cattle.txt" >> ngsshort_summary_cattle.sh
done

chmod 755 ngsshort_summary_cattle.sh
./ngsshort_summary_cattle.sh

# Pending to add header 

# Sort samples in ngsshort_cattle.txt
sort ngsshort_cattle.txt >>sorted_ngsshort_summary.txt


# Transfer ngsShoRT summary to laptop via SCP.
scp -r acampos@rodeo.ucd.ie://home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/sorted_ngsshort_summary.txt \
Dropbox/Liver-Fluke-RNAseq/ngsshort/


################################################
# FastQC quality check of filtered FASTQ files #
################################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering
cd /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering

# Create a bash script to perform FastQC quality check on all filtered
# FASTQ files:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ \
-name *_L00*_R*.fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering $file" \
>> fastqc_filt.sh; done;

# Split and run all scripts on Stampede:
split -d -l 100 fastqc_filt.sh fastqc_filt.sh.
chmod 755 fastqc_filt.sh.*

for script in `ls fastqc_filt.sh.*`; \
do nohup ./$script > ${script}.nohup & 
done

# Check if all the files were processed:
for file in `ls fastqc_filt.sh.*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering/tmp; \
done;

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering/tmp \
-name summary.txt`; do more $file >> reports_post-filtering.txt; done

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post-filtering.txt; \
done

# Remove temporary folder and its files:
rm -rf /home/workspace/acampos/flukeRNAseq/quality_check/post-filtering/tmp 

# Copy files in laptop
scp -r acampos@rodeo.ucd.ie://home/workspace/acampos/flukeRNAseq/quality_check/post-filtering .

