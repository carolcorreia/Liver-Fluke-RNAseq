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
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.001.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L002_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.002.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L003_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.003.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L004_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.004.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L005_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.005.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L006_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.006.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L007_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
-pe1 $file -pe2 $file2 \
-o /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/$sample \
-methods 5adpt_lqr -5a_f i-p -5a_mp 90 \
-5a_del 0 -5a_ins 0 -5a_fmi 100 -5a_axn kr -lqs 20 -lq_p 25 -gzip" \
>> filtering.007.sh; done;

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data/ \
-name *L008_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/R1(_001.fastq.gz)$/R2$1/'`; \
sample=`basename $file | perl -p -e 's/_R1_001.fastq.gz$//'`; \
echo "perl /usr/local/src/ngsShoRT_2.2/ngsShoRT.pl -t 20 -mode trim -min_rl 100 \
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


