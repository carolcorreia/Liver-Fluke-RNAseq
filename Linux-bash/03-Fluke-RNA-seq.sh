#####################################################
#      Cattle Liver Fluke RNA-seq Time Course       #
#####################################################

# Repository DOI badge: 
# Author: Andres Garcia Campos
# Last updated on: 09/11/2017

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

# Copy files in personal laptop in shared folder
scp -r acampos@rodeo.ucd.ie://home/workspace/acampos/flukeRNAseq/quality_check/post-filtering .

##############################################################################
# Alignment of FASTQ files against the Bos taurus reference genome with STAR #
##############################################################################

# Required software is STAR 2.5.1b, consult manual/tutorial for details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# NCBI changed the location of several genome and annotation files in their FTP
# server on December 2016, which means that the original links listed below
# do not work anymore. The same files used in this analysis can now be found at:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1

# Download Bos taurus reference genome, version UMD3.1.1 from NCBI:  Folder last version: 11/06/2016
mkdir /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file
cd /workspace/storage/genomes/bostaurus/UMD3.1.1_NCBI/source_file
nohup wget nohup wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz &

gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz

# Download annotation file for UMD3.1.1 NCBI Bos taurus genomic annotation (GCF_000003055.6): Last version: 21/09/2017

cd /home/workspace/acampos/flukeRNAseq/bostaurus/UMD3.1.1_NCBI/annotation_file_Sept_2017
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz

# Generate genome indexes files using annotations:
mkdir /workspace/acampos/flukeRNAseq/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index
cd /workspace/acampos/flukeRNAseq/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index

nohup STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_124 \
--genomeFastaFiles \
/home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/source_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna \
--sjdbGTFfile /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 124 \
--outFileNamePrefix \
/home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_124 &

# Create and enter alignment working directory:
mkdir /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment
cd /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment

# Mapping reads from one FASTQ file to the indexed genome,
# to check if it works well:
nohup STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_124/ \
--readFilesIn \
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L001/trimmed_A01_S23_L001_R1_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L002/trimmed_A01_S23_L002_R1_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L003/trimmed_A01_S23_L003_R1_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L004/trimmed_A01_S23_L004_R1_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L005/trimmed_A01_S23_L005_R1_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L006/trimmed_A01_S23_L006_R1_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L007/trimmed_A01_S23_L007_R1_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L008/trimmed_A01_S23_L008_R1_001.fastq.gz \
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L001/trimmed_A01_S23_L001_R2_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L002/trimmed_A01_S23_L002_R2_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L003/trimmed_A01_S23_L003_R2_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L004/trimmed_A01_S23_L004_R2_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L005/trimmed_A01_S23_L005_R2_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L006/trimmed_A01_S23_L006_R2_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L007/trimmed_A01_S23_L007_R2_001.fastq.gz,\
/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/A01_S23_L008/trimmed_A01_S23_L008_R2_001.fastq.gz \
--readFilesCommand gunzip -c --outFilterMultimapNmax 10 \
--outFilterMismatchNmax 10 --outFileNamePrefix ./A01_S23\
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx &

# Put files of the one-sample try run in one new folder

mkdir mkdir /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/mock_one_sample
cd mkdir /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/mock_one_sample

mv *--outSAMtype* mock_one_sample/
mv nohup.out mock_one_sample/

# Create a bash script to perform alignment of paired FASTQ files:
for file in `find /home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/*/ \
-name *_L001_R1_001.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/L001/L002/g'`; \
file3=`echo $file | perl -p -e 's/L001/L003/g'`; \
file4=`echo $file | perl -p -e 's/L001/L004/g'`; \
file5=`echo $file | perl -p -e 's/L001/L005/g'`; \
file6=`echo $file | perl -p -e 's/L001/L006/g'`; \
file7=`echo $file | perl -p -e 's/L001/L007/g'`; \
file8=`echo $file | perl -p -e 's/L001/L008/g'`; \
read1=`echo $file | perl -p -e 's/R1/R2/'`; \
read2=`echo $file2 | perl -p -e 's/R1/R2/'`; \
read3=`echo $file3 | perl -p -e 's/R1/R2/'`; \
read4=`echo $file4 | perl -p -e 's/R1/R2/'`; \
read5=`echo $file5 | perl -p -e 's/R1/R2/'`; \
read6=`echo $file6 | perl -p -e 's/R1/R2/'`; \
read7=`echo $file7 | perl -p -e 's/R1/R2/'`; \
read8=`echo $file8 | perl -p -e 's/R1/R2/'`; \
sample=`basename $file | perl -p -e 's/_L001_R1_001.fastq.gz//'`; \
foldername=`basename $sample | perl -p -e 's/trimmed\_//'`; \
echo "mkdir /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/$foldername; \
cd /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/$foldername; \
STAR --runMode alignReads --runThreadN 1 --genomeLoad LoadAndRemove \
--genomeDir /home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/STAR-2.5.2b_index_124/ \
--readFilesIn $file,$file2,$file3,$file4,$file5,$file6,$file7,$file8 \
$read1,$read2,$read3,$read4,$read5,$read6,$read7,$read8 --readFilesCommand gunzip -c \
--outFilterMultimapNmax 10 --outFilterMismatchNmax 10 \
--outFileNamePrefix ./${foldername}_ --outSAMtype BAM Unsorted \
--outSAMattrIHstart 0 --outSAMattributes Standard --outReadsUnmapped Fastx" \
>> alignment.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 22 alignment.sh alignment.sh.
chmod 755 alignment.sh.*

for script in `ls alignment.sh.*`; \
do
nohup ./$script > ${script}.nohup &
done

# Check nohup.out file to see how many jobs finished successfully:
grep -c 'finished successfully' alignment.sh.00.nohup
# 22: It worked
grep -c 'finished successfully' alignment.sh.01.nohup
# 22: It worked
grep -c 'finished successfully' alignment.sh.02.nohup
# 22: It worked

# Merge all STAR log.final.out files into a single file:
for file in `find /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/*_S*/ \
-name *Log.final.out`; \
do perl /home/workspace/acampos/flukeRNAseq/star_report_opener.pl -report $file; done;

#############################################
# FastQC quality check of aligned BAM files #
#############################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and go to working directory:
mkdir /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment
cd /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment

# I try one sample first:

fastqc --noextract --nogroup -t 2 -o /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment \
/home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/A01_S23/A01_S23_Aligned.out.bam

# Analysis is complete. Transfer to own laptop and check
mkdir /Users/andresgarciacampos/Dropbox/Liver-Fluke-RNAseq/Quality_check/post_alignment_mock
scp -r acampos@rodeo.ucd.ie://home/workspace/acampos/flukeRNAseq/quality_check/post_alignment/A01_S23_Aligned.out_fastqc.zip .

unzip A01_S23_Aligned.out_fastqc.zip

# File looks ok. Then I make script ru run all files

# Create a bash script to perform FastQC quality check on aligned SAM files:
for file in `find /home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/*_S*/ \
-name *.bam`; do echo "fastqc --noextract --nogroup -t 2 \
-o /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment $file" >> \
fastqc_aligned.sh; done;

# Split and run all scripts on Stampede
split -d -l 22 fastqc_aligned.sh fastqc_aligned.sh.
chmod 755 fastqc_aligned.sh.*

for script in `ls fastqc_aligned.sh.*`;
do
nohup ./$script > ${script}.nohup &
done

# Delete all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment/tmp; \
done

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment/tmp/*_fastqc/ \
-name summary.txt`; do more $file >> reports_post-alignment.txt; done

for file in \
`find /home/workspace/acampos/flukeRNAseq/quality_check/post_alignment/tmp/*_fastqc/ \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_post_alignment.txt; \
done

# Check if all files were processed:
grep -c '##FastQC' basic_stats_post_alignment.txt
#66
grep -c 'Basic Statistics' reports_post-alignment.txt
#66
grep -c 'Analysis complete' fastqc_aligned.sh.00.nohup
#22
grep -c 'Analysis complete' fastqc_aligned.sh.01.nohup
#22
grep -c 'Analysis complete' fastqc_aligned.sh.02.nohup
#22

# Secure copy tmp folder to laptop (share folder) and rename folder as post_alignment.
# I use my own laptop (out of Rodeo)
cd /Users/andresgarciacampos/Dropbox/Liver-Fluke-RNAseq/Quality_check
scp -r acampos@rodeo.ucd.ie://home/workspace/acampos/flukeRNAseq/quality_check/post_alignment/tmp .

# Remove temporary folder in Rodeo:
rm -r tmp/

###################################################################
# Summarisation of gene counts with featureCounts for sense genes #
###################################################################

# Required package is featureCounts, which is part of Subread 1.5.1 software,
# consult manual for details:
# http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf

# Create working directories:
cd /home/workspace/acampos/flukeRNAseq
mkdir -p Count_summarisation/sense
cd /home/workspace/acampos/flukeRNAseq/Count_summarisation/sense

###### Pending to do ###########

# Run featureCounts with one sample to check if it is working fine:
featureCounts -a \
/home/workspace/genomes/bostaurus/UMD3.1.1_NCBI/annotation_file/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff \
-B -p -C -R -s 1 -T 15 -t gene -g Dbxref -o ./counts.txt \
/home/workspace/acampos/flukeRNAseq/STAR-2.5.2b_alignment/A01_S23/A01_S23_Aligned.out.bam


