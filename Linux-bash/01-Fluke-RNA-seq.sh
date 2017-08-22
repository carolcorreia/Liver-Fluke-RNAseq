#####################################################
#      Cattle Liver Fluke RNA-seq Time Course       #
#####################################################

# Repository DOI badge: 
# Author: Carolina N. Correia
# Last updated on: 22/08/2017

##############################
# Download raw data from BYU #
##############################

# Create and enter working directory on Rodeo:
mkdir /home/workspace/acampos/flukeRNAseq/fastq
cd !$

# Download files:
screen -D -R fluke_download
scp -r fslcollab161@scp.fsl.byu.edu:/fslgroup/fslg_dnasc/compute/170803_D00723_0206_ACBDG4ANXX/Unaligned-8lanes/Project1 .

# Detach the screen session by pressing Ctrl+A then d.

# Terminate the ssh session
exit

# Re-attach the screen session:
screen -D -R fluke_download

# When the download is complete, rename the directory:
mv Project1 raw_data

# Modify file permissions to read and execute only:
chmod 555 all.md5

for file in `find /home/workspace/acampos/flukeRNAseq/fastq/raw_data \
-name *.fastq.gz`; \
do chmod 555 $file; \
done

########################
# Perform MD5 checksum #
########################

# Enter working directory:
cd /home/workspace/acampos/flukeRNAseq/fastq/raw_data

# Perform md5sum check:
md5sum -c all.md5 >> \
/home/workspace/acampos/flukeRNAseq/fastq/raw_data/md5check_UCD.txt

# Check that all files passed the check:
grep -c 'OK' md5check_UCD.txt

# Modify folder permissions to read and execute only:
cd /home/workspace/acampos/flukeRNAseq/fastq
chmod -R 555 raw_data

# File names correspond to:
A01_S23_L001_R1_001.fastq.gz


###############################
# Quality check and filtering #
###############################

# Please go to file: 02-Fluke-RNA-seq.sh
# for subsequent analysis.

