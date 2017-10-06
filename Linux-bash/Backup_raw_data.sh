##############################################
##    BACKUP DATA  in external harddriver   ##
##############################################

#The purpose is to create a secure copy of the raw data.

# Create and enter working directory on EXTERNAL HARD DRIVER:
mkdir /Andres/Volumes/Seagate_Backup_Plus_Drive/backup_RNA
cd !$

# I make a try with the folder created in the afternoon on Rodeo 
scp -r acampos@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/fastq/filt_fastq/ .

# Passphrase for key is requested
# Verify that folder is in the new location
ls -l

# Download files:
scp -r acampos@rodeo.ucd.ie:/home/workspace/acampos/flukeRNAseq/fastq/raw_data/ .

# Write passphrase for key again

# Once files are copied, terminate the ssh session

# check file permissions (only read and execute)
ls -l 