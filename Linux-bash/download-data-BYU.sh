# Downloading RNA-seq data from BYU to Rodeo
# Liver fluke project
# Author: Carolina N. Correia
# 21/08/2017

# Create and enter working directory on Rodeo:
mkdir /home/workspace/acampos/flukeRNAseq/raw_data
cd !$

# Download files:
screen -D -R fluke_download
scp -r username@remoteserver:remote_location .

# Detach the screen session by pressing Ctrl+A then d.

# Terminate the ssh session
exit

# Re-attach the screen session:
screen -D -R fluke_download

# Change directory permissions to read and execute only:
chmod -R 555 raw_data

# Create md5sum working folder:
mkdir /home/workspace/acampos/flukeRNAseq/raw_data/md5_check
cd !$

# Create bash script to perform md5sum check:
for file in `find /home/workspace/acampos/flukeRNAseq/raw_data -name md5.txt`; \
do echo "cd `dirname $file` && md5sum -c `basename $file` >> \
$HOME/scratch/miRNAseqValidation/md5check/md5_UCD.txt" >> md5sum.sh; \
done

# Run script on Stampede:
chmod 755 ./md5sum.sh
nohup ./md5sum.sh &

# Check that all files passed the check:
grep -c 'OK' md5_UCD.txt