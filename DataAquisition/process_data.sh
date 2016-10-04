#! /bin/bash

# Download SRA project and convert to fastq
mkdir ../Reads
cd ../Reads
sh ../DataAquisition/download_from_sra.sh -f SraAccList.txt

# create list of sample names
python ../ProcessingCode/process_xml.py 

# rename files
python ../ProcessingCode/rename_files.py 

# map reads
map() {
    # check if there are two fastq files
    if [ -f "${1}_2.fastq.gz"  ]; then
	nice -n 15 bowtie2 -k 10 -p 20 -R 5 --rg-id $1 \
	     -q -x tair10 -1 "${1}_1.fastq.gz" -2 "${1}_2.fastq.gz" \
	    | samtools view -b - > "${1}.bam"
    else
	nice -n 15 bowtie2 -k 10 -p 20 -R 5 --rg-id $1 \
	     -q -x tair10 -U "${1}_1.fastq.gz" \
            | samtools view -b - > "${1}.bam"
    fi
}


for directory in ./*; do
    if [ -d $directory ]; then
	cd $directory
	printf "Mapping ${directory}\n"
	map $directory
	current_time="`date`"
	printf "\tFinished mapping ${directory} at ${current_time}\n"
	cd ..
    fi
done

# Call peaks
for directory in ./*; do
    if [ -d $directory ]; then
	cd $directory
	for bamfile in $(ls *.bam); do
	    gem --f SAM --k_min 6 --kmax 20 --k_seqs 600 --k_neg_dinu_shuffle
	done
    fi
done
