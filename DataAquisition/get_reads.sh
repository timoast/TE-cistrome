#! /bin/bash

# Download SRA project and convert to fastq
cd ../RawData/
sh ../DataAquisition/download_from_sra.sh -f SraAccList.txt

# Map each fastq file to the genome allowing lots of multimappers
map() {
    bowtie2 
}

for directory in $(ls -d .); do
    if [ -d $directory ]; then
	cd $directory
	map
	cd ..
    fi
done

