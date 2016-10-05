#! /bin/bash

#################################################
### Download SRA project and convert to fastq ###
#################################################

mkdir ../Reads
cd ../Reads
sh ../DataAquisition/download_from_sra.sh -f ../RawData/SraAccList.txt

# create list of sample names
python ../ProcessingCode/process_xml.py 

# rename files
python ../ProcessingCode/rename_files.py 

################################
### map reads and call peaks ###
################################

# Get the genome and build an index

mkdir ../genome
cd ../genome
wget
ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/arabidopsis_thaliana/dna/Arabidopsis_thalia\
na.TAIR10.dna.toplevel.fa.gz

# need to change chromosome names slightly
gzip -dc Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
| sed 's/>/>chr/g' - > tair10.fa
rm Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# make the chromosome sizes file for gem
samtools faidx tair10.fa
cut -f 1,2 tair10.fa.fai > tair10.chr.sizes
rm tair10.fa.fai

# Build index
# As we want to look for and report potentially many alignments,
# we set the offrate to 3, and this will speed up alignment
# but create a larger index
bowtie2-build --threads 10 --offrate 3 tair10.fa tair10
genome_path=`pwd`
cd ../Reads

# define a function to do the mapping
map() {
    # check if there are two fastq files
    if [ -f "${1}_2.fastq.gz"  ]; then
	nice -n 15 bowtie2 -k 10 -p 20 -R 5 --no-unal \
	     -q -x $2 -1 "${1}_1.fastq.gz" -2 "${1}_2.fastq.gz" \
	    | samtools view -b - > $1.bam
    else
	nice -n 15 bowtie2 -k 10 -p 20 -R 5 --no-unal \
	     -q -x $2 -U "${1}_1.fastq.gz" \
            | samtools view -b - > $1.bam
    fi
}

# loop over directories, moving each to the scratch drive to map
wd=`pwd` # remember where we are
for directory in ./*; do
    if [ -d $directory ]; then
	printf "Mapping ${directory}\n"
	mv $directory ~/scratch
	cd ~/scratch/$directory
	map $directory $genome_path/tair10
	current_time=`date`
	printf "\tFinished mapping ${directory} at ${current_time}\n"
	printf "Calling peaks\n"
	gem --f SAM --k_min 6 --kmax 20 --k_seqs 600 \
	    --k_neg_dinu_shuffle --t 10 --q 2 \
	    --d ~/working_data/Tools/gem/Read_Distribution_default.txt \
	    --g $genome_path/tair10.chr.sizes \
	    --expt $directory.bam
	cd ..
	mv $directory $wd
	cd $wd
    fi
done
