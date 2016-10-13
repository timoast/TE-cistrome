#! /bin/bash

cd ../Reads

for directory in ./*; do
    if [ -d $directory ]; then
	cd $directory
	gemfile="${directory}_GEM_events.txt"
	python ../../ProcessingCode/reformat_gem.py $gemfile \
	       > "${directory}_GEM_events.bed"
	cd ..
    fi
done
