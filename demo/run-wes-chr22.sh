#!/bin/bash
input=NA12878.exome.hg19.chr22 #without extension
ref=hg19
sites1=Mills_and_1000G_gold_standard.indels.hg19
sites2=dbsnp_137.hg19
bed=nexterarapidcapture_expandedexome_targetedregions

getdata() {
    if [ ! -f $1 ]
    then
        echo "Input file $1 not present."
        if command -v curl > /dev/null
        then
	    echo "Downloading input file $1."
	    curl -O http://www.exascience.com/public-files/elprep-demo/$1
        elif command -v wget > /dev/null
        then
	    echo "Downloading input file $1."
	    wget http://www.exascience.com/public-files/elprep-demo/$1
        else
	    echo "Please download the input file from http://www.exascience.com/public-files/elprep-demo/$1 before running this script."
	    exit
        fi
   fi
}

getdata $input.bam
getdata $ref.elfasta
getdata $sites1.elsites
getdata $sites2.elsites
getdata $bed.bed

if command -v elprep > /dev/null
then
	    elprep filter $input.bam $input.filter.out.bam --mark-duplicates --mark-optical-duplicates $input.filter.out.metrics --sorting-order coordinate --bqsr $input.filter.out.recal --known-sites $sites1.elsites,$sites2.elsites --bqsr-reference $ref.elfasta --filter-non-overlapping-reads $bed.bed
else
    echo 'elprep not found. Please download it from https://github.com/exascience/elprep and make sure that its binary is present in your PATH.'
fi
