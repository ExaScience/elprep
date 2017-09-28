input=NA12878_S1 #without extension
threads=`grep -c processor /proc/cpuinfo`

if [ ! -f $input.bam ]
then
    echo "Input file not present."
    if command -v curl > /dev/null
    then
        echo "Please download the input file. See http://www.ebi.ac.uk/ena/data/view/ERS189474"
        exit
    fi
fi

if command -v elprep > /dev/null
then
    if command -v samtools > /dev/null
    then
        elprep-sfm.py $input.bam $input.only_mapped.reordered-contigs.sorted.deduplicated.read-group.bam --filter-unmapped-reads --replace-reference-sequences ucsc.hg19.dict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --intermediate-files-output-type bam --nr-of-threads $threads
    else
        echo 'samtools not found. Please download it from http://samtools.sourceforge.net and make sure that its binary is present in your PATH.'
    fi
else
    echo 'elprep not found. Please download it from https://github.com/exascience/elprep and make sure that its binary is present in your PATH.'
fi
