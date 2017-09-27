for input in NA12878-chr22 NA12878-chr22-10pct SRR1611184 
do
if [ -f $input.only_mapped.reordered-contigs.sorted.deduplicated.read-group.sam ]
then rm $input.only_mapped.reordered-contigs.sorted.deduplicated.read-group.sam
fi
if [ -f $input.only_mapped.reordered-contigs.sorted.deduplicated.read-group.bam ]
then rm $input.only_mapped.reordered-contigs.sorted.deduplicated.read-group.bam
fi
if [ -f $input.only_mapped.reordered-contigs.read-group.sam ]
then rm $input.only_mapped.reordered-contigs.read-group.sam
fi
if [ -f $input.only_mapped.reordered-contigs.read-group.bam ]
then rm $input.only_mapped.reordered-contigs.read-group.bam
fi
done
