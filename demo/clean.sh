#!/bin/bash
for input in NA12878.exome.hg19 NA12878.exome.hg19.chr22 
do
if [ -f $input.sfm.out.bam ]
then rm $input.sfm.out.bam
fi
if [ -f $input.filter.out.bam ]
then rm $input.filter.out.bam
fi
if [ -f $input.sfm.out.metrics ]
then rm $input.sfm.out.metrics
fi
if [ -f $input.filter.out.metrics ]
then rm $input.filter.out.metrics
fi
if [ -f $input.sfm.out.recal ]
then rm $input.sfm.out.recal
fi
if [ -f $input.filter.out.recal ]
then rm $input.filter.out.recal
fi
done
