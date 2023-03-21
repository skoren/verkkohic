#!/bin/bash

set -e -o pipefail

echo "Aligning reads from $3..."
if [ -e $3/hic ]; then
   HIC1=`ls $3/hic/*R1_001.fast[aq].gz|tr '\n' ' ' |awk '{print substr($0, 1, length($0)-1)}'`
   HIC2=`ls $3/hic/*R2_001.fast[aq].gz|tr '\n' ' ' |awk '{print substr($0, 1, length($0)-1)}'`
   echo "Hic pairs 1: '$HIC1'"
   echo "Hic pairs 2: '$HIC2'"

   if [ ! -s $2/hic_to_assembly.sorted_by_read.bam ]; then
      echo "Mapping Illumina HiC w/BWA"
      bwa index $2/unitigs.fasta && bwa mem -t $4 -5 -S -P $2/unitigs.fasta <(zcat $HIC1) <(zcat $HIC2) | samtools view -bh -@ $4 -q 1 - | samtools sort -n -@ $4 - -o $2/hic_to_assembly.sorted_by_read.bam
   fi
elif [ -e $3/porec ]; then
   POREC=`ls $3/porec/*.fast[aq].gz |tr '\n' ' ' |awk '{print substr($0, 1, length($0)-1)}'`
   echo "Porec file: $POREC"

   if [ ! -s $2/hic_to_assembly.sorted_by_read.bam ]; then
      echo "Mapping PoreC w/Minimap"
      minimap2 -a -x map-ont -k 17 -t $4 -K 10g -I 8g $2/unitigs.fasta <(zcat $POREC) | samtools view -bh -@ $4 -q 1 - > $2/hic_to_assembly.sorted_by_read.bam
   fi
else
   echo "Error no reads found, expected hic or porec folder"
   exit -1
fi
