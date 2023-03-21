#!/bin/bash

set -e -o pipefail

# this script converts an unphased verkko assembly and HPC-compressed gfa to an uncompressed assembly and GFA
if [ ! -s $2/unitigs.gfa ]; then
   echo "Generating uncompressed gfa..."

   cp $1/assembly.homopolymer-compressed.noseq.gfa $2/unitigs.homopolymer-compressed.noseq.gfa

   # create the appropriate input for mapping
   ## create renaming map...
   cat $1/6-*/*scfmap \
       | grep utig4 \
       | awk '{print $2"\t"$NF}' \
       > $2/unitigs.rename.map

   ## rename fasta from ">unassigned-0000001" to names like ">utig4-0" (etc)
   python $VERKKO/lib/verkko/scripts/fasta_combine.py rename $2/unitigs.fasta $2/unitigs.rename.map $1/assembly.fasta

   ## align to get homopolymer uncompressed GFA?
   $VERKKO/lib/verkko/bin/alignGFA \
       -V -e 0.30 \
       -gfa \
       -i $1/assembly.homopolymer-compressed.gfa \
       -T $2/unitigs.fasta 0 \
       -t $3 \
       -o $2/unitigs.gfa
fi
