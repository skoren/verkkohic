#!/bin/bash

# build mappings
# this used the homoplymer compressed contigs but could probably use regular  as well
mashmap -r unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.fasta -q unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.fasta -t 8 -f none --pi 95 -s 10000
cat mashmap.out |awk '{if ($NF > 99 && $4-$3 > 500000 && $1 != $6) print $1"\t"$6}'|sort |uniq > unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches

#run clustering
python3 cluster.py unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches hic_mapping.byread.output > cluster.out 2> cluster.err

#convert to rukki inputs 
sh convert.sh 

# run rukki
sh rukki.sh
