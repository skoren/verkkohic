#!/usr/bin/python3

import sys
import os
import random
from os import listdir
from os.path import isfile, join

if len(sys.argv) < 3:
    print(f'Usage: {sys.argv[0]} <input_dir> <output_dir>')
    exit()
# build mappings
# this used the homoplymer compressed contigs but could probably use regular  as well

#run clustering
#python3 cluster.py unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches hic_mapping.byread.output > cluster.out 2> cluster.err

#convert to rukki inputs
#sh convert.sh

# run rukki
#sh rukki.sh
cur_dir = os.path.abspath(os.path.dirname(__file__))
input_dir = sys.argv[1]
output_dir = sys.argv[2]
#here should be mashmap running and parsing
#mashmap -r unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.fasta -q unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.fasta -t 8 -f none --pi 95 -s 10000
#cat mashmap.out |awk '{if ($NF > 99 && $4-$3 > 500000 && $1 != $6) print $1"\t"$6}'|sort |uniq > unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches
matches_file = os.path.join(input_dir, "unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches")
hic_file = os.path.join(input_dir, "hic_mapping.byread.output")
compressed_hic = os.path.join(output_dir, "hic.byread.compressed")
if os.path.exists(compressed_hic):
    hic_file = compressed_hic

noseq_gfa = os.path.join(input_dir, "unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa")
clustering_output = os.path.join(output_dir, "cluster.out")
os.system(f'python3 cluster.py {noseq_gfa} {matches_file} {hic_file} > {clustering_output}')
csv_output = os.path.join(output_dir, "unitig-popped-unitig-normal-connected-tip.colors.csv")
#Parsing clustering output
#echo -e "node\tmat\tpat\tmat:pat\tcolor" > unitig-popped-unitig-normal-connected-tip.colors.csv
csv_file = open(csv_output, 'w')
csv_file.write("node\tmat\tpat\tmat:pat\tcolor\n")

contig_names = set()
for line in open(noseq_gfa, 'r'):
    arr = line.split()
    if arr[0] == "S":
        contig_names.add(arr[1])

for line in open (clustering_output, 'r'):
    if len(line) > 3 and line[0:3] == "RES":
        right = False
        for lines in line.split('}, {'):
            arr = lines.split('\'')
            for contig in arr:
                if contig in contig_names:
                    if right:
                        csv_file.write(f'{contig}\t0\t100000\t0:100000\t#8888FF\n')
                    else:
                        csv_file.write(f'{contig}\t100000\t0\t100000:0\t#FF8888\n')
            right = True
csv_file.close()

#cat cluster.out |grep -A 1 Seed|grep -v Initial | grep -v Seed |awk -F "}," '{alen=split($1, a, ","); blen=split($2, b, ","); for (i = 1; i<=alen; i++) { print a[i]"\t0\t100000\t0:100000\t#8888FF"} for (i = 1; i<= blen; i++) {print b[i]"\t100000\t0\t100000:0\t#FF8888"} }'|sed 's/({//g' |sed 's/})//g' |sed s/\'//g|sed s/\ //g |sed 's/{//g' | sort |uniq |grep -w -v -f unassigned >> unitig-popped-unitig-normal-connected-tip.colors.csv
#hi-c gfa(noseq) trio_colors


#../../devel/rukki/target/release/rukki trio --graph unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa --markers unitig-popped-unitig-normal-connected-tip.colors.csv -p out.path.tsv
rukki_line = f'rukki trio --graph {noseq_gfa} --markers {csv_output} -p {os.path.join(output_dir, "unitig-popped-unitig-normal-connected-tip.paths.tsv")}'
#what about rukki options?
os.system(rukki_line)
