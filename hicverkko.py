#!/usr/bin/python3
import shutil
import sys
import os
import random
from os import listdir
from os.path import isfile, join
import evaluate_rukki
import cluster
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} <input_dir> <output_dir>')
        exit()
    #run clustering
    #python3 cluster.py unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches hic_mapping.byread.output > cluster.out 2> cluster.err

    #convert to rukki inputs
    #sh convert.sh

    # run rukki
    #sh rukki.sh
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    eval_file = ""
    if len(sys.argv) > 3:
        eval_file = os.path.join(sys.argv[2], sys.argv[3])

    os.makedirs(output_dir, exist_ok=True)
    #here should be mashmap running and parsing
    #TODO
    #mashmap -r unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.fasta -q unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.fasta -t 8 -f none --pi 95 -s 10000
    #cat mashmap.out |awk '{if ($NF > 99 && $4-$3 > 500000 && $1 != $6) print $1"\t"$6}'|sort |uniq > unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches
    matches_file = os.path.join(output_dir, "assembly.matches")
    hic_file = os.path.join(output_dir, "hic_mapping.byread.output")
    if not os.path.exists(hic_file):
        hic_file = os.path.join(output_dir, "hic.byread.compressed")
    compressed_hic = os.path.join(output_dir, "hic.byread.compressed")
    if os.path.exists(compressed_hic):
        hic_file = compressed_hic

    noseq_gfa = os.path.join(output_dir, "assembly.hpc.noseq.gfa")
    clustering_output = os.path.join(output_dir, "hicverkko.out")

    cluster.run_clustering(noseq_gfa, matches_file, hic_file, output_dir)
