#!/usr/bin/env python3


import sys
import re
import telomere_check
import os

data_file = "/data/antipovd2//devel/verkkohic/human_only.txt"
data_dir = "/data/antipovd2/data/"
data_dir = "/Users/antipovd2/work/verkkohic/data/"
res_dir = "/Users/antipovd2/work/verkkohic/res/"


def check_unclassified(rukki_res, t2t, gfa):
    lengths = {}
    total_len = 0
    for line in open(gfa, 'r'):
        arr = line.split()
        if arr[0] == "S":
            lengths[arr[1]] = int(arr[3].split(':')[-1])
            total_len += lengths[arr[1]]
    unalinged_count = 0
    unalinged_len = 0
    for line in open (rukki_res, 'r'):
        arr = line.split()
        if arr[0][0:3] == "na_":
            edges = arr[1].split(',')
            for e in edges:
                unalinged_count += 1
                tr_e = e[:-1]

                if tr_e in lengths:
                    unalinged_len += lengths[tr_e]
    print (f"Total count {len(lengths)}, unaligned edge fraction {unalinged_count/len(lengths)}")
    print (f"Total len {total_len}, unaligned length fraction {unalinged_len/total_len}")

for dir in os.listdir(data_dir):
    workdir = os.path.join(data_dir, dir)
    if os.path.isdir(workdir):
        print (f'dataset: {dir}')
        t2t = os.path.join(workdir, "assembly_graph", "assembly_graph.windows.0.4.50kb.ends.bed")
        if os.path.exists(t2t):
            rukki_res = os.path.join(res_dir, dir, sys.argv[1], "unitig-popped-unitig-normal-connected-tip.paths.tsv")
            gfa = os.path.join(workdir, "unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa")
            if os.path.exists(rukki_res):
                check_unclassified(rukki_res, t2t, gfa)


