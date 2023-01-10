#!/usr/bin/env python3


import sys
import re
import telomere_check
import os

data_file = "/data/antipovd2//devel/verkkohic/human_only.txt"
data_dir = "/data/antipovd2/data/"
data_dir = "/Users/antipovd2/work/verkkohic/data/"
res_dir = "/Users/antipovd2/work/verkkohic/res/"


for dir in os.listdir(data_dir):
    workdir = os.path.join(data_dir, dir)
    if os.path.isdir(workdir):
        print (f'dataset: {dir}')
        t2t = os.path.join(workdir, "assembly_graph", "assembly_graph.windows.0.4.50kb.ends.bed")
        if os.path.exists(t2t):
            rukki_res = os.path.join(res_dir, dir, sys.argv[1], "unitig-popped-unitig-normal-connected-tip.paths.tsv")
            gfa = os.path.join(workdir, "unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa")
            if os.path.exists(rukki_res):
                telomere_check.evaluate_telomers(rukki_res, t2t, gfa)


