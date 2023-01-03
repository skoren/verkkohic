#!/usr/bin/env python3


import sys
import re
import telomere_check
import os

data_file = "/data/antipovd2//devel/verkkohic/human_only.txt"
data_dir = "/data/antipovd2/data/"
data_dir = "/Users/antipovd2/work/verkkohic/data/"

for dir in os.listdir(data_dir):
    workdir = os.path.join(data_dir, dir)
    if os.path.isdir(workdir):
        print (f'dataset: {dir}')
        t2t = os.path.join(workdir, "assembly_graph", "assembly_graph.windows.0.4.50kb.ends.bed")
        if os.path.exists(t2t):
            telomere_check.evaluate_telomers(os.path.join(workdir, "unitig-popped-unitig-normal-connected-tip.paths.tsv"), t2t)


