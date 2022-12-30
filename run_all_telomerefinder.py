#!/usr/bin/env python3


import sys
import re
import telomere_check
import os

submit_telomere = "/data/antipovd2/devel/vgp-assembly/pipeline/telomere/_submit_telomere.sh"
data_file = "/data/antipovd2//devel/verkkohic/human_only.txt"
data_dir = "/data/antipovd2/data/"

for dir in os.listdir(data_dir):
    workdir = os.path.join(data_dir, dir)
    if os.path.isdir(workdir):
        graph = os.path.join(workdir, "assembly_graph.gfa")
        if os.path.exists(graph):
            run_line = f"{submit_telomere} {graph}"
            print(run_line)
            #os.system(run_line)
            shutil.move("assembly_graph", workdir)