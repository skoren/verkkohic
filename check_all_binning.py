#!/usr/bin/env python3


import sys
import re
import telomere_check
import os
import evaluate_binning

data_file = "/data/antipovd2//devel/verkkohic/human_only.txt"
data_dir = "/data/antipovd2/data/"
data_dir = "/Users/antipovd2/work/verkkohic/data/"
res_dir = "/Users/antipovd2/work/verkkohic/res/"


for dir in os.listdir(data_dir):
    workdir = os.path.join(data_dir, dir)
    if os.path.isdir(workdir):
        print (f'dataset: {dir}')
        hic_file = os.path.join(res_dir, dir, sys.argv[1], "unitig-popped-unitig-normal-connected-tip.colors.csv")
        trio_file = os.path.join(workdir, "unitig-popped-unitig-normal-connected-tip.colors.csv")
        gfa = os.path.join(workdir, "unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa")
        if os.path.exists(hic_file):
            evaluate_binning.evaluate_dataset(hic_file, gfa, trio_file, "")


