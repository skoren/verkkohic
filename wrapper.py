#!/usr/bin/python3

import sys
import os
import random
from os import listdir
from os.path import isfile, join
out_prefix = "/data/antipovd2/data/"
out_prefix = "/Users/antipovd2/work/verkkohic/res/"
if len(sys.argv) < 3:
    print(f'Usage: {sys.argv[0]} <input_file> <eval_name>')
    exit()
for line in open(sys.argv[1], 'r'):
    arr = line.split()
    print (arr)
    if len(arr)==2 and os.path.exists(arr[0]):
        os.system(f'python3 pipeline.py {arr[0]} {os.path.join(out_prefix, arr[1])} {sys.argv[2]}')
        print (arr[1])
total_count = 0
hiconly_count = 0
datasets = 0
for line in open(sys.argv[1], 'r'):
    arr = line.split()
    if len(arr)==2 and os.path.exists(arr[0]):
        found = False
        datasets += 1
        for nline in open(os.path.join(out_prefix, arr[1], sys.argv[2])):
            narr = nline.strip().split()
            if len(narr) > 1 and narr[-1] == "errors":
                if found:
                    hiconly_count += int(narr[-2])
                else:
                    total_count += int(narr[-2])
                found = True
print (f'Total color switches: {total_count} using all triobinned, {hiconly_count} using only hic colored, totally in {datasets} datasets')
