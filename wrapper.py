#!/usr/bin/python3

import sys
import os
import random
from os import listdir
from os.path import isfile, join
out_prefix = "/data/antipovd2/data/"
for line in open(sys.argv[1], 'r'):
    arr = line.split()
    if len(arr)==2 and os.path.exists(arr[0]):
        os.system(f'python3 pipeline.py {arr[0]} {os.path.join(out_prefix, arr[1])} {sys.argv[2]}')
