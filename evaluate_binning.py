#!/usr/bin/python3

import sys
import os
import random
from os import listdir
from os.path import isfile, join

#hi-c gfa(noseq) trio_colors
comp_file = sys.argv[1]
gfa_file = sys.argv[2]
trio_file = sys.argv[3]

lengths = {}
for line in open (gfa_file,'r'):
    arr = line.split()
    if arr[0] == "S":
        lengths[arr[1]] = int(arr[3].split(':')[-1])
colors = {}
for line in open (trio_file, 'r'):
    arr = line.split()
    name = arr[0]
    color = "a"
    if arr[4] == "#8888FF":
        color = "m"
    elif arr[4] == "#FF8888":
        color = "p"
    colors[name] = color
for line in open (comp_file, 'r'):
    if line[0] != "(":
        continue
    print (line.strip())
    for lines in line.split('}, {'):
        arr = lines.split('\'')
        m_len = 0
        p_len = 0
        u_len = 0
        for contig in arr:
            if contig in lengths:
                l = lengths[contig]
                if not (contig in colors):
                    u_len += l
                else:
                    if colors[contig] == 'm':
                        m_len +=l
                    elif colors[contig] == 'p':
                        p_len +=l
        total_l = m_len + p_len + u_len
        if total_l > 1000000:
            print (f'{total_l}    {m_len/total_l:.3f}/{p_len/total_l:.3f}/{u_len/total_l:.3f}\n')
#    exit()
