#!/usr/bin/python3

import sys
import os
import random
from os import listdir
from os.path import isfile, join

def evaluate_set(contig_set, lengths, colors):
    m_len = 0
    p_len = 0
    u_len = 0
    for contig in contig_set:
        if contig in lengths:
            l = lengths[contig]
            if not (contig in colors):
                u_len += l
            else:
                if colors[contig] == 'm':
                    m_len += l
                elif colors[contig] == 'p':
                    p_len += l
    total_l = m_len + p_len + u_len
    if total_l > 1000000:
        print(f'{total_l}    {m_len / total_l:.3f}/{p_len / total_l:.3f}/{u_len / total_l:.3f}')
        if m_len > 0 and p_len > 0:
            print('BAD')

if len(sys.argv) < 4:
    print(f'Usage: {sys.argv[0]} <hi-c phasing results> <graph.gfa> <trio coloring> [mashmap based chromosome assignment]')
    exit()
#hi-c gfa(noseq) trio_colors
comp_file = sys.argv[1]
gfa_file = sys.argv[2]
trio_file = sys.argv[3]
chrs = {}
if len(sys.argv) > 4:
    chromosomal_file = sys.argv[4]
    for line in open (chromosomal_file):
        arr = line.split()
        contig = arr[0]
        chr = arr[1]
        if not (chr in chrs):
            chrs[chr] = set()
        chrs[chr].add(contig)
lengths = {}
for line in open(gfa_file, 'r'):
    arr = line.split()
    if arr[0] == "S":
        lengths[arr[1]] = int(arr[3].split(':')[-1])
colors = {}
for line in open(trio_file, 'r'):
    arr = line.split()
    name = arr[0]
    color = "a"
    if arr[4] == "#8888FF":
        color = "m"
    elif arr[4] == "#FF8888":
        color = "p"
    colors[name] = color
if len(chrs) == 0:
    chrs['ALL'] = set()
    for name in colors.keys():
        chrs['ALL'].add(name)
for line in open (comp_file, 'r'):
    if line[0] != "(":
        continue
    print()
    print(line.strip())
    for lines in line.split('}, {'):
        arr = lines.split('\'')
        cont_set = set(arr)
        nonempty = 0
        for chr in chrs:
            chr_comp = chrs[chr].intersection(cont_set)
            if len (chr_comp) > 0:
                nonempty+=1
        for chr in chrs:
            chr_comp = chrs[chr].intersection(cont_set)
            if len (chr_comp) > 0:
                print (f'Non empty component on chr {chr}')
                if nonempty == 1:
                    evaluate_set(cont_set, lengths, colors)
                else:
                    print (chr_comp)
                    evaluate_set(chr_comp, lengths, colors)
#    exit()
