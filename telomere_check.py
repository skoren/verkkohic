#!/usr/bin/env python3


import sys
import re

def evaluate_telomers(paths, t2t_contigs, gfa):
    telomere_contigs = {}
    for line in open(t2t_contigs, 'r'):
        forward = (line.split()[0]+"+")
        backward = (line.split()[0]+"-")
        if not (forward in telomere_contigs):
            telomere_contigs[forward] = 0
            telomere_contigs[backward] = 0
        telomere_contigs[forward] += 1
        telomere_contigs[backward] += 1
    lens = {}
    for line in open (gfa, 'r'):
        arr = line.split()
        if arr[0] == 'S':
            lens[arr[1]+'+'] = int(arr[3].split(':')[2])
            lens[arr[1]+'-'] = int(arr[3].split(':')[2])
    t2t_paths = 0

    for line in open(paths, 'r'):
        path = line.split()[1].split(',')
        start = path[0]
        if path[0] in telomere_contigs:
            if len(path) > 1 and path[-1] in telomere_contigs:
                t2t_paths += 1
            elif len(path) == 1 and telomere_contigs[start] >= 2 and lens[start] > 20000000:
                print(f'1-edge component {path[0]}')
                t2t_paths += 1
        for edge in path:
            if edge in telomere_contigs and edge != path[0] and edge != path[-1]:
                print (f'Edge {edge} is reported as telomere but is in the middle of path {line.split()[0]}')
    print(f'Total {t2t_paths} t2t paths')

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

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f'Usage: {sys.argv[0]} <path_file> <telomere_file> <gfa>')
        exit()
    evaluate_telomers(sys.argv[1], sys.argv[2], sys.argv[3])
    check_unclassified(sys.argv[1], sys.argv[2], sys.argv[3])
