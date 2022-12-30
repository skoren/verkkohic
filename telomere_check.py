#!/usr/bin/env python3


import sys
import re

def evaluate_telomers(paths, t2t_contigs):
    telomere_contigs = set()
    for line in open(sys.argv[2], 'r'):
        telomere_contigs.add(line.split()[0]+"+")
        telomere_contigs.add(line.split()[0]+"-")

    t2t_paths = 0

    for line in open(sys.argv[1], 'r'):
        path = line.split()[1].split(',')
        if len(path) > 1 and path[0] in telomere_contigs and path[-1] in telomere_contigs:
            t2t_paths += 1
        for edge in path:
            if edge in telomere_contigs and edge != path[0] and edge != path[-1]:
                print (f'Edge {edge} is reported as telomere but is in the middle of path {line.split()[0]}')
    print(f'Total {t2t_paths} t2t paths')
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} <path_file> <telomere_file>')
        exit()
    evaluate_telomers(sys.argv[1], sys.argv[2])
