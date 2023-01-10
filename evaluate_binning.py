#!/usr/bin/python3

import sys
import os
import random
from os import listdir
from os.path import isfile, join
import networkx as nx
import graph_functions
def evaluate_set(contig_set, lengths, colors):
    m_len = 0
    p_len = 0
    u_len = 0
    m_count = 0
    p_count = 0
    for contig in contig_set:
        if contig in lengths:
            l = lengths[contig]
            if not (contig in colors):
                u_len += l
            else:
                if colors[contig] == 'm':
                    m_len += l
                    m_count += 1
                elif colors[contig] == 'p':
                    p_len += l
                    p_count += 1
    total_l = m_len + p_len + u_len
    if total_l > 0:
        print(f'Total & haplo lengths {total_l}    {m_len / total_l:.3f}/{p_len / total_l:.3f}/{u_len / total_l:.3f}')
        if m_len > 0 and p_len > 0:
            print('BAD')
            print(contig_set)
            if m_len < p_len:
                return m_len, m_count
            else:
                return p_len, p_count
    return 0,0
def evaluate_dataset(hic_file, gfa_file, trio_file, chromosomal_file):
    # hi-c gfa(noseq) trio_colors
    chrs = {}
    if len(chromosomal_file) > 0:
        for line in open(chromosomal_file):
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
    trio_colors = {}
    for line in open(trio_file, 'r'):
        arr = line.split()
        name = arr[0]
        color = "a"
        if arr[4] == "#8888FF":
            color = "m"
        elif arr[4] == "#FF8888":
            color = "p"
        trio_colors[name] = color
    hic_colors = {}
    for line in open(hic_file, 'r'):
        arr = line.split()
        name = arr[0]
        color = "a"
        if arr[4] == "#8888FF":
            color = "m"
        elif arr[4] == "#FF8888":
            color = "p"
        hic_colors[name] = color
    if len(chrs) == 0:
        chrs['ALL'] = set()
        for name in trio_colors.keys():
            chrs['ALL'].add(name)
    G = nx.Graph()
    graph_functions.load_indirect_graph(gfa_file, G)
    graph_functions.remove_large_tangles(G, 200000, 100)
    wrong_len = 0
    wrong_int = 0
    for c in sorted(nx.connected_components(G), key=len, reverse=True):
        #    print("Connected component with %d nodes is: %s" % (len(c), c))
        mset = set()
        pset = set()
        count = 0
        for e in c:
            if e in hic_colors:
                if hic_colors[e] == 'p':
                    pset.add(e)
                else:
                    mset.add(e)
                count += 1
        if len(mset) + len (pset) > 1:
#            print(f"Evaluating component of {len(mset)} and {len(pset)} edges of 1/2 haplotypes")
            w_len, w_int = evaluate_set(mset, lengths, trio_colors)
            wrong_len += w_len
            wrong_int += w_int
            w_len, w_int = evaluate_set(pset, lengths, trio_colors)
            wrong_len += w_len
            wrong_int += w_int
    print(f'length and count of erroneous edges {wrong_len} {wrong_int}')
    return wrong_len, wrong_int
'''
    for line in open(comp_file, 'r'):
        if line.split()[0] != "RES":
            continue
        print()
        print(line.strip())
        for lines in line.split('}, {'):
            arr = lines.split('\'')
            cont_set = set(arr)
            nonempty = 0
            for chr in chrs:
                chr_comp = chrs[chr].intersection(cont_set)
                if len(chr_comp) > 0:
                    nonempty += 1
            for chr in chrs:
                chr_comp = chrs[chr].intersection(cont_set)
                if len(chr_comp) > 0:
                    print(f'Non empty component on chr {chr}')
                    if nonempty == 1:
                        evaluate_set(cont_set, lengths, colors)
                    else:
                        print(chr_comp)
                        evaluate_set(chr_comp, lengths, colors)
    #    exit()
'''



if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            f'Usage: {sys.argv[0]} <hi-c phasing results> <graph.gfa> <trio coloring> [mashmap based chromosome assignment]')
        exit()
    chromosomal_file = ""
    if len (sys.argv) == 5:
        chromosomal_file = sys.argv[4]
    evaluate_dataset(sys.argv[1], sys.argv[2], sys.argv[3], chromosomal_file)