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
#pat_from_utig4-2246     utig4-835+,utig4-2245+,utig4-2246+,utig4-2520+,utig4-2521+      PATERNAL

def get_phased_edges(phasedfile):
    phased = set()
    for line in open(phasedfile, 'r'):
        arr = line.split()
        name = arr[0]
        color = "a"
        if arr[4] == "#8888FF" or arr[4] == "#FF8888":
            phased.add(name)
    return phased

#Evaluate paths using only phased edges. If empty set passed - all edges are used.
def evaluate_rukki(rukkifile, triofile, phased_edges, out_f):
    colors = {}
    for line in open(triofile, 'r'):
        arr = line.split()
        name = arr[0]
        color = "a"
        if arr[4] == "#8888FF":
            color = "m"
        elif arr[4] == "#FF8888":
            color = "p"
        colors[name] = color

    if len(phased_edges) == 0:
        for c in colors.keys():
            phased_edges.add(c)
    unassigned = 0
    assigned = 0
    errors = 0
    for line in open(rukkifile):
        strpath = line.split()[1]
        path = strpath.split(',')
        state = "0"
        prev_contig = ""
        for sp in path:
            p = sp[:-1]
            if p in colors:
                if colors[p] == "a" or not p in phased_edges:
                    unassigned += 1
                    continue
                else:
                    if colors[p] != state and state != "0":
                        out_f.write(f"Discordant colors between {prev_contig} {p} !!!\n")
                        out_f.write(strpath +"\n")
                        errors += 1
                    assigned += 1            
                    prev_contig = p
                    state = colors[p]
    out_f.write(f"Among contigs in paths, using uncolored/colored {unassigned}/{assigned} edges, we see {errors} errors\n")
    out_f.close()

if __name__ == "__main__":                
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} <rukkifile> <triofile>')
        exit()
    evaluate_rukki(sys.argv[1], sys.argv[2], set(), sys.stdout)
