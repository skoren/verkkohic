#!/usr/bin/python3
import sys
import random
import networkx as nx
import math
from networkx.algorithms import community

MIN_LEN=200000 # best result so far with 200000
FIXED_WEIGHT=100000 # best result so far with 100000
MAX_GRAPH_DIST = 10000000 #hi-c links over 10M are believed to be useless
MAX_COV = 100 #temporary coverage cutoff, should be rewritten
KLIN_STARTS = 100 #number of different starts of kernighan lin
KLIN_ITER = 10000 #number of iterations inside kernighan lin
print(nx.__version__)
print(nx.__file__)
G = nx.Graph()

if len(sys.argv)!=3:
    print (f'Usage: {sys.argv[0]} graph.gfa shasta_clustering')
    exit ()
#load the assembly gfa
translate = open(sys.argv[1], 'r')
for line in translate:
   if "#" in line:
      continue
   line=line.strip().split()

   if line[0] == "S":
      G.add_node(line[1], length=int(line[3][5:]), coverage=float(line[5][5:]))
   elif  line[0] == "L":
      if line[1] not in G or line[3] not in G:
         sys.stderr.write("Warning, skip link between nodes not in graph:%s"%(line))
         sys.exit(1)
      G.add_edge(line[1], line[3])
colors = {}
for line in open (sys.argv[2]):
    arr = line.split()
    contig = arr[0]
    if arr[1] == str(100000):
        haplo = 'p'
    else:
        haplo = 'm'
    colors[contig] = haplo
count = 0

for c in sorted(nx.connected_components(G), key=len, reverse=True):
#    print("Connected component with %d nodes is: %s" % (len(c), c))
    mset = set()
    pset = set()
    for e in c:
        if e in colors:
            if colors[e] == 'p':
                pset.add(e)
            else:
                mset.add(e)
            count += 1
    if len(pset) > 0 and len(mset) > 0:
        print(f'({pset}, {mset})')
print (f'Total {count} edges classified')
