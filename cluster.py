import sys
import random
import networkx as nx
import math
from networkx.algorithms import community

MIN_LEN=200000 # best result so far with 200000
FIXED_WEIGHT=1000000 # best result so far with 100000
MAX_COV = 100 #temporary coverage cutoff, should be rewritten
KLIN_STARTS = 1000 #number of different starts of kernighan lin
KLIN_ITER = 10000 #number of iterations inside kernighan lin
print(nx.__version__)
print(nx.__file__)
G = nx.Graph()

if len(sys.argv)!=4:
    print (f'Usage: {sys.argv[0]} graph.gfa homologous_nodes.matches hic_byread' )
    exit()
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
degrees = [val for (node, val) in G.degree()]
mean = sum(degrees) / G.number_of_nodes()
variance = sum([((x - mean) ** 2) for x in degrees]) / G.number_of_nodes()
res = variance ** 0.5
sys.stderr.write("Loaded a graph with %d nodes and %d edges avg degree %f and stdev %f max is %f\n"%(G.number_of_nodes(), G.number_of_edges(), mean, res, mean+5*res))
translate.close()

#no reason for this to be a grpah but meh, why not
# load pairs of matching nodes based on self-similarity
matchGraph = nx.Graph()
translate = open(sys.argv[2], 'r')
for line in translate:
   if "#" in line:
      continue
   line=line.strip().split()

   if len(line) < 2:
      continue

   if line[0] == line[1]:
      continue

   matchGraph.add_edge(line[0], line[1])
sys.stderr.write("Loaded match info with %d nodes and %d edges\n"%(matchGraph.number_of_nodes(), matchGraph.number_of_edges()))
translate.close()

#load hic connections based on mappings, weight corresponds to number of times we see a connection
hicGraph = nx.Graph()
translate = open(sys.argv[3], 'r')
for line in translate:
   if "#" in line:
      continue
   line=line.strip().split()

   if len(line) < 3:
      continue

   if line[1] == line[2]:
      continue

   hicGraph.add_node(line[1])
   hicGraph.add_node(line[2])
   if len(line) > 3:
      add_w = int(line[3])
   else:
      add_w = 1
   w = hicGraph.get_edge_data(line[1], line[2], 0)
   if w == 0:
      hicGraph.add_edge(line[1], line[2], weight=add_w)
   else:
      w = w['weight'] + add_w
      hicGraph[line[1]][line[2]]['weight'] = w
compressed_file = open("hic.byread.compressed", 'w')
for node1, node2 in hicGraph.edges():
   compressed_file.write(f'X {node1} {node2} {hicGraph[node1][node2]["weight"]}\n')
sys.stderr.write("Loaded hic info with %d nodes and %d edges\n"%(hicGraph.number_of_nodes(), hicGraph.number_of_edges()))
for node1, node2 in hicGraph.edges():
    sys.stderr.write("Edge from %s to %s of weight %s\n"%(node1, node2, hicGraph.get_edge_data(node1, node2)))
translate.close()

# connected components decomposition and log the IDs and partition each one
# currently not going to do the right thing on rDNA component
for c in sorted(nx.connected_components(G), key=len, reverse=True):
   print("Connected component with %d nodes is: %s"%(len(c), c))

   C = nx.Graph()
   # rebuild the graph from this component using only edges in the hic graph
   C.add_nodes_from(c)

   # first we ignore any nodes that are too short
   short = []
   for n in C.nodes():
      if n not in G:
         sys.stderr.write("Error got a node not in original graph %s !"%(n))
         sys.exit()
      if G.nodes[n]['length'] < MIN_LEN:
         sys.stderr.write("While partitoning dropping node %s its too short\n"%(n))
         short.append(n)
      elif G.nodes[n]['coverage'] > MAX_COV:
         sys.stderr.write("While partitoning dropping node %s coverage too high\n"%(n))
         short.append(n)

   C.remove_nodes_from(short)
    #Adding some auxilary vertices to allow slightly unbalanced partitions
    #sqrt/2 is quite arbitrary and subject to change
   aux_nodes = int(math.sqrt (C.number_of_nodes() // 2))
   for i in range (0, aux_nodes):
        C.add_node ("Aux" + str(i))
   for e in hicGraph.edges(c):
      # currently only added edges if these nodes are in the component and not matches (homologous) but should allow links to singletons too (to phase disconnected nodes correctly)
      if e[0] in C and e[1] in C and matchGraph.get_edge_data(e[0], e[1]) == None:
         C.add_edge(e[0], e[1], weight=hicGraph[e[0]][e[1]]['weight'])
         # print("Adding edge between %s and %s of weight %s"%(e[0], e[1], hicGraph[e[0]][e[1]]['weight']))
   for e in G.edges(c):
      if e[0] in C and e[1] in C and matchGraph.get_edge_data(e[0], e[1]) == None:
         C.add_edge(e[0], e[1], weight=FIXED_WEIGHT)
   if C.number_of_nodes() > 1:
      for n in C.nodes():
         if n in matchGraph:
            for ec in matchGraph.edges(n):
               # while we build the initial partition give a big bonus edge for putting the homologous nodes into different partitions
               if ec[0] in C and ec[1] in C:
                  C.add_edge(ec[0], ec[1], weight=-10 * FIXED_WEIGHT)
      best_score = FIXED_WEIGHT* C.number_of_nodes() * C.number_of_nodes()
      for seed in range(0, KLIN_STARTS):  # iterate on starting partition
         random.seed(seed)
         p1 = []
         p2 = []
         # make an initail guess at partitioning by putting homologous nodes into opposite clusters
         for n in C.nodes():
            if n in matchGraph:
               for ec in matchGraph.edges(n):
                  if ec[0] == n and ec[1] in p1:
                    if n not in p1:
                       p2.append(n)
                  elif ec[0] == n and ec[1] in p2:
                     if n not in p2:
                       p1.append(n)
                  elif ec[1] == n and ec[0] in p1:
                      if n not in p1:
                         p2.append(n)
                  elif ec[1] == n and ec[0] in p2:
                     if n not in p2:
                        p1.append(n)
                  elif n not in p1 and n not in p2:
                     p1.append(n)
            else:
               if random.random() <= 0.5:
                  p1.append(n)
               else:
                  p2.append(n)
         print("Initial partitions are %s and %s"%(set(p1), set(p2)))
         part=community.kernighan_lin.kernighan_lin_bisection(C, partition=[set(p1), set(p2)], max_iter=KLIN_ITER, weight='weight', seed=seed)
         sum = 0
         for i in part[0]:
           for j in part[1]:
               if [i,j] in C.edges():
                   #print (f'{i} {j} edge   {C.edges[i,j]}')
                   sum += C.edges[i,j]['weight']

         if (sum < best_score):
            print(f'Seed {seed} score {sum} improved over {best_score}')
            best_part = part
            best_score = sum
         #for ()
         #cut_value, part = nx.stoer_wagner(C)
      print(best_part)
