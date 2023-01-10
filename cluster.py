import sys
import random
import networkx as nx
import math
import os
import graph_functions
from networkx.algorithms import community

def check_non_empty(part, G):
    for p in part:
        if p in G.nodes:
            return True
    return False


MIN_LEN = 200000  # best result so far with 200000
FIXED_WEIGHT = 100000  # best result so far with 100000 
#currently replaced with max pairwise weight among datasets

MAX_GRAPH_DIST = 10000000  # hi-c links over 10M are believed to be useless
MAX_COV = 100  # temporary coverage cutoff, currently replaced by median coverage from gfa

KLIN_STARTS = 1000  # number of different starts of kernighan lin
KLIN_ITER = 10000  # number of iterations inside kernighan lin

MAX_SHORT_COMPONENT = 50 # we remove from consideration each connected compoents of edges < MIN_LEN that is larger than

MIN_WEIGHT = 10 # Ignore edges with few links
#this cutoff
print(nx.__version__)
print(nx.__file__)


def revnode(n):
    assert len(n) >= 2
    assert n[0] == "<" or n[0] == ">"
    return (">" if n[0] == "<" else "<") + n[1:]

#Currently not in use anymore
def IsTip(node, edges):
    for pref in ['>', '<']:
        ornode = pref + node
        if ornode not in edges or len(edges[ornode]) == 0:
            if revnode(ornode) not in edges:
                continue
#rc to the tip predeccor
            for edge in edges[revnode(ornode)]:
#alternative edge to tip that should be not a deadend
                for alt_tip in edges[revnode(edge)]:
                    if alt_tip in edges and len(edges[alt_tip])> 0:
                        print (f'Tip {node}')
                        return True
    return False


if len(sys.argv) != 5:
    print(f'Usage: {sys.argv[0]} graph.gfa homologous_nodes.matches hic_byread output_dir')
    exit()
# load the assembly gfa
G = nx.Graph()
graph_functions.load_indirect_graph(sys.argv[1], G)

degrees = [val for (node, val) in G.degree()]
mean = sum(degrees) / G.number_of_nodes()
variance = sum([((x - mean) ** 2) for x in degrees]) / G.number_of_nodes()
res = variance ** 0.5
sys.stderr.write("Loaded a graph with %d nodes and %d edges avg degree %f and stdev %f max is %f\n" % (
G.number_of_nodes(), G.number_of_edges(), mean, res, mean + 5 * res))

#Store rDNA component, not to add additional links from matchGraph
largest_component = max(nx.connected_components(G), key=len)
sys.stderr.write(f"Found an rDNA huge component of {len(largest_component)} edges\n")

#Here we remove large connected components of short edge, just to exclude rDNA cluster
delete_set = graph_functions.remove_large_tangles(G, MIN_LEN, MAX_SHORT_COMPONENT)

filtered_graph = open(os.path.join(sys.argv[4], "filtered.gfa"), 'w')

translate = open(sys.argv[1], 'r')

for line in translate:
    arr = line.split()
    if arr[0] == "S":
        if not (arr[1] in delete_set):
            filtered_graph.write(line)
    if arr[0] == "L":
        if not(arr[1] in delete_set) and not(arr[3] in delete_set):
            filtered_graph.write(line)
#loading oriented graph
#TODO: only place which used it is likely outdated
nodelines = []
nodelens = []
edges = {}
for l in open(sys.argv[1], 'r'):
    parts = l.strip().split('\t')
    if parts[0] == 'S':
        nodelines.append((parts[1], l.strip()))
#        nodelens[parts[1]] = len(parts[2])
    elif parts[0] == 'L':
        fromnode = (">" if parts[2] == "+" else "<") + parts[1]
        tonode = (">" if parts[4] == "+" else "<") + parts[3]
#        edgelines.append((fromnode, tonode, l.strip()))
        if fromnode not in edges:
            edges[fromnode] = set()
        if revnode(tonode) not in edges:
            edges[revnode(tonode)] = set()
        edges[fromnode].add(tonode)
        edges[revnode(tonode)].add(revnode(fromnode))

#dirty calculation median coverage, considering that most of the phasing is done
all_lengths = {}
total_length = 0
for node in G.nodes():
    total_length += G.nodes[node]['length']
sum = 0
med_cov = 0

for node in sorted(G.nodes(data=True), key=lambda node: node[1]['length']):
    print (node[1])
    sum += node[1]['length']
    if 2*sum > total_length:
        med_cov = node[1]['coverage']
        print (f'Median coverage is {med_cov}')
        break
MAX_COV = med_cov * 1.5
# no reason for this to be a grpah but meh, why not
# load pairs of matching nodes based on self-similarity
matchGraph = nx.Graph()
translate = open(sys.argv[2], 'r')

component_colors = {}
current_color = 0
for c in sorted(nx.connected_components(G), key=len, reverse=True):
    print (f" component {current_color} len {len(c)}")
    for e in c:
        component_colors[e] = current_color
    current_color += 1
for line in translate:
    if "#" in line:
        continue
    line = line.strip().split()
    if len(line) < 2:
        continue
    if line[0] == line[1]:
        continue
    matchGraph.add_edge(line[0], line[1])
    #Adding link between matched edges to include separated sequence to main component

    if line[0] in G.nodes and line[1] in G.nodes:
        if component_colors[line[0]] != component_colors[line[1]]:
            if line[0] in largest_component and line[1] in largest_component:
                print(f"Attempt to restore removed link in former rDNA component between {line[0]} {line[1]} forbidden")
            else:
                print(f"Currently not adding graph link between homologous {line[0]} {line[1]}, components {component_colors[line[0]]} "
                      f" and {component_colors[line[1]]}")
#                G.add_edge(line[0], line[1])

sys.stderr.write("Loaded match info with %d nodes and %d edges\n" % (matchGraph.number_of_nodes(), matchGraph.number_of_edges()))
translate.close()

# load hic connections based on mappings, weight corresponds to number of times we see a connection
hicGraph = nx.Graph()
translate = open(sys.argv[3], 'r')
for line in translate:
    if "#" in line:
        continue
    line = line.strip().split()

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
compressed_file = open(os.path.join(sys.argv[4], "hic.byread.compressed"), 'w')
max_w = 0
for node1, node2 in hicGraph.edges():
    if max_w < hicGraph[node1][node2]["weight"]:
        max_w = hicGraph[node1][node2]["weight"]
    compressed_file.write(f'X {node1} {node2} {hicGraph[node1][node2]["weight"]}\n')
FIXED_WEIGHT = max_w
sys.stderr.write(f'Constant for neighbouring edges set to be  {FIXED_WEIGHT}, for homologous edges {-10 * FIXED_WEIGHT} \n')


sys.stderr.write("Loaded hic info with %d nodes and %d edges\n" % (hicGraph.number_of_nodes(), hicGraph.number_of_edges()))
# for node1, node2 in hicGraph.edges():
#    sys.stderr.write("Edge from %s to %s of weight %s\n"%(node1, node2, hicGraph.get_edge_data(node1, node2)))
translate.close()

dists = dict(nx.all_pairs_dijkstra_path_length(G, weight=lambda u, v, d: G.nodes[v]['length']))
sys.stderr.write("Distances counted\n")
# connected components decomposition and log the IDs and partition each one
# currently not going to do the right thing on rDNA component

for c in sorted(nx.connected_components(G), key=len, reverse=True):
    print("Connected component with %d nodes is: %s" % (len(c), c))

    C = nx.Graph()
    # rebuild the graph from this component using only edges in the hic graph
    C.add_nodes_from(c)
#    if "utig4-1014" not in C.nodes():
#        continue
    Subgraph = G.subgraph(c).copy()


    # first we ignore any nodes that are too short or have too few links.
    short = []
    tips = set()
    for n in C.nodes():
        if n not in G:
            sys.stderr.write("Error got a node not in original graph %s !" % (n))
            sys.exit()
#        if not (n in matchGraph):

        if G.nodes[n]['length'] < MIN_LEN:
#            sys.stderr.write("While partitoning dropping node %s its too short\n" % (n))
            short.append(n)

        elif G.nodes[n]['coverage'] > MAX_COV:
            sys.stderr.write("While partitoning dropping node %s coverage too high\n" % (n))
            short.append(n)
        else:
            good = False
            for e in hicGraph.edges(n):

                if (e[0] != n or e[1] != n) and hicGraph[e[0]][e[1]]["weight"] > MIN_WEIGHT \
                        and (e[0] in C and e[1] in C):
                    good = True
            if not good:
                sys.stderr.write("While partitoning dropping node %s low links count\n" % (n))
                short.append(n)

    C.remove_nodes_from(short)
    # Adding some auxilary vertices to allow slightly unbalanced partitions
    # sqrt/2 is quite arbitrary and subject to change
    aux_nodes = int(math.sqrt(C.number_of_nodes() // 2))
    if C.number_of_nodes() <= 3:
        aux_nodes = 0
    for i in range(0, aux_nodes):
        C.add_node("Aux" + str(i))
    for e in hicGraph.edges(c):
        # currently only added edges if these nodes are in the component and not matches (homologous) but should allow links to singletons too (to phase disconnected nodes correctly)
        if e[0] in C and e[1] in C and matchGraph.get_edge_data(e[0], e[1]) == None:
            # if edges are too distant in graph, hi-c info is trash
            # using homologous edges when counting distances to help with coverage gaps
            similar_edges = [set(), set()]
            for ind in range(0, 2):
                for match_edge in matchGraph.edges(e[ind]):
                    similar_edges[ind].add(match_edge[0])
                    similar_edges[ind].add(match_edge[1])
                similar_edges[ind].add(e[ind])
#TODO: likely this is not needed anymore since we already added links between homologous edges.
            for e0like in similar_edges[0]:
                for e1like in similar_edges[1]:
                    if e0like in dists and e1like in dists[e0like] and dists[e0like][e1like] < MAX_GRAPH_DIST + G.nodes[e1like]['length']:
                        C.add_edge(e[0], e[1], weight=hicGraph[e[0]][e[1]]['weight'])
                        break

    # neighboring edges are encouraged to belong to the same component
    #Currently not in use
    for e in G.edges(c):
        if e[0] in C and e[1] in C and matchGraph.get_edge_data(e[0], e[1]) == None:
            w = hicGraph.get_edge_data(e[0], e[1], 0)
            add_w = 0
            if w != 0:
                add_w = w['weight']
#            C.add_edge(e[0], e[1], weight=FIXED_WEIGHT + add_w)
            C.add_edge(e[0], e[1], weight=0 + add_w)
#            aux_nodes = aux_nodes


    if C.number_of_nodes() > 1:
#TODO: why not just iterate on matchGraph.edges()?
        for n in C.nodes():
            if n in matchGraph:
                for ec in matchGraph.edges(n):
                    # while we build the initial partition give a big bonus edge for putting the homologous nodes into different partitions
                    if ec[0] in C and ec[1] in C and ec[0] < ec[1]:
                        C.add_edge(ec[0], ec[1], weight=-10 * FIXED_WEIGHT)
                        C.add_edge(ec[1], ec[0], weight=-10 * FIXED_WEIGHT)
        for edge in C.edges:
            print(f'EEEDGE {edge} {C.edges[edge]}')
        best_score = FIXED_WEIGHT * C.number_of_nodes() * C.number_of_nodes()
        
        for seed in range(0, KLIN_STARTS):  # iterate on starting partition
            random.seed(seed)
            p1 = []
            p2 = []
            #debug code
            if False: #"utig4-1014" in C.nodes():
                shastas_set = {'utig4-1448', 'utig4-2694', 'utig4-1015', 'utig4-1111', 'utig4-631', 'utig4-341', 'utig4-1444', 'utig4-1449', 'utig4-872', 'utig4-1558', 'utig4-1607', 'utig4-342', 'utig4-1850', 'utig4-1113', 'utig4-1531', 'utig4-1076', 'utig4-1684', 'utig4-371', 'utig4-874', 'utig4-1062'}
                for n in C.nodes():
                    if n in shastas_set:
                        p1.append(n)
                    else:
                        p2.append(n)
            # make an initail guess at partitioning by putting homologous nodes into opposite clusters
            else:
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
                                if random.random() <= 0.5:
                                    p1.append(n)
                                else:
                                    p2.append(n)
                    else:
                        if random.random() <= 0.5:
                            p1.append(n)
                        else:
                            p2.append(n)
#            print("Initial partitions are %s and %s" % (set(p1), set(p2)))
            if len(p1) * len(p2) > 0:
                part = community.kernighan_lin.kernighan_lin_bisection(C, partition=[set(p1), set(p2)], max_iter=KLIN_ITER,
                                                                       weight='weight', seed=seed)
                sum = 0
                for i in part[0]:
                    for j in part[1]:
                        if [i, j] in C.edges():
                            # print (f'{i} {j} edge   {C.edges[i,j]}')
                            sum += C.edges[i, j]['weight']

                if (sum < best_score):
#lets forbid Aux only components
                    if check_non_empty(part[0], G) and check_non_empty(part[1], G):
                        print(f'Seed {seed} score {sum} improved over {best_score}')
                        best_part = part
                        best_score = sum
            # for ()
            # cut_value, part = nx.stoer_wagner(C)
        print(f'RES\t{best_part}')
        print(best_score)
