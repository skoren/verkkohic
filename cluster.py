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


#collapse a node, add links from
def collapseOrientedNode (edges, node):
    for pref in ['>', '<']:
        ornode = pref + node
        if not ornode in edges:
            continue
        if not revnode(ornode) in edges:
            continue
        outgoing = set(edges[ornode])
        incoming = set()
        for n in edges[revnode(ornode)]:
            incoming.add(revnode(n))
        for inc in incoming:
            for outg in outgoing:
                edges[inc].add(outg)
    for pref in ['>', '<']:
        ornode = pref + node
        if not revnode(ornode) in edges:
            continue
        incoming = set()
        for n in edges[revnode(ornode)]:
            incoming.add(revnode(n))
        for inc in incoming:
            edges[inc].discard(ornode)
    for pref in ['>', '<']:
        ornode = pref + node
        if ornode in edges:
            edges.pop(ornode)

def fixUnbalanced(part, C, G):
    auxs = [0, 0]
    lens = [0, 0]
    for ind in range (0, 2):
        for node in part[ind]:
            if node.startswith('Aux'):
                auxs[ind] += 1
            elif node in G.nodes:
                lens[ind] += G.nodes[node]['length']
            else:
                print(f"Someting went wrong with node {node}")
#    print (f"In fix unbalanced {lens[0]} {lens[1]}")
    #just _a_ limitation on the swap length is required
    maxswap = min(lens[0]/10, lens[1]/10)
    curswap = 0
    edge_swapped =0
    for ind in range(0, 2):
        #do not move edges _to_ the part where aux edge is present -otherwise it could be swapped with aux
        if auxs[1 - ind] == 0:
#            print(f'processing components {part[0]} {part[1]}')
            changed = True
            while changed:
                changed = False
                for node in part[ind]:
                    if node in G:
                        cost = 0
                        for i in part[ind]:
                            if [i, node] in C.edges:
                                cost += C.edges[i, node]['weight']
                        for j in part[1 - ind]:
                            if [j, node] in C.edges:
                                cost -= C.edges[j, node]['weight']
                        if cost < 0 and curswap + G.nodes[node]['length'] < maxswap:
                            changed = True
                            curswap += G.nodes[node]['length']
                            part[1 - ind].append(node)
                            part[ind].remove(node)
                            edge_swapped +=1
                            print (f"SWAPPED {node}")
                            break
    if curswap > 0:
        print (f"Fixing uneven component, moved {curswap} bases and {edge_swapped} nodes")

def run_clustering (graph_gfa, homologous_nodes, hic_byread, output_dir):
    #TODO: move constants to some more relevant place
    MIN_LEN = 200000  # best result so far with 200000

    MAX_GRAPH_DIST = 10000000  # hi-c links over 10M are believed to be useless

    KLIN_STARTS = 1000  # number of different starts of kernighan lin
    KLIN_ITER = 10000  # number of iterations inside kernighan lin

    MAX_SHORT_COMPONENT = 100 # we remove from consideration each connected compoents of edges < MIN_LEN that is larger than
    MIN_WEIGHT = 10 # Ignore edges with few links
    SIGNIFICANT_MAJORITY = 2.5

    RESULT_FILENAME = "hicverkko.colors.tsv"
    FILTERED_GFA_FILENAME = "filtered.gfa"
    HIC_COMPRESSED_FILENAME = "hic.byread.compressed"
    LOGGING_FILENAME = "hicverkko.log"

    MAX_COV = 100  # temporary coverage cutoff, currently replaced by median coverage from gfa
    FIXED_WEIGHT = 100000  # best result so far with 100000 #currently replaced with max pairwise weight among datasets


    # load the assembly gfa
    G = nx.Graph()
    logging_f = open (os.path.join(output_dir, LOGGING_FILENAME), 'w')
    graph_functions.load_indirect_graph(graph_gfa, G)

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

    filtered_graph = open(os.path.join(output_dir, FILTERED_GFA_FILENAME), 'w')
    tsv_output = os.path.join(output_dir, RESULT_FILENAME)
    tsv_file = open(tsv_output, 'w')
    tsv_file.write("node\tmat\tpat\tmat:pat\tcolor\n")

    translate = open(graph_gfa, 'r')

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
    for l in open(graph_gfa, 'r'):
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

#let's find that nodes that have multiple extension
    multiple_ext = set()
    for node in edges:
        if len(node) > 1:
            multiple_ext.add(node)
    #dirty calculation median coverage, considering that most of the phasing is done
    all_lengths = {}
    total_length = 0
    for node in G.nodes():
        total_length += G.nodes[node]['length']
    sum_l = 0
    med_cov = 0

    for node in sorted(G.nodes(data=True), key=lambda node: node[1]['length']):
        sum_l += node[1]['length']
        if 2*sum_l > total_length:
            med_cov = node[1]['coverage']
            sys.stderr.write (f'Median coverage is {med_cov}\n')
            break
    MAX_COV = med_cov * 1.5
    # no reason for this to be a grpah but meh, why not
    # load pairs of matching nodes based on self-similarity
    matchGraph = nx.Graph()
    translate = open(homologous_nodes, 'r')

    component_colors = {}
    current_color = 0
    for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
        for e in current_component:
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
                    logging_f.write(f"Attempt to restore removed link in former rDNA component between {line[0]} {line[1]} forbidden\n")
                else:
                    logging_f.write(f"Currently not adding graph link between homologous {line[0]} {line[1]}, components {component_colors[line[0]]} "
                          f" and {component_colors[line[1]]}\n")
    #                G.add_edge(line[0], line[1])

    sys.stderr.write("Loaded match info with %d nodes and %d edges\n" % (matchGraph.number_of_nodes(), matchGraph.number_of_edges()))
    translate.close()

    # load hic connections based on mappings, weight corresponds to number of times we see a connection
    hicGraph = nx.Graph()
    hic_file = open(hic_byread, 'r')
    for line in hic_file:
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
    compressed_file = open(os.path.join(output_dir, HIC_COMPRESSED_FILENAME), 'w')
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
    compressed_file.close()

    dists = dict(nx.all_pairs_dijkstra_path_length(G, weight=lambda u, v, d: G.nodes[v]['length']))
    sys.stderr.write("Distances counted\n")


    # connected components decomposition and log the IDs and partition each one
    # currently not going to do the right thing on rDNA component
    long_colors = [set(), set()]
    for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
        logging_f.write("Connected component with %d nodes is: %s\n" % (len(current_component), current_component))

        C = nx.Graph()
        # rebuild the graph from this component using only edges in the hic graph
        C.add_nodes_from(current_component)
    #    if "utig4-1014" not in C.nodes():
    #        continue

        # first we ignore any nodes that are too short or have too few links.
        short = []
        tips = set()
        for n in C.nodes():
            if n not in G:
                sys.stderr.write("Error got a node not in original graph %s !" % (n))
                sys.exit()
    #        if not (n in matchGraph):

            if G.nodes[n]['coverage'] > MAX_COV:
                logging_f.write("While partitoning dropping node %s coverage too high\n" % (n))
                short.append(n)
            elif G.nodes[n]['length'] < MIN_LEN:
    #            sys.stderr.write("While partitoning dropping node %s its too short\n" % (n))
#let's not delete but collapse all short nodes. to "extend" homology information.
                '''neighbours = set()
                for edge in G.edges(n):
                    if edge[0] != n:
                        neighbours.add(edge[0])
                    if edge[1] != n:
                        neighbours.add(edge[1])
#Do we need a separate subgraph for these manipulations?
                for n1 in neighbours:
                    for n2 in neighbours:
                        if n1 != n2:
                            G.add_edge(n1, n2)
#                C.remove_nodes_from(n)
                G.remove_nodes_from(n) '''
                collapseOrientedNode(edges, n)
                short.append(n)

            else:
                good = False
                for e in hicGraph.edges(n):

                    if (e[0] != n or e[1] != n) and hicGraph[e[0]][e[1]]["weight"] > MIN_WEIGHT \
                            and (e[0] in C and e[1] in C):
                        good = True
                if not good:
                    logging_f.write("While partitoning dropping node %s low links count\n" % (n))
                    short.append(n)

        C.remove_nodes_from(short)
        logging_f.write(f'Currently {C.number_of_nodes()} nodes\n')
        # Adding some auxilary vertices to allow slightly unbalanced partitions
        # sqrt/2 is quite arbitrary and subject to change
        aux_nodes = int(math.sqrt(C.number_of_nodes() // 2))
        if C.number_of_nodes() <= 3:
            aux_nodes = 0
        for i in range(0, aux_nodes):
            C.add_node("Aux" + str(i))
        for e in hicGraph.edges(current_component):
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
        for e in C.nodes:
            for pref in ['>', '<']:

                ornode = pref + e
                if not (ornode in edges.keys()):
                    continue
                for neighbour in edges[ornode]:
                    suff = neighbour[1:]
                    if suff in C and matchGraph.get_edge_data(e, suff) == None:
#dirty hack to kill erroneous z-connections issues
                        connection_bonus = FIXED_WEIGHT // 2
                        if len (matchGraph.edges(e)) > 0 and len (matchGraph.edges(suff)) > 0:
                            connection_bonus = 0
#connections with more than one extension are ambiguous and do not deserve bonus
                        if ornode in multiple_ext or revnode(suff) in multiple_ext:
                            connection_bonus = 0

                        w = hicGraph.get_edge_data(e, suff, 0)
                        add_w = 0
                        if w != 0:
                            add_w = w['weight']
                        C.add_edge(e, suff, weight= connection_bonus + add_w)
        #            C.add_edge(e[0], e[1], weight=0 + add_w)
        #            aux_nodes = aux_nodes
        '''for e in G.edges(current_component):
            if e[0] in C and e[1] in C and matchGraph.get_edge_data(e[0], e[1]) == None:
                w = hicGraph.get_edge_data(e[0], e[1], 0)
                add_w = 0
                if w != 0:
                    add_w = w['weight']
                C.add_edge(e[0], e[1], weight=FIXED_WEIGHT + add_w)
    #            C.add_edge(e[0], e[1], weight=0 + add_w)
    #            aux_nodes = aux_nodes
'''
        logging_f.write(f'Currently {C.number_of_nodes()} in current component\n')

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
                logging_f.write(f'HIC edge: {edge} {C.edges[edge]}\n')
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
                    sum_w = 0
                    for i in part[0]:
                        for j in part[1]:
                            if [i, j] in C.edges():
                                # print (f'{i} {j} edge   {C.edges[i,j]}')
                                sum_w += C.edges[i, j]['weight']

                    if (sum_w < best_score):
    #lets forbid Aux only components
                        if check_non_empty(part[0], G) and check_non_empty(part[1], G):
                            logging_f.write(f'Seed {seed} score {sum_w} improved over {best_score}\n')
                            best_part = part
                            best_score = sum_w
                # for ()
                # cut_value, part = nx.stoer_wagner(C)
            #try to move relatively short edges to fix case of unbalanced haplo sizes (complex repeats on only one haplo)
            fixUnbalanced(best_part, C, G)
            logging_f.write(f'RES\t{best_part}\n')

            add_part = [set(), set()]
            only_weights = {}
            for n in current_component:
                if G.nodes[n]['length'] < MIN_LEN and G.nodes[n]['coverage'] < MAX_COV:
                    #            sys.stderr.write("While partitoning dropping node %s its too short\n" % (n))

                    weights = [0, 0]
                    all_weights = [0, 0]
                    total_w = 0
                    for e in hicGraph.edges(n):
                        if (e[0] != n or e[1] != n) and (e[0] in current_component and e[1] in current_component):
                            logging_f.write(f'DEBUG: edge {e} to short {hicGraph[e[0]][e[1]]["weight"]}\n')
                            for ind in range(0, 2):
                                if e[0] in best_part[ind] or e[1] in best_part[ind]:
                                    if hicGraph[e[0]][e[1]]["weight"] > MIN_WEIGHT:
                                        logging_f.write(f'DEBUG: {e} edge to partition {ind}!\n')
                                        weights[ind] += hicGraph[e[0]][e[1]]["weight"]
                                    all_weights[ind] += hicGraph[e[0]][e[1]]["weight"]
                            total_w += hicGraph[e[0]][e[1]]["weight"]
                    if weights[0] > SIGNIFICANT_MAJORITY * weights[1]:
                        add_part[0].add(n)
                        logging_f.write(f"Edge {n} assigned color 0\n")
                    elif weights[1] > SIGNIFICANT_MAJORITY * weights[0]:
                        add_part[1].add(n)
                        logging_f.write(f"Edge {n} assigned color 1\n")
                    else:
                        logging_f.write(f"Edge {n} weights {weights[0]} and {weights[1]} coverage {G.nodes[n]['coverage']}\n")
                        logging_f.write(f"Edge {n} real weights {all_weights[0]} and {all_weights[1]} coverage {G.nodes[n]['coverage']} total weight {total_w}\n")
                        if (all_weights[0] > SIGNIFICANT_MAJORITY * all_weights[1] or all_weights[1] > SIGNIFICANT_MAJORITY * all_weights[0]) and (all_weights[0] > MIN_WEIGHT or all_weights[1]> MIN_WEIGHT):
                            logging_f.write(f"Edge {n} real weights allows to make a decision forbidden by weights!\n")
                            only_weights[n] = all_weights
                elif G.nodes[n]['length'] < MIN_LEN:
                    logging_f.write(f"Edge {n} did not pass coverage check\n")

            logging_f.write(f'Added {len(add_part[0]) + len (add_part[1])} short edges to bipartition\n')
            logging_f.write(f'RES\t{add_part}\n')
            best_part[0].update(add_part[0])
            best_part[1].update(add_part[1])
            right = False
            for ind in range (0, 2):
                for contig in best_part[ind]:
                    if contig in current_component:
                        if ind == 1:
                            tsv_file.write(f'{contig}\t0\t100000\t0:100000\t#8888FF\n')
                        else:
                            tsv_file.write(f'{contig}\t100000\t0\t100000:0\t#FF8888\n')
            for contig in only_weights.keys():
                tsv_file.write(f'{contig}\t{only_weights[contig][0]}\t{only_weights[contig][1]}\t{only_weights[contig][0]}:{only_weights[contig][1]}\t#88FF88\n')
    #all_weights should be added here!

    #        for ind in [0, 1]:
    #            for j in best_part[ind]:
    #               long_colors[ind].add(j)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f'Usage: {sys.argv[0]} graph.gfa homologous_nodes.matches hic_byread output_dir')
        exit()
    run_clustering(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
