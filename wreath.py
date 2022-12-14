import itertools
import networkx as nx

# def wreath_product(G,H):
#     #Undirected definition from https://www.sciencedirect.com/science/article/pii/S0096300318304181#bib0029
#     Gn = G.number_of_nodes()
#     Hn = H.number_of_nodes()
#     V = list(itertools.product(itertools.product(list(range(Hn)), repeat=Gn),list(range(Gn))))
#     E = []
#     for e in itertools.combinations(V, 2):
#         if (e[0][1] == e[1][1] and
#            e[0][0][:e[0][1]]+e[0][0][e[0][1]+1:]==e[1][0][:e[1][1]]+e[1][0][e[1][1]+1:] and
#            H.has_edge(list(H)[e[0][0][e[0][1]]],list(H)[e[1][0][e[1][1]]])):
#             E.append(e)
#         elif e[0][0] == e[1][0] and G.has_edge(list(G)[e[0][1]],list(G)[e[1][1]]):
#             E.append(e)
#     return (V,E)

def wreath_product(G,H):
    #directed definition from https://www.sciencedirect.com/science/article/pii/S0195669880800072
    Gn = G.number_of_nodes()
    Hn = H.number_of_nodes()
    V = list(itertools.product(list(range(Gn)),list(range(Hn))))
    E = []
    for e in itertools.combinations(V, 2):
        if (G.has_edge(list(G)[e[0][0]],list(G)[e[1][0]])
           or (e[0][0] == e[1][0] and H.has_edge(list(H)[e[0][1]],list(H)[e[1][1]]))):
            E.append(e)
        if (G.has_edge(list(G)[e[1][0]],list(G)[e[0][0]])
           or (e[0][0] == e[1][0] and H.has_edge(list(H)[e[1][1]],list(H)[e[0][1]]))):
            E.append((e[1],e[0]))
    return (V,E)

def make_graph(X,directed=False):
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    G.add_nodes_from(X[0])
    G.add_edges_from(X[1])
    return G

def edge_file(X,address):
    D={X[0][i]:i for i in range(len(X[0]))}
    f = open(address,'w')
    for e in X[1]:
        f.write(str(D[e[0]])+' '+str(D[e[1]])+'\n')
    f.close()

def run_A2XK3():
    path_poset_address = '/home/phys3smithj/Research/Software/path_poset/'
    G=nx.DiGraph()
    G.add_nodes_from([0,1,2])
    G.add_edges_from([(1,0),(1,2)])
    H=nx.DiGraph()
    H.add_nodes_from([0,1,2])
    H.add_edges_from([(0,1),(1,2),(2,0)])
    X=wreath_product(G,H)
    edge_file(X,'A2XK3.edges')
    print(path_poset_address+'path_poset_homology 9 8 A2XK3.edges A2XK3.homology')
    print('cat A2XK3.homology')
