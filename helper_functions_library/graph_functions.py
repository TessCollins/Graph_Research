from sage.all import *
from .matrix_functions import *

##################################################

def lib_make_group_graph(p,n,conn_set):
    """
    Will make a graph on Z_p^n with specified connection set.
    
    Args:
    - p for Z_p^n
    - n for Z_p^n
    - Connection set of vectors
    
    Returns:
    - A nice graph :-)
    """
    
    if len(set(list(map(len,conn_set)))) != 1: #check if vectors in connection set are valid
        raise ValueError("Vectors in the connection set must be of the same length")
    C = list(map(vector,conn_set)) #make the connection set into vectors
    
    S=Permutations([i for i in range(p)]*n,n).list() #permuatations of possibilities as a list
    veclist = [] #empty list
    
    for i in range(len(S)):
        veclist.append(vector(Integers(p),S[i])) #make a list of the vectorified possibilities
                       
    G = Graph() #empty graph
    
    for vec in veclist:
        for c in C:
            G.add_edge([tuple(vec),tuple(vec+c)]) #add edges between the vertices based on the connection set

    return G
    
##################################################

def lib_make_Hadamard(N):
    """
    Makes Hadamard graph G_N
    
    Args:
    -N, for N = 4*k
    
    Returns:
    -Hadamard graph G_N
    """
    if N%4 != 0: #check to make sure there's the right number of vertices
        raise ValueError("Graph needs to have a number of vertices who are a multiple of 4.")
    VG = Permutations([i for i in range(2)]*N,N).list() #vector names as a list
    G = Graph()
    
    veclist = [] #empty list
    for i in range(len(VG)):
        veclist.append(vector(VG[i])) #make a list of the vectorified possibilities

    for j in range(2**N):
        for k in range(j,2**N):
            if sum((veclist[j]+veclist[k])%2) == N/2: #if vectors have Hamming distance N/2, put an edge between them
                G.add_edge([tuple(veclist[j]),tuple(veclist[k])])
    return G

##################################################

def lib_make_all_distance_graphs(G):
    """
    Will give all distance graphs of graph G, but not same vertex set.
    
    Args:
    - Graph G
    
    Returns:
    - List of all distance graphs of G, for all possible distances.
    """

    graphlist = []
    d = G.diameter()
    dm = G.distance_matrix()
    n = len(G.vertices())
    for j in range(d+1):
        dm1=copy(dm)
        for j1 in range(n):
            for j2 in range(n):
                if dm[j1,j2] == j:
                    dm1[j1,j2] = 1
                else:
                    dm1[j1,j2] = 0
        graphlist.append(Graph(dm1))
    
    return graphlist

##################################################

def lib_make_distance_graph(G,d,is_mats):
    """
    Makes a distance graph of a particular graph, with the same vertex set.
    
    Args:
    - Graph G to make a distance graph of.
    - Desired distance d for distance graph.
    - 0 or 1 value to indicate if the vertices should be strings or tuples (resp.)
        (tuples can be matrices very easy, strings are harder)
        
    Returns:
    - Distance graph G_d whose vertices are named the same as G
    """
    
    def matify(G):
        G_d = Graph()
        for j1 in range(n):
            for j2 in range(j1,n):
                v1=Matrix(map(sage_eval,G.vertices()[j1]))
                v2=Matrix(map(sage_eval,G.vertices()[j2]))
                if sum(list((v1+v2)%2)[i in range(5)]) == d:
                    G_d.add_edge([tuple(vector(v1)),tuple(vector(v2))])
        return G_d
    
    n = len(G.vertices())
    
    if is_mats == 1:
        G_d = matify(G)
    else:
        G_d = Graph()
        for j1 in range(n):
            for j2 in range(j1,n):
                if sum(list((Matrix(map(sage_eval,G.vertices()[j1]))
                             +Matrix(map(sage_eval,G.vertices()[j2])))%2)[i in range(5)]) == d:
                    G_d.add_edge([G.vertices()[j1],G.vertices()[j2]])
    return G_d

##################################################