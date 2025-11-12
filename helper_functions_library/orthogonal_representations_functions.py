from sage.all import *

from .matrix_functions import *
from .graph_functions import *

##################################################

def lib_degen_find(G1):
    """
    Finds minimum degeneracy of a graph and the vertex ordering that admits that degeneracy
    
    Args:
    - A graph
    
    Returns:
    - Minimum degeneracy of the graph using the two different methods
    - The vertex ordering that admits that degeneracy
    """
    def degen_find_min(G1):
        """
        Finds the degeneracy of a given graph by plucking off minimum degree vertices,
        and placing them at the end of the optimal ordering

        Args:
        - The graph

        Returns:
        - Degeneracy of the graph
        - Ordering of vertices that admits that degeneracy
        """
        ### Finding the optimal ordering ###
        G=G1.copy()
        n = len(G.vertices());
        ordering = [None for i in range(n)];
        j = 0; #counts how many vertices we've deleted so far

        while j < n:
            degs = G.degree(i for i in G.vertices()) #degrees of all vertices in graph
            mind = min(degs); #mininum degree
            for k in G.vertices():
                if G.degree(k) == mind:
                    ordering[n-1-j] = k; #if the vertex has min degree, put it at end of ordering
                    G.delete_vertex(k); #and remove that vertex from the graph
                    j = j+1
                    break

        ### Finding the degeneracy of that ordering###
        npnv=[] #empty list of the number of previous neighbors

        for i in range(n):
            vi = ordering[i] #the vertex associated with index i in the optimal ordering
            npn=0 #number of previous neighbors to vertex vi
            for j in range(i): #for all previous vertices in the optimal ordering
                vj = ordering[j]
                if vj in G1.neighbors(vi): #if the vertex is adjacent to vertex vi, increase the number of previous neighbors of i
                    npn+=1
            npnv.append(npn) #put the number of vi's previous neighbors into the list of number of previous neighbors
        degen = max(npnv) #the degeneracy is the max of all of those previous neighbors

        return ordering,degen

    def degen_find_max(G1):
        """
        Finds the degeneracy of a given graph by ordering the vertices from highest to lowest degree for the optimal ordering

        Args:
        - The graph

        Returns:
        - Degeneracy of the graph
        - Ordering of vertices that admits that degeneracy
        """

        ### Finding the optimal ordering ###
        G=G1.copy()
        n = len(G.vertices());
        ordering=[];
        j = 0; #counts how many vertices we've deleted so far
        degs = G.degree(i for i in G.vertices()) #degrees of all vertices in graph

        while j < n:
            maxd=max(degs)
            for k in G.vertices():
                if G1.degree(k) == maxd:
                    ordering.append(k)
                    degs.remove(maxd)
                    G.delete_vertex(k)
                    j=j+1
                    break

        ### Finding the degeneracy of that ordering###
        npnv=[] #empty list of the number of previous neighbors
        for i in range(n):
            vi = ordering[i] #the vertex associated with index i in the optimal ordering
            npn=0 #number of previous neighbors to vertex vi
            for j in range(i): #for all previous vertices in the optimal ordering
                vj = ordering[j]
                if vj in G1.neighbors(vi): #if the vertex is adjacent to vertex vi, increase the number of previous neighbors of i
                    npn+=1
            npnv.append(npn) #put the number of vi's previous neighbors into the list of number of previous neighbors
        degen = max(npnv) #the degeneracy is the max of all of those previous neighbors

        return ordering,degen
    
    degss = [degen_find_max(G1), degen_find_min(G1)]

    mindegen = min(degss[i][1] for i in range(1))
    minorder = [degss[i][0] for i in range(1) if degss[i][1] == mindegen]
    
    return mindegen, minorder

##################################################

def lib_modified_LSS(G,Vorder,d):
    """
    Modified LSS, finds orthogonal representation of graph using its degeneracy
    Finds an orthogonal representation of the graph
    
    Args:
    - The graph
    - An ordering of the vertices
    - Desired dimension of the orthogonal representation
    
    Returns:
    - An orthogonal representation of the graph
    """
    
    n = len(Vorder) #number of vertices
    us = [None for i in range(n)] #empty list in which we can put uniformly random vectors
    fs=[] #empty list into which we shall put orthogonal vectors

    for i in range(n): # assign uniformly random vectors to everything
        us[i] = vector(lib_make_uniform_rand_vec(d)) 
        
    fs.append(Matrix(us[0])); #first one is just the regular random guy
 
    for i in range(1,n):
        prev_f = []; #gather the previously assigned f vectors (the orthogonal ones)
        vi = Vorder[i]; #the vertex in the ith place
        
        nnp=0; #number of neighbors previous in the list
        for j in range(i): #for all the vertices previous in the listed order
            vj = Vorder[j]; #the vertex in the jth place
            
            if vj in G.neighbors(vi): #if the previous vertex is a neighbor of this one
                nnp+=1; #add one to the count of previous neighbors
                prev_f.append(vector(fs[j])) #collect its vector in the list of previous f vectors
                bassisor,b = Matrix(prev_f).gram_schmidt() #gram schmidtify the previous f vectors into an orthogonal basis
        
        newv = Matrix(us[i]) #the new vector will be the randomly generated u and then we'll mess with it
        if nnp != 0: #if there were neighbors previous in the list
    
            for v1 in bassisor: #for all vectors in that basis of prev fs
                v1vec=Matrix(v1); #make a matrix out of it
                newv -= Matrix((newv*v1vec.T)*v1vec) #subtract the projection of u onto that vec from u
            
        fs.append(newv) #normalize and add it to the list of f vectors
    return fs

##################################################

def lib_reconstructing(vecs,ordering):
    """
    Reconstructs a graph from vector adjacencies
    
    Args:
    - List of vectors
    - Ordering of vertices
    
    Returns:
    - A graph with as many vertices as there are vectors, where vertices are adjacent if vectors are orthogonal
    """
    
    n = len(vecs) #how many vectors we have
    G_r = Graph() #empty graph
    G_r.add_vertices([i for i in range(n)]) #add vertices to the empty graph
    for i in range(n): #for all vertices
        for j in range(i+1,n):
            b=Matrix(vecs[i])*Matrix(vecs[j]).conjugate_transpose() #find the dot product of the vectors in the orthogonal representation
            if abs(b[0])<1e-9: #if that's about zero
                vi=ordering[i]; #name vi and vj as the ones in the order
                vj=ordering[j];
                G_r.add_edge([vi,vj]) #add an edge between vi and vj
    return G_r

##################################################

def lib_degen_OR_recon(G): #puts previous guys together:
    """
    Finds degeneracy of graph and optimal vertex-ordering, constructs orthogonal representation, and reconstructs the graph
    
    Args:
    - The graph
    
    Returns:
    - A graph reconstructed from vector adjacencies
    """
    
    [ordering, degen] = lib_degen_find(G)
    vecs = lib_modified_LSS(G,ordering,degen+1)
    G_r = lib_reconstructing(vecs,ordering)
    
    return G_r

##################################################