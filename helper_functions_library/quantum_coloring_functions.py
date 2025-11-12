from sage.all import *
from .matrix_functions import *
from .graph_functions import *
from .orthogonal_representations_functions import *

##################################################

def lib_quat_tophats(vec):
    """
    Makes the four tophat vectors for quaternion construction from one 4D vector, and the associated projectors
    
    Args:
    - A 4D vector
    
    Returns:
    - A list of four orthogonal projectors
    """
    
    phi0 = Matrix(vec[0])
    phi1 = Matrix([-vec[0][1],vec[0][0],-vec[0][3],vec[0][2]])
    phi2 = Matrix([-vec[0][2],vec[0][3],vec[0][0],-vec[0][1]])
    phi3 = Matrix([-vec[0][3],-vec[0][2],vec[0][1],vec[0][0]])

    #now outerproduct dotify them to make projection matrices
    proj0 = phi0.T * phi0
    proj1 = phi1.T * phi1
    proj2 = phi2.T * phi2
    proj3 = phi3.T * phi3
    
    return proj0, proj1, proj2, proj3 #tadah
    
##################################################
    
def lib_oct_tophats(vec):
    """
    Makes the eight tophat vectors for octonion construction from one 8D vector, and the associated projectors
    
    Args:
    - An 8D vector
    
    Returns:
    - A list of eight orthogonal projectors
    """
        
    tophat=[None for i in range(8)];
    tophat[0] = Matrix([ vec[0][0],  vec[0][1],  vec[0][2],  vec[0][3],  vec[0][4],  vec[0][5],  vec[0][6],  vec[0][7] ])
    tophat[1] = Matrix([-vec[0][1],  vec[0][0], -vec[0][3],  vec[0][2], -vec[0][5],  vec[0][4],  vec[0][7], -vec[0][6] ])
    tophat[2] = Matrix([-vec[0][2],  vec[0][3],  vec[0][0], -vec[0][1], -vec[0][6], -vec[0][7],  vec[0][4],  vec[0][5] ])
    tophat[3] = Matrix([-vec[0][3], -vec[0][2],  vec[0][1],  vec[0][0], -vec[0][7],  vec[0][6], -vec[0][5],  vec[0][4] ])
    tophat[4] = Matrix([-vec[0][4],  vec[0][5],  vec[0][6],  vec[0][7],  vec[0][0], -vec[0][1], -vec[0][2], -vec[0][3] ])
    tophat[5] = Matrix([-vec[0][5], -vec[0][4],  vec[0][7], -vec[0][6],  vec[0][1],  vec[0][0],  vec[0][3], -vec[0][2] ])
    tophat[6] = Matrix([-vec[0][6], -vec[0][7], -vec[0][4],  vec[0][5],  vec[0][2], -vec[0][3],  vec[0][0],  vec[0][1] ])
    tophat[7] = Matrix([-vec[0][7],  vec[0][6], -vec[0][5], -vec[0][4],  vec[0][3],  vec[0][2], -vec[0][1],  vec[0][0] ])
    
    cowboy=[];
    for i in range(8):
        cowboy.append(tophat[i].conjugate_transpose() * tophat[i])
    
    return cowboy

##################################################

def lib_quat_const(vecs):
    """
    Constructs quantum coloring using quaternion construction
    
    Args:
    - List of 4D vectors in a real orthogonal representation
    
    Returns:
    - A quantum 4-coloring
    """
    
    n = len(vecs)
    qch = [];
    
    for j in range(n):
        qch.append(lib_quat_tophats(vecs[j]/vecs[j].norm()))
    return qch

##################################################

def lib_oct_const(vecs):
    """
    Constructs quantum coloring using octonion construction
    
    Args:
    - List of 8D vectors in a real orthogonal representation
    
    Returns:
    - A quantum 8-coloring
    """
        
    n=len(vecs)
    qch = [];
    for j in range(n):
         qch.append(lib_oct_tophats(vecs[j]/vecs[j].norm()))
    return qch
    
##################################################

def lib_check_qch(qch,G,ordering):
    """
    Checks whether or not it's actually a quantum coloring
    
    Args:
    - Purported quantum coloring
    - The graph
    
    Returns:
    - IF FAILED, returns errors
    """
    
    n = len(qch);
    c = len(qch[0])
    d = len(qch[0][0][0])
    
    comp_errs = 0; #initialize number of completeness errors
    orth_errs = 0; #initialize number of orthogonality errors
    #checking completeness
    for vi in range(n):
        summ = sum(qch[vi]) #sum them all up
        for i in range(d):
            for j in range(d):
                if ((i == j) and (abs(summ[i,j]) >= 1+1e-5)) or ((i != j) and (abs(summ[i][j] > 1e-5))):
                    print('Failure on vertex ', G.vertices()[vi]) #if ON diagonal and NOT ~1, or OFF diagonal and NOT ~0, print where it failed
                    comp_errs +=1;

    #checking orthogonality
    #if they're different vertices and their respective qch matrices multiply to have trace not ~0, and they're adjacent
    #then print where they failed
    for i in range(n):
        vi = ordering[i]
        for j in range(i):
            vj = ordering[j];
            for ci in range(c):
                #if the vertices are different AND the trace of their product isn't zero AND they're neighbors, then it's failure
                #if ((vi != vj) and (abs((qch[vi][ci]*qch[vj][ci].conjugate_transpose()).trace()) > 1e-5) and (G.vertices()[vi] in G.neighbors(G.vertices()[vj]))):
                #if (  abs( (qch[i][ci]*qch[j][ci].conjugate_transpose() ).trace() ) > 1e-5  ) and ( vi in G.neighbors(vj)):
                vjMAT = qch[j][ci]
                viMAT = qch[i][ci]
                thetrace = (vjMAT*viMAT).trace()
                if (vi in G.neighbors(vj)) and (abs(thetrace) > 1e-5) :
                    print('Failure on vertices ', G.vertices()[vi], G.vertices()[vj],' on color ',ci)
                    print('With trace of ', thetrace)
                    #print('With trace of ',((qch[vi][ci]*qch[vj][ci].conjugate_transpose()).trace()))
                    #print('And matrices of ',((qch[vi][ci]*qch[vj][ci].conjugate_transpose())))
                    orth_errs +=1;
    
    if ((comp_errs == 0) and (orth_errs == 0)):
        print('All quiet on the errors front')
    else:
        print('Wow, whoops!\nCompleteness errors:', comp_errs,'\nOrthogonality errors: ',orth_errs)

    return

##################################################