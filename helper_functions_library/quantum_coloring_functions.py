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

def lib_quantum_coloring_to_gameplay(qch,a,b):
    """
    Takes a quantum coloring of a graph (in matrix form and all) and plays the quantum coloring "game" with Alice and Bob and a ref.

    Args:
    - G: A graph.
    - qch: A quantum coloring of that graph.
    - vecs: List of vectors for orthogonal representation of the graph.
    - a: Query vertex for Alice.
    - b: Query vertex for Bob.

    Returns:
    - Alice's answer color.
    - Bob's answer color.
    """

    #Gets the measurements associated with the query vertices.
    E = qch[a] #The set of POVMs associated to Alice's vertex
    F = qch[b] #And for Bob

    d = len(qch[0][0][0]) #dimension of coloring
    c = len(qch[0]) #quantum chromatic number

    #The "maximally entangled state" is the normalized all ones vector. Calling this |ψ>
    S0 = lib_make_standard_basis_vectors(d)
    st = sum(S0)/(sum(S0)).norm()

    
    #Measure probabilities of outcomes. For each measure E_m in {E}, apply that to the entangled state.
    #Then you get that the probability of getting outcome m when measuring the given state |ψ> is <ψ|E_m|ψ>.
    #Measuring a quantum state inherently changes it, and that new state is (E_m|ψ>) / (sqrt(<ψ|E_m|ψ>)).
    
    #We'll first have Alice measure the entangled state with her POVMs, and gather a list of the possible outcome states, as well as their probabilities. Then we'll generate a random number and use it to decide which of the possible outcome states is the one that "actually happened" (since we don't have a quantum computer, this is the best I can do right now). The m associated with the E_m|ψ> that "actually happened" is Alice's answer.
    #Then we'll take that outcome state E_m|ψ> and have Bob measure it with his POVMs. Similar to Alice's case, we'll collect lists of the possible outcome states and their probabilities, and use a generated random number to decide which of those states is the one that "actually happened". The m associated with the F_m|ψ> that "actually happened" is Bob's answer.

    #For all m, gather the list of E_m|ψ> and <ψ|E_m|ψ>
    newst1s = [] 
    pr_newst1s = []
    for j in range(c):
        a_pr = st.transpose()*E[j]*st 
        newst = E[j]*st
        pr_newst1s.append(a_pr)
        newst1s.append(newst)
    

    #Decide which state we want (randomly)
    #Cumulative probabilities for Alice and Bob's states. This is not physics or anything, this is me not really knowing how to code.
    rand_a = uniform(0,1)
    a_cum_probs = []
    for j in range(4):
        a_cum_probs.append(sum(pr_newst1s[k] for k in range(j+1)))
    for j in range(4):
        if rand_a <= a_cum_probs[j]:
            newst1 = (newst1s[j])/(newst1s[j].norm())
            alice_ans = j
            del a_cum_probs
            break

    #That chosen newst1 = E_m|ψ> is Alice and Bob's entangled state after Alice applies her unitary measurement to it.
    #Now Bob will apply his unitary measurement to that new state to get his answer.
    newst2s = []
    pr_newst2s = []
    for j in range(4):        
        b_pr = newst1.transpose()*F[j]*newst1   
        newst2 = F[j]*newst1
        
        pr_newst2s.append(b_pr)
        newst2s.append(newst2)

    #That thing again where we make a random number and use it to choose Bob's state, but I can't code so I have to do it like this.
    rand_b = uniform(0,1)
    b_cum_probs = []
    for j in range(4):
        b_cum_probs.append(sum(pr_newst2s[k] for k in range(j+1)))
    for j in range(4):
        if rand_b <= b_cum_probs[j]:
            newst2 = (newst2s[j])/(newst2s[j].norm())
            bob_ans = j
            break
    
    return alice_ans,bob_ans

##################################################