from sage.all import *

from .matrix_functions import *
from .graph_functions import *
from .orthogonal_representations_functions import *
from .quantum_coloring_functions import *
from .quantum_algorithm_functions import *

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister, ClassicalRegister
from qiskit.circuit.library import QFTGate

import numpy as np

##################################################

def libp_spectral_lower_bounds_part1(G):
    """
    Takes the inequality from Spectral Lower Bounds Part 1 and applies it to a graph.

    Args:
    - A graph

    Returns:
    - Lower bound for quantum chromatic number of the graph
    """
    n = len(G.vertices())
    m = len(G.edges())
    
    A = G.adjacency_matrix()
    D = diagonal_matrix(vector(G.degree()))
    Q = D + A #signless Laplacian of G
    L = D - A #Laplacian of G
    
    mu1 = max(A.eigenvalues()); mun = min(A.eigenvalues()) #max and min eigenvalues of A
    delt1 = max(Q.eigenvalues()); deltn = min(Q.eigenvalues()) #max and min eigenvalues of Q
    thet1 = max(L.eigenvalues()); thetn = min(L.eigenvalues()) #max and min eigenvalues of L
    
    eigpos = [j for j in A.eigenvalues() if j > 0] #positive eig vals of A
    eigneg = [j for j in A.eigenvalues() if j < 0] #negative eig vals of A
    
    nplus = len(eigpos); nneg = len(eigneg) #number of positive and negative eigenvalues of A
    
    splus = sum([j**2 for j in eigpos]); sneg = sum([j**2 for j in eigneg]) #sum of squares of pos and neg eig vals of A
    
    ## VARIOUS BOUNDS ON CHROMATIC NUMBER IN EQN. (3) ##
    bound_1 = 1 + mu1 / abs(mun);                   #Hoffman [8]
    bound_2 = 1 + 2*m / (2* m - n*deltn);           #Lima et al[11]
    bound_3 = 1 + mu1 / (mu1 - delt1 + thet1);      #Kolotilina [10]
    bound_4 = 1 + max(nplus / nneg, nneg / nplus);  #Elphick and Wocjan [7]
    bound_5 = 1 + max(splus / sneg, sneg / splus);  #Ando and Lin [1]
    
    low_bound = max(bound_1, bound_2, bound_3, bound_4, bound_5); #Implementing Eqn. (3)
    
    ## TK TO DO: find out how to make the lovasz theta function part work. Section 6 of the paper.
    ### b1 \leq vector chromatic number (G) \leq lovasz theta function (G complement) \leq quantum chromatic number (G)
    return low_bound, bound_1, bound_2, bound_3, bound_4, bound_5

##################################################

def libp_spectral_lower_bounds_part2(G):
    """
    Takes the inequality from Spectral Lower Bounds Part 2 and applies it to a graph.

    Args:
    - A graph

    Returns:
    - Lower bound for quantum chromatic number of the graph
    """
    
    A = G.adjacency_matrix()
    spec_A = A.eigenvalues() #eig vals of A
    spec_A_inc = spec_A.copy(); spec_A_inc.sort() #eig vals from lowest to highest
    spec_A_dec = spec_A.copy(); spec_A_dec.sort(reverse=1) #eig vals from highest to lowest
    
    spec_A_unique = list(set(spec_A)); spec_A_unique.sort(reverse = 1); #eigs of A w/o multiplicities from big to small
    
    umax = max(spec_A) #largest eigenvalue of A
    u2 = spec_A_unique[1] #second largest eigenvalue
    un = min(spec_A) #smallest eigenvalue
    
    ## FIRST BOUND ##
    summ = umax #implementing the inequality in Eqn. (20)
    for j in range(len(spec_A)):
        if summ <= 0:
            k=j
            break
        else:
            summ = summ + spec_A_inc[j]
        
    bound_1 = 1 + k #the lower bound in Eqn. (21)
    
    
    ## SECOND BOUND ##
    for (l,E) in A.right_eigenspaces(): #finding multiplicity of smallest eig val
        if l == un:
            g = E.dimension() #multiplicity of smallest eig val!
            break

    bound_2 = 1 + min(g,abs(un)/u2) #the lower bound in Eqn. (34)
    
    lower_bound = max(bound_1, bound_2)

    return lower_bound,bound_1, bound_2

##################################################

def libp_all_lower_bounds(G):
    """
    Gives information on the lower bounds of the quantum chromatic number of the graph, as well as its chromatic number
    
    Args:
    - A graph
    
    Returns:
    - Print statement of chromatic number
    - Print statement of clique number
    - Print statement of spectral lower bounds part 1
    - Print statement of spectral lower bounds part 2
    """
    
    sp1 = [float(j) for j in libp_spectral_lower_bounds_part1(G)]
    sp2 = [float(j) for j in libp_spectral_lower_bounds_part2(G)]
    
    print('Chromatic number:\t\t',G.chromatic_number(),'\n')
    
    print('Lower bound:\t\t\t',max(G.clique_number(),max(sp1),max(sp2)),'\n')
    
    print('Clique number:\t\t\t',G.clique_number(),'\n')
    
    print('Spectral bounds 1:\t\t',ceil(sp1[0]))
    print('\tHoffman\t\t\t', sp1[1])
    print('\tLima\t\t\t', sp1[2])
    print('\tKolotilina\t\t', sp1[3])
    print('\tElphick and Wocjan\t', sp1[4])
    print('\tAndo and Lin\t\t', sp1[5])
    
    print('\nSpectral bounds 2:\t\t', ceil(sp2[0]))
    print('\tBound 1:\t\t',sp2[1])
    print('\tBound 2:\t\t',sp2[2])
    
    return

##################################################

def libp_avis_QFTN(k0):
    """
    Performs the quantum Fourier transform in dimension N, given in equation (1) of Avis et al.
    
    Args:
    -A tall skinny vector
    
    Returns:
    -The quantum Fourier transform of that vector
    """
    k1 = Matrix(k0) #if it's not a matrix, it is now
    if k1.nrows() < k1.ncols(): #if it's short fat, now it's tall skinny
        k = k1.transpose()
    else:
        k = k1
    N=k.nrows() #how long is the vector
    QFTN = lib_make_Fourier_Matrix(N)*sqrt(N) #make a Fourier matrix of that same dimension
    return QFTN*k#QFTN multiplies the Fourier matrix by the vector to do the thing

##################################################

def libp_avis_phaseshift(L0):
    """
    Performs the phaseshift given in equation (2) of Avis et al.
    
    Args:
    -A tall skinny vector
    
    Returns:
    -The phaseshifted vector
    """
    L = Matrix(L0) #if it's not a matrix, it is now
    N = max(L.nrows(),L.ncols()) #length of vector
    PL=[]
    for j in range(N):
        PL.append((-1)**(L[j])*L[j])
    return PL

##################################################

def libp_avis_CNOT(a0, b0):
    """
    Performs CNOT operations on target vector (b0) using control vector (a0), as given in equation (4) of Avis et al.
    
    Args:
    -Target vector b0
    -Control vector a0
    
    Returns:
    -CNOT'd version of the b0 vector
    """
    
    a1 = Matrix(a0) #if it's not a matrix, it is now
    b1 = Matrix(b0)
    
    if a1.nrows() < a1.ncols():a = a1.transpose() #if it's short fat, now it's tall skinny
    else: a = a1
    if b1.nrows() < b1.ncols():b = b1.transpose()
    else: b = b1
        
    N = a.nrows()
    CNOT_B = [None for i in range(N)];
    for j in range(N):
        CNOT_B[j] = (RR(a[j][0])+RR(b[j][0]))%2
    return Matrix(CNOT_B).transpose()

##################################################

def libp_avis_protocol(G,a,b):
    #takes in Hadamard graph and query vertices from ref
    N = log(len(G.vertices()),2)
    a = Matrix(a)
    b = Matrix(b)
    c = (a+b)%2
    S0 = lib_make_standard_basis_vectors(N)
    
    #Step 3
    recm=[[]for i in range(N)] #empty list to store things in
    for jA in range(N):
        for jB in range(N):
            summ3 = 0
            for k in range(N):
                summ3 = summ3 + exp(2*pi*i/N)**(k*(jA-jB))*(-1)**(c[0][k])
            summ3 = CC((1/(sqrt(N)))**3*summ3)
            recm[jA].append(summ3.norm())
    
    #for alice:
    alicer0 = [sum(recm[j])for j in range(N)]
    alicer = []
    for j in range(N):
        alicer.append(sum(alicer0[k] for k in range(j+1)))
    
    randa = uniform(0,1) #random uniform number between 0 and 1
    for j in range(N):
        if randa <= alicer[j]:
            alice = S0[j]
            alices = j
            break

    #for bob:
    bobr = []
    for j in range(N):
        bobr.append(sum(recm[alices][k] for k in range(j+1)))

    randb = uniform(0,sum(recm[alices]))
    for j in range(N):
        if randb <= bobr[j]:
            bob = S0[j]
            break
    
    return alice, bob,recm

##################################################

def libp_avis_game(G,numit):
    #takes in Hadamard graph and plays the game desired number of times
    N = len(G.vertices())
    it = 0
    bads = 0
    goods = 0
    while it < numit:
        ei = randrange(0,N)
        a = Matrix(G.edges()[ei][randrange(0,2)])
        b = Matrix(G.edges()[ei][randrange(0,2)])

        [alice, bob,recm] = libp_avis_protocol(G,a,b)

        if ((a+b)%2 == 0):
            if alice != bob:
                bads += 1
                print('FAILURE: same vertex',ei)
                print('alice:',alice,'\nbob:',bob)
                print(recm)
                break
            else:
                goods += 1
        if (a+b)%2 != 0:
            if alice == bob:
                bads += 1
                print('FAILURE: edge',ei)
                print('alice:',alice,'\nbob:',bob)
                print(recm)
                break
            else:
                goods += 1
        it += 1
    return bads, goods

##################################################

def libp_quantum_avis_algorithm(va,vb): #for Hadamard graph G_{4N}
    """
    Runs the algorithm from Avis et al. on vertices given to Alice and Bob, and outputs the related quantum circuit.

    Args:
    - va: Alice's vertex as a vector (list format)
    - vb: Bob's vertex as a vector (list format)

    Returns:
    - The circuit qc that is created by the steps in the paper.
    - A list of 3-item lists: Alice's answer, Bob's answer, and the number of times that occurrence happened in the simulation.
    """
    def make_list(L,n1):
        bitslist = []
        for j1 in range(L):
            bitslist.append(format(j1,'0'+str(n1)+'b'))
        return bitslist
    
    def Hadamard_game_measure_answers(c2):
        #takes results from the simulation and puts them in form of Alice and Bob's answers
        c3=[]
        for j in range(len(c2)):
            c3.append([int(c2[j][0][:n],2),int(c2[j][0][n:],2),c2[j][1]])
        return c3

    N=len(va) #N will be the length of the vector associated with a vertex
    n = int(np.log2(N)) #and big N = 2^n so find little n
    
    vc=[((va[j1]+vb[j1])%2) for j1 in range(N)] #add the two vectors together mod 2
    

    player_reg = QuantumRegister(2*n,name="player bits") #quantum register for the players (first half are Alice's, second half are Bob's)
    qc = QuantumCircuit(player_reg,name = "Avis Game") #set up the circuit

    qft_a = QFTGate(n) #set up the Quantum Fourier Transform as a gate on n dimensions

    #Step 1: Prepare initial state
    qc.append(qft_a,player_reg[:n]) #Alice applies qft to her qubits

    for j1 in range(n//2+1): #Bob applies CNOT gate to his qubits with Alice's as the control
        qc.cx(player_reg[j1],player_reg[j1+n])

    qc.barrier() #visual separation
    statelist = make_list(N,n)
    #Step 2: Apply phaseshifts
    for j1 in range(len(vc)):
        if va[j1] == 1:
            for j2 in range(len(statelist[j1])):
                if statelist[j1][j2] == '0':
                    qc.x(player_reg[j2])
            # qc.cz(player_reg[0],player_reg[1])
            qc.h(player_reg[n-1])
            qc.mcx(player_reg[0:n-1],player_reg[n-1])
            qc.h(player_reg[n-1])
            for j2 in range(len(statelist[j1])):
                if statelist[j1][j2] == '0':
                    qc.x(player_reg[j2])
        
        if vb[j1] == 1:
            for j2 in range(len(statelist[j1])):
                if statelist[j1][j2] == '0':
                    qc.x(player_reg[j2+n])
            # qc.cz([player_reg[j] for j in range(n+1,N)])
            qc.h(player_reg[N-1])
            qc.mcx(list(range(n,N-1)),player_reg[N-1])
            # qc.cz(player_reg[2],player_reg[3])
            qc.h(player_reg[N-1])
            for j2 in range(len(statelist[j1])):
                if statelist[j1][j2] == '0':
                    qc.x(player_reg[j2+n])

    #Step 3: Another Fourier transform
    qc.append(qft_a,player_reg[:n]) #Alice applies qft to her qubits
    qc.append(qft_a.inverse(),player_reg[n:]) #bob applies inverse qft to his qubits

    #Step 4: Everyone measures.
    results = Hadamard_game_measure_answers(libq_simulate(qc))

    return qc, results