from sage.all import *

from .matrix_functions import *
from .graph_functions import *
from .orthogonal_representations_functions import *
from .quantum_coloring_functions import *

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister, ClassicalRegister
from qiskit.quantum_info import Statevector, Operator, partial_trace
from qiskit.circuit.library.standard_gates import RYGate, RZGate
from qiskit.circuit.library import QFTGate

from qiskit_aer import AerSimulator
from qiskit import transpile

##################################################

def libq_simulate(circuit):
    """
    Runs an experiment on a circuit to determine outputs.
    
    Args:
    - circuit: The quantum circuit to be experimented on.
    
    Returns:
    - A list of pairs: the output states and the times they were simulated.
    """

    simulator = AerSimulator()
    # Transpile the circuit for the backend
    compiled_circuit = transpile(circuit, simulator)

    # Run the circuit
    job = simulator.run(compiled_circuit, shots=1000)

    # Get the measurement counts
    counts = job.result().get_counts()
    
    k=list(counts.keys())
    v=list(counts.values())
    v=[v[j]/1000 for j in range(len(v))]
    c2=[]
    for j in range(len(v)):
        m=v.index(max(v))
        c2.append([k[m],max(v)])
        v.remove(max(v))
        k.remove(k[m])
    return c2

##################################################

def libq_results_from_statevector(output_state):
    """
    Takes a Statevector item and turns it into a list of lists of the states and their probabilities.
    
    Args:
    - output_state: A Statevector item of a quantum circuit
    
    Returns:
    - List of pairs: a possible output state, and the probability of it occurring.
    """

    t=output_state.to_dict()

    #get the keys
    k=list(t.keys())
    v=list(t.values())
    L=[]
    v1=[] #make the 'results' into probabilities so they're useful
    for j in range(len(v)):
        L.append([k[j],v[j]*v[j].conjugate()])
        v1.append(v[j]*v[j].conjugate())
    return L

##################################################

def libq_quantum_avis_algorithm(va,vb): #for Hadamard graph G_{4N}
    """
    Runs the algorithm from Avis et al. on vertices given to Alice and Bob, and outputs the related quantum circuit.

    Args:
    - va: Alice's vertex as a vector (list format)
    - vb: Bob's vertex as a vector (list format)

    Returns:
    - The circuit U that is created by the steps in the paper.
    - A list of 3-item lists: Alice's answer, Bob's answer, and the number of times that occurrence happened in the simulation.
    """

    N=len(va) #N will be the length of the vector associated with a vertex
    n = int(np.log2(N)) #and big N = 2^n so find little n
    
    vc=[((va[j1]+vb[j1])%2) for j1 in range(N)] #add the two vectors together mod 2
    

    player_reg = QuantumRegister(2*n,name="player bits") #quantum register for the players (first half are Alice's, second half are Bob's)
    U = QuantumCircuit(player_reg,name = "Avis Game") #set up the circuit

    qft_a = QFTGate(n) #set up the Quantum Fourier Transform as a gate on n dimensions

    #Step 1: Prepare initial state
    U.append(qft_a,player_reg[:n]) #Alice applies qft to her qubits

    for j1 in range(n//2+1): #Bob applies CNOT gate to his qubits with Alice's as the control
        U.cx(player_reg[j1],player_reg[j1+n])

    U.barrier()

    #Step 2: Apply phaseshifts
    for j1 in range(len(vc)): #For each entry in the vc vector, add a Z gate to the circuit if the entry is 1.
        if vc[j1] == 1:
            U.z(player_reg[j1])

    #Step 3: Another Fourier transform
    U.append(qft_a,player_reg[:n]) #Alice applies qft to her qubits
    U.append(qft_a.inverse(),player_reg[n:]) #bob applies inverse qft to his qubits

    #Step 4: Everyone measures
    U.measure_all(inplace=True)

    #Step Me: I get the measurements and the answers
    c2=libq_simulate(U)

    c3=[]
    for j in range(len(c2)):
        c3.append([int(c2[j][0][:n],2),int(c2[j][0][n:],2),c2[j][1]])

    return U, c3

##################################################