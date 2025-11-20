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

import copy

from qiskit.circuit.library.standard_gates import RYGate, RZGate, CPhaseGate
import numpy as np

##################################################

def libq_prepare_quantum_state(psi):
    """
    Takes a vector |ψ> with complex entries and builds a circuit that will prepare that state.
    
    Args:
    - psi: A vector input of a desired state |ψ> (with dimension 2^n for some n)
    
    Returns:
    - U_circ: The quantum circuit that prepares the state |ψ>
    """
    def make_list(L,n1):
        bitslist = []
        for j1 in range(L):
            bitslist.append(format(j1,'0'+str(n1)+'b'))
        return bitslist
    
    def get_all_angles(state):
        """
        Iteratively build a list of angles of rotation.
        #Get the angles of rotation to be used lastly in the iterative building process, as well as the coefficients of the state |ψ>_{n-1} (the state we must prepare in order to iteratively build |ψ>_n)
        
        Args:
        - state: A vector of coefficients of the state.
        
        Returns:
        - A list of all the angles of rotation for the circuit.
        """
        n = int(np.log2(len(state)))
        def get_angles(state_0): #Finds the rotation angles of a given state. Outputs a list of angles as well as a new state
            thetas = []
            phis = []
            gammas = []
            newstate = []
            for j1 in range(len(state_0)//2):
                globalphase = (np.conj(state_0[2*j1])) / (np.linalg.norm(state_0[2*j1]))
                theta = 2*np.arccos(globalphase*state_0[2*j1]/np.sqrt(np.linalg.norm(state_0[2*j1])**2+np.linalg.norm(state_0[2*j1+1])**2))


                gamma = np.log(globalphase)/1j
                phi = np.log((globalphase*state_0[2*j1+1])/np.sin(theta/2))/1j

                if np.imag(theta) <= 1e-5:
                    thetas.append(np.real(theta))

                phis.append(np.real(phi))

                if np.imag(gamma) <= 1e-5:
                    gammas.append(np.real(2*gamma))
                
                newstate.append(np.sqrt(np.linalg.norm(state_0[2*j1])**2+np.linalg.norm(state_0[2*j1+1])**2))
            return thetas, phis, gammas, newstate

        tot_thetas = []; tot_phis = []; tot_gammas = []; #Start with empty lists
        [thetas, phis, gammas, newstate] = get_angles(state)
        tot_thetas.append(thetas); tot_phis.append(phis); tot_gammas.append(gammas); #Add these angles to the list
        for j1 in range(n-1):
            #Get the angles of rotation to be used on the new state, as well as the coefficients of the next state down in the iterative process
            [thetas, phis, gammas, newstate] = get_angles(newstate)
            tot_thetas.append(thetas); tot_phis.append(phis); tot_gammas.append(gammas)
        
        #We will perform this iterative process starting from |ψ>_1 and going up to state |ψ>_n, so we must reverse the order of the angles in the list.
        tot_thetas = tot_thetas[::-1]; tot_phis = tot_phis[::-1]; tot_gammas = tot_gammas[::-1]
        return tot_thetas, tot_phis, tot_gammas

    L = len(psi)
    n = int(np.log2(L))

    #Get the angles and phases to be used in the circuit
    [thetas,phis,gammas] = get_all_angles(psi)
    #Make a list of the basis states for n-1,...,1 (example for n = 3: [[00, 01, 10, 11] , [0, 1], [0]]). This is how we will determine where to put X gates in the circuit
    basis_states = []
    for j1 in range(1,n+1):
        basis_states.append(make_list(2**(n-j1),n-j1))

    #Again, since we are performing this iterative process starting with |ψ>_1 and going up to state |ψ>_n, we reverse the order of this list.
    basis_states=list(reversed(basis_states))

    #Initialize the circuit.
    q_reg = QuantumRegister(n, name = "x") #We need as many bits as there are in |ψ>_n
    U_circ = QuantumCircuit(q_reg, name = "Quantum State Preparation")

    #Perform the first Y rotation on the first qubit
    U_circ.ry(thetas[0][0],q_reg[n-1])    

    #For the remainder of the rotations, we will go through the list of angles and apply them to the circuit with a controlled Y rotation gate.
    for j1 in range(1,len(thetas)):
        curr_thetas = thetas[j1]
        curr_phis = phis[j1]
        curr_gammas = gammas[j1]
        for j2 in range(len(curr_thetas)):
            now_theta = curr_thetas[j2]
            now_phi = curr_phis[j2]
            now_gamma = curr_gammas[j2]

            #If the state corresponding to this angle has a zero qubit in it, we apply an X gate to that qubit both before and after the rotation, to ensure that the Y rotation applies correctly.
            for j3 in range(len(basis_states[j1][j2])):
                if basis_states[j1][j2][j3] == '0':
                    U_circ.x(q_reg[n-1-j3])

            
            #Multi-controlled Y gate, Z, gate, and two phase gates
            U_circ.append(RYGate(now_theta).control(j1),q_reg[::-1][0:j1+1]) #puts theta on both |0> and |1>
            U_circ.append(RZGate(now_gamma).control(j1),q_reg[::-1][0:j1+1]) #puts global phase shift on both |0> and |1>
            U_circ.append(CPhaseGate(now_phi).control(j1-1),q_reg[::-1][0:j1+1]) #puts phi phase on |1>
            U_circ.append(CPhaseGate(-now_gamma).control(j1-1),q_reg[::-1][0:j1+1]) #takes global phase off of |1>

            #Uncompute previously applied X gates
            for j3 in range(len(basis_states[j1][j2])):
                if basis_states[j1][j2][j3] == '0':
                    U_circ.x(q_reg[n-1-j3])
            U_circ.barrier()

    return U_circ

##################################################
##################################################

def libq_simulate(inputcirc):
    """
    Runs an experiment on a circuit to determine outputs.
    
    Args:
    - circuit: The quantum circuit to be experimented on.
    
    Returns:
    - A list of pairs: the output states and the times they were simulated.
    """
    circuit = copy.deepcopy(inputcirc)
    circuit.measure_all(inplace=True)
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