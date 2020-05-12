# Importing standard Qiskit libraries and configuring account
from qiskit import QuantumCircuit, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
from qiskit.tools.jupyter import *
from qiskit.visualization import *
import numpy as np
from qiskit import compiler
from qiskit.circuit import QuantumCircuit
from qiskit.transpiler.passes import Unroller
from qiskit.transpiler import PassManager
from qiskit import BasicAer
from qiskit.providers.aer import noise
import random
from qiskit import QuantumCircuit
from qiskit import QuantumRegister
from qiskit import ClassicalRegister
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import gridspec
from qiskit.providers.aer.noise import NoiseModel
import itertools
import scipy.special
from scipy.special import binom as binomial
from itertools import permutations
import random 

####################################################################
### Utility functions ###
def apply_extra_cnots(circuit, a, b, num_pairs):
    circuit.cx(a, b)
    for n in range(num_pairs):
        circuit.cx(a, b)
        circuit.cx(a, b)
####################################################################
def apply_u1(circuit, lam, t):
    circuit.u1(lam, t)
####################################################################
def apply_u2(circuit, phi, lam, t):
    circuit.u2(phi, lam, t)
####################################################################
def apply_u3(circuit, theta, phi, lam, t):
    circuit.u3(theta, phi, lam, t)
####################################################################
def generate_combinations(num_qubits):
    result = []
    for n in range(num_qubits):
        for m in range(num_qubits):
            if n!=m:
                result.append([n,m])
    return result
####################################################################
def kbits(n, k):
    result = []
    for bits in itertools.combinations(range(n), k):
        s = ['0'] * n
        for bit in bits:
            s[bit] = '1'
        result.append(''.join(s))
    return result
####################################################################
def place_permutation(ref_str, perm):
    result = []
    count = 0
    for s in ref_str:
        if s == '1':
            result.append(perm[count])
            count += 1
        else:
            result.append(0)
    return result
####################################################################
def generate_cx_error_noise_model(num_qubits, param):
    noise_model = noise.noise_model.NoiseModel()
    pairs = generate_combinations(num_qubits)
    for pair in pairs:
        err = noise.errors.standard_errors.depolarizing_error(param, 2, True)
        noise_model.add_quantum_error(err, "cx", pair)
    return noise_model
####################################################################
def generate_combinations(num_qubits):
    result = []
    for n in range(num_qubits):
        for m in range(num_qubits):
            if n!=m:
                result.append([n,m])
    return result
####################################################################
def p_to_lam(p):
    return p*((4**2)/(4**2-1))

####################################################################
def q1_q2_obs(counts, shots):
    total = 0
    for key in counts.keys():
        total +=  counts[key]*int(key, 2)
    return total/shots

####################################################################
def FIIM_generate_circs_helper(circuit, n, basis_gates=['u1', 'u2', 'u3', 'cx']):
    unroller = Unroller(basis=basis_gates)  
    p_m = PassManager(passes=[unroller])
    transpile_result = compiler.transpile(
        circuit,
        basis_gates=basis_gates,
        optimization_level=1
    )
    ops = transpile_result.data

    qc = QuantumCircuit()
    qregs = transpile_result.qregs
    qubits = []
    for qreg in qregs:
        if not qc.has_register(qreg):
            qc.add_register(qreg)
        qubits.extend(qreg)
    cregs = circuit.cregs
    clbits = []
    for creg in cregs:
        if not qc.has_register(creg):
            qc.add_register(creg)
        clbits.extend(creg)
    for op in ops:
        if op[0].name == 'cx':
            apply_extra_cnots(qc, *op[1], n)
        elif op[0].name == 'u1':
            apply_u1(qc, *op[0].params, op[1][0])
        # elif op[0].name == 'h':
        #    apply_h(qc, *op[0].params, op[1][0])
        elif op[0].name == 'u2':
            apply_u2(qc, *op[0].params, op[1][0])
        elif op[0].name == 'u3':
            apply_u3(qc, *op[0].params, op[1][0])
        elif op[0].name == 'barrier':
            qc.barrier()
        elif op[0].name == 'measure':
            qc.measure(op[1][0], op[2][0])
        else:
            print("ERROR:", op[0].name)
    return qc

####################################################################
def FIIM_generate_circs(n, circ):
    """
    Inserts extra CNOT gate pairs into circuit for zero-noise extrapolation
    on every CNOT in the input circuit.
    :param circuit: (QuantumCircuit) The quantum circuit to add CNOTs to.
    :param num_pairs_per_cnot: (int) number of CNOTs pairs to add to each original CNOT in the circuit
    :return: (QuantumCircuit) new Quantum Circuit object which is equivalent to the original with extra CNOTs inserted.
    """
    FIIM_circs = []
    r_vals = [1 + 2 * i for i in range(n + 1)]
    for i in range(n + 1):
        FIIM_circs.append(FIIM_generate_circs_helper(circ, i))
    return FIIM_circs, r_vals

####################################################################
def RIIM_coeffs(n, num_cnots):
    """
    Generates coefficients a_n for RIIM extrapolation.
    :param n: (int) The order of RIIM to be performed.
    :param num_cnots: (int) number of CNOTs in the circuit.
    :return: (list) A list of RIIM coefficients
    """
    
    Nc = num_cnots

    a3_1 = -1 / 2
    a_1 = -Nc * a3_1 + 1

    a3_2 = -(Nc + 4) * (1 / 4)
    a5_2 = 3 / 8
    a33_2 = 1 / 4
    a_2 = 1 - Nc * a3_2 - Nc * a5_2 - binomial(Nc, 2) * a33_2

    a3_3 = -(Nc ** 2 + 10 * Nc + 24) * (1 / 16)
    a5_3 = (1 / 16) * (3 * (Nc + 6))
    a33_3 = (Nc + 6) * (1 / 8)
    a7_3 = -5 / 16
    a35_3 = -3 / 16
    a333_3 = -1 / 8
    a_3 = 1 - Nc * a3_3 - Nc * a5_3 - binomial(Nc, 2) * a33_3 - Nc * a7_3 - 2 * binomial(Nc, 2) * a35_3 - binomial(Nc,
                                                                                                                   3) * a333_3

    a3_4 = -(Nc ** 3 + 18 * Nc ** 2 + 104 * Nc + 192) * (1 / 96)
    a5_4 = (3 * Nc ** 2 + 32 * Nc + 154) * (1 / 64)
    a33_4 = (Nc ** 2 + 14 * Nc + 58) * (1 / 32)
    a7_4 = -45 / 32
    a35_4 = -(3 * Nc + 29) * (1 / 32)
    a333_4 = -(Nc + 8) * (1 / 16)
    a9_4 = 35 / 128
    a55_4 = 29 / 64
    a335_4 = 3 / 32
    a3333_4 = 1 / 16
    a73_4 = 0
    a_4 = 1 - Nc * a3_4 - Nc * a5_4 - Nc * a7_4 - Nc * a9_4 - (binomial(Nc, 4) * a3333_4) - (
                binomial(Nc, 3) * a333_4) - (3 * binomial(Nc, 3) * a335_4) - (binomial(Nc, 2) * a33_4) - (
                      2 * binomial(Nc, 2) * a35_4) - (binomial(Nc, 2) * a55_4)

    coeffs_1 = [
        a_1,
        a3_1
    ]

    coeffs_2 = [
        a_2,
        a3_2,
        a5_2,
        a33_2
    ]

    coeffs_3 = [
        a_3,
        a3_3,
        a5_3,
        a33_3,
        a7_3,
        a35_3,
        a333_3
    ]

    coeffs_4 = [
        a_4,
        a3_4,
        a5_4,
        a33_4,
        a7_4,
        a35_4,
        a333_4,
        a9_4,
        a73_4,
        a55_4,
        a335_4,
        a3333_4
    ]

    coeffs = [coeffs_1, coeffs_2, coeffs_3, coeffs_4]
    return coeffs[n - 1]

####################################################################
def RIIM_arrange(combination, num_cnots):
    positions = kbits(num_cnots, len(combination))
    perms = set(list(permutations(combination)))
    arrangements = []
    for p in positions: 
        for q in perms: 
            arrangements.append(place_permutation(p, q))
    return arrangements

####################################################################
def RIIM_arrange(combination, num_cnots):
    positions = kbits(num_cnots, len(combination))
    perms = set(list(permutations(combination)))
    arrangements = []
    for p in positions:  
        for q in perms:
            arrangements.append(place_permutation(p, q))
    return arrangements

####################################################################
def sample_from_list(n, input_list):
    max_index = len(input_list)
    output_list = []
    for i in range(n): 
        output_list.append(input_list[random.randint(0, max_index-1)]) 
    return output_list    

####################################################################
def RIIM_gate_combinations(n, num_cnots):
    combinations = [
        [[1]],
        [[2], [1, 1]],
        [[3], [2, 1], [1, 1, 1]],
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    ]
    all_combinations = []
    all_combinations.append([[0 for i in range(num_cnots)]])
    for i in range(n):
        for j in combinations[i]:
            all_combinations.append(RIIM_arrange(j, num_cnots))
    return all_combinations

####################################################################
def RIIM_gate_combinations_sampled(n, num_cnots, num_max): 
    combinations = [
        [[1]],
        [[2], [1, 1]],
        [[3], [2, 1], [1, 1, 1]],
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    ]
    all_combinations = []
    normalizations = []
    all_combinations.append([[0 for i in range(num_cnots)]]) 
    normalizations.append(1)
    for i in range(n):
        for j in combinations[i]:
            arranged = RIIM_arrange(j, num_cnots)
            all_combinations.append(sample_from_list(num_max, arranged))
            normalizations.append(len(arranged)/num_max)
    return all_combinations, normalizations

####################################################################
def RIIM_generate_circs(n, circuit, basis_gates=['u1', 'u2', 'u3', 'cx']):
    """
    Inserts extra CNOT gate pairs into circuit for RIIM.
    :param n: (int) order of RIIM to be performed 
    :param circuit: (QuantumCircuit) The quantum circuit to add CNOTs to.
    :return: (list) list of quantum circuits for RIIM extrapolation
    """
    unroller = Unroller(basis=basis_gates)
    p_m = PassManager(passes=[unroller])
    transpile_result = compiler.transpile(
        circuit,
        basis_gates=basis_gates,
        optimization_level=0
    )
    ops = transpile_result.data
    num_cnots = (circuit.count_ops()['cx'])
    combinations = RIIM_gate_combinations(n, num_cnots)
    new_circs = []
    for extra_set in combinations:
        tmp_circs = []
        for j in extra_set:
            qc = QuantumCircuit()
            qregs = transpile_result.qregs
            qubits = []
            for qreg in qregs:
                if not qc.has_register(qreg):
                    qc.add_register(qreg)
                qubits.extend(qreg)
            cregs = circuit.cregs
            clbits = []
            for creg in cregs:
                if not qc.has_register(creg):
                    qc.add_register(creg)
                clbits.extend(creg)
            cx_count = 0

            for op in ops:
                if op[0].name == 'cx':
                    if j[cx_count] == 0:
                        apply_extra_cnots(qc, *op[1], 0)
                        cx_count += 1
                    else:
                        apply_extra_cnots(qc, *op[1], j[cx_count])
                        cx_count += 1
                elif op[0].name == 'u1':
                    apply_u1(qc, *op[0].params, op[1][0])
                # elif op[0].name == 'h':
                #    apply_h(qc, *op[0].params, op[1][0])
                elif op[0].name == 'u2':
                    apply_u2(qc, *op[0].params, op[1][0])
                elif op[0].name == 'u3':
                    apply_u3(qc, *op[0].params, op[1][0])
                elif op[0].name == 'barrier':
                    qc.barrier()
                elif op[0].name == 'measure':
                    qc.measure(op[1][0], op[2][0])
                else:
                    print("ERROR:", op[0].name)
            tmp_circs.append(qc)
        new_circs.append(tmp_circs)
    return new_circs, RIIM_coeffs(n, num_cnots)

####################################################################
def RIIM_generate_circs_sampled(n, circuit, num_max, basis_gates=['u1', 'u2', 'u3', 'cx']):
    """
    Inserts extra CNOT gate pairs into circuit for RIIM, but with randomized sampling.
    :param n: (int) order of RIIM to be performed 
    :param circuit: (QuantumCircuit) The quantum circuit to add CNOTs to.
    :return: (list) list of quantum circuits for RIIM extrapolation
    """
    unroller = Unroller(basis=basis_gates)
    p_m = PassManager(passes=[unroller])
    transpile_result = compiler.transpile(
        circuit,
        basis_gates=basis_gates,
        optimization_level=0
    )
    ops = transpile_result.data
    num_cnots = (circuit.count_ops()['cx'])
    combinations, normalizations = RIIM_gate_combinations_sampled(n, num_cnots, num_max)
    new_circs = []
    for extra_set in combinations:
        tmp_circs = []
        for j in extra_set:
            qc = QuantumCircuit()
            qregs = transpile_result.qregs
            qubits = []
            for qreg in qregs:
                if not qc.has_register(qreg):
                    qc.add_register(qreg)
                qubits.extend(qreg)
            cregs = circuit.cregs
            clbits = []
            for creg in cregs:
                if not qc.has_register(creg):
                    qc.add_register(creg)
                clbits.extend(creg)
            cx_count = 0

            for op in ops:
                if op[0].name == 'cx':
                    if j[cx_count] == 0:
                        apply_extra_cnots(qc, *op[1], 0)
                        cx_count += 1
                    else:
                        apply_extra_cnots(qc, *op[1], j[cx_count])
                        cx_count += 1
                elif op[0].name == 'u1':
                    apply_u1(qc, *op[0].params, op[1][0])
                # elif op[0].name == 'h':
                #    apply_h(qc, *op[0].params, op[1][0])
                elif op[0].name == 'u2':
                    apply_u2(qc, *op[0].params, op[1][0])
                elif op[0].name == 'u3':
                    apply_u3(qc, *op[0].params, op[1][0])
                elif op[0].name == 'barrier':
                    qc.barrier()
                elif op[0].name == 'measure':
                    qc.measure(op[1][0], op[2][0])
                else:
                    print("ERROR:", op[0].name)
            tmp_circs.append(qc)
        new_circs.append(tmp_circs)
    return new_circs, RIIM_coeffs(n, num_cnots), normalizations

####################################################################
### Functions for performing FIIM and RIIM
####################################################################
def run_FIIM_extrapolation(n, circuit, observable, backend, shots, noise_model, coupling_map):
    """
    Performs FIIM with user-defined observable. 
    :param n: (int) order of FIIM to be performed 
    :param circuit: (QuantumCircuit) The quantum circuit to add CNOTs to.
    :param observable: (function) Function which takes in Qiskit experiment output dictionary and returns a float.
    :param backend: (Qiskit Backend) Backend to run circuits on. 
    :param shots: (int) Number of shots for each experiment.
    :param noise_model: (NoiseModel) Noise model used for simulation.
    :param coupling_map: (list) Coupling map pulled from Qiskit device.
    :return: (float) extrapolated value
    """
    FIIM_circs, r_vals = FIIM_generate_circs(n, circuit)
    FIIM_results = [execute(c,
                    backend,
                    shots=shots,
                    noise_model=noise_model,
                    coupling_map = coupling_map,
                    optimization_level=0).result().get_counts() for c in FIIM_circs]
    FIIM_obs = [observable(r, shots) for r in FIIM_results]
    fitted = np.poly1d(np.polyfit(r_vals[:n+1], FIIM_obs[:n+1], n))(0.0)
    return fitted

####################################################################
def run_RIIM_extrapolation(n, circuit, observable, backend, shots, noise_model, coupling_map):
    """
    Performs RIIM with user-defined observable. 
    :param n: (int) order of FIIM to be performed 
    :param circuit: (QuantumCircuit) The quantum circuit to add CNOTs to.
    :param observable: (function) Function which takes in Qiskit experiment output dictionary and returns a float.
    :param backend: (Qiskit Backend) Backend to run circuits on. 
    :param shots: (int) Number of shots for each experiment.
    :param noise_model: (NoiseModel) Noise model used for simulation.
    :param coupling_map: (list) Coupling map pulled from Qiskit device.
    :return: (float) extrapolated value
    """
    RIIM_circs, a_coeffs = RIIM_generate_circs(n, circuit)
    RIIM_obs = [[observable(execute(c,
                    backend,
                    shots=shots,
                    noise_model=noise_model,
                    coupling_map = coupling_map,
                    optimization_level=0).result().get_counts(), shots) for c in x] for x in RIIM_circs]
    RIIM_collapse = [sum(x) for x in RIIM_obs]
    fitted = sum([a*b for a,b in zip(a_coeffs, RIIM_collapse)])
    return fitted

def run_RIIM_extrapolation_sampled(n, circuit, num_max, normalize_shots, observable, backend, shots, noise_model, coupling_map):
    RIIM_circs, a_coeffs, normalizations = RIIM_generate_circs_sampled(n, circuit, num_max)
    """
    Performs RIIM with user-defined observable, with randomized sampling. 
    :param n: (int) order of FIIM to be performed 
    :param circuit: (QuantumCircuit) The quantum circuit to add CNOTs to.
    :param num_max: (int) maximum number of possible circuit permutations to be included for each group. 
    :param normalize shots: (bool) maximum number of possible circuit permutations to be included for each group. 
    :param observable: (function) Function which takes in Qiskit experiment output dictionary and returns a float.
    :param backend: (Qiskit Backend) Backend to run circuits on. 
    :param shots: (int) Number of shots for each experiment.
    :param noise_model: (NoiseModel) Noise model used for simulation.
    :param coupling_map: (list) Coupling map pulled from Qiskit device.
    :return: (float) extrapolated value
    """
    if not normalize_shots: 
        RIIM_obs = [[observable(execute(c,
                        backend,
                        shots=shots,
                        noise_model=noise_model,
                        coupling_map = coupling_map,
                        optimization_level=0).result().get_counts(), shots) for c in x] for x in RIIM_circs]
    else: 
        RIIM_obs = [[observable(execute(c,
                    backend,
                    shots=int(round(shots*y)),
                    noise_model=noise_model,
                    coupling_map = coupling_map,
                    optimization_level=0).result().get_counts(), int(round(shots*y))) for c in x] 
                    for x, y in zip(RIIM_circs, normalizations)]
    RIIM_collapse = [sum(x) for x in RIIM_obs]
    RIIM_collapse_normed = [a*b for a,b in zip(RIIM_collapse, normalizations)]
    fitted = sum([a*b for a,b in zip(a_coeffs, RIIM_collapse_normed)])
    return fitted