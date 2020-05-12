# RIIM Tools

This folder contains implementations for the Fixed Identity Insertion Method (FIIM) and Random Identity Insertion Method (RIIM) proposed in "Resource Efficient Zero Noise Extrapolation with Identity Insertions" (https://arxiv.org/abs/2003.04941). 

## RIIM_tools.py

Contains utility functions for performing RIIM on a quantum circuit in IBM's Qiskit. The main functions included to be used by the user are 

```python
run_FIIM_extrapolation(n, circuit, observable, backend, shots, noise_model, coupling_map)
run_RIIM_extrapolation(n, circuit, observable, backend, shots, noise_model, coupling_map)
run_RIIM_extrapolation_sampled(n, circuit, num_max, normalize_shots, 
observable, backend, shots, noise_model, coupling_map)
```

run_FIIM_extrapolation takes a quantum circuit and and observable function which converts the output of the quantum circuit to a float, along with other input specific to Qiskit's Aer simulator, and executes the FIIM protocol to estimate the extrapolated value of the observable. 

run_RIIM_extrapolation takes a quantum circuit and and observable function which converts the output of the quantum circuit to a float, along with other input specific to Qiskit's Aer simulator, and executes the RIIM protocol to estimate the extrapolated value of the observable. The number of shots in this function refers to the number of shots in each permutation represented in each component of O({e1,e2,e3...}). 

run_RIIM_extrapolation_sampled executes the RIIM protocol for estimating the observable, but with the modification that the number of permutations sampled for each O({e1,e2,e3...}) is bounded by a user-defined n_max. Importantly, there is a boolean flag normalize_shots which can be set to TRUE in order to retain the same total number of shots in an experiment as with normal RIIM, or set to FALSE where the number of shots in O({e1,e2,e3...}) is dependent on n_max. This is done so that both implementations of RIIM can be compared with the same statistical precision. 

## Usage

See the included Jupyter notebook applying_RIIM_tools.py for usage examples.

