# Zero Noise Extrapolation

In addition to readout errors, two-qubit gate noise is the main challenge for complex quantum algorithms on noisy intermediate-scale quantum (NISQ) computers. These errors are a significant challenge for making accurate calculations for quantum chemistry, nuclear physics, high energy physics, and other emerging scientific and industrial applications. There are two proposals for mitigating two-qubit gate errors: error-correcting codes and zero-noise extrapolation. This repo contains the code for the latter, studying it in detail and proposing modifications to existing approaches.  In particular, we propose a random identity insertion method (RIIM) that can achieve competitive asymptotic accuracy with far fewer gates than the traditional fixed identity insertion method (FIIM). For example, correcting the leading order depolarizing gate noise requires n<sub>CNOT</sub>+2 gates for RIIM instead of 3n<sub>CNOT</sub> gates for FIIM. This significant resource saving may enable more accurate results for state-of-the-art calculations on near term quantum hardware.

Details of the methods can be found in [2003.04941 [quant-ph]](https://arxiv.org/abs/2003.04941).

