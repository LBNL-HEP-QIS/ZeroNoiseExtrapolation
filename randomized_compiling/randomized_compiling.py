"""randomize_compiling.py
"""
import qiskit
import random


def rc(basis_gates=None, optimization_level=1):
  """Convert randomized compiling circuits to gates.

  Returns:
    (list) Gates for randomized compiling.
  """
  gates = []
  for c, i in zip(_rc_circs(), range(len(_rc_circs()))):
    g = qiskit.compiler.transpile(
        c, basis_gates=basis_gates, optimization_level=optimization_level
      ).to_gate()
    g.name = 'RC%d' % i
    gates.append(g)
  return gates


def apply_rc(qc):
  """Apply randomized compiling to a given quantum circuit.

  For each `cx` gate in the given QuantumCircuit, four randomized compiling
  gates are chosen at random and applied to the control and target qubits
  such that,

    q0 ---- ... ---- RC1 ---- CTL ---- RC3 ---- ...
    q1 ---- ... ---- RC2 ---- TGT ---- RC4 ---- ...
    ...
    qn

  Args:
    qc (QuantumCircuit): The circuit to apply randomized compiling to.

  Returns:
    (QuantumCircuit) A circuit with randomized compiling.
  """
  # The data to work with from the input circuit
  data = qc.decompose().data
  qregs = qc.qregs
  cregs = qc.cregs

  # The circuit and its qubits to be returned
  qcrc = qiskit.QuantumCircuit()
  # Copy the input circuit q/cregs to the new circuit
  for qreg in qregs:
    if not qcrc.has_register(qreg):
      qcrc.add_register(qreg)
  for creg in cregs:
    if not qcrc.has_register(creg):
      qcrc.add_register(creg)

  # Loop over the data operators and add to the new circuit
  for op in data:
    if op[0].name == 'cx':
      _apply_rc(qcrc, op[1])
    else:
      pass
  return qcrc


def _rc_circs():
  """Randomized compiling circuits.

  Returns:
    (list) Circuits for randomized compiling.
  """
  # pylint: disable=multiple-statements
  c0 = qiskit.QuantumCircuit(1)
  c0.id(0); c0.id(0); c0.id(0); c0.id(0)

  c1 = qiskit.QuantumCircuit(1)
  c1.id(0); c1.x(0); c1.id(0); c1.x(0)

  c2 = qiskit.QuantumCircuit(1)
  c2.id(0); c2.y(0); c2.z(0); c2.y(0)

  c3 = qiskit.QuantumCircuit(1)
  c3.id(0); c3.z(0); c3.z(0); c3.z(0)

  c4 = qiskit.QuantumCircuit(1)
  c4.x(0); c4.id(0); c4.x(0); c4.x(0)

  c5 = qiskit.QuantumCircuit(1)
  c5.x(0); c5.x(0); c5.x(0); c5.id(0)

  c6 = qiskit.QuantumCircuit(1)
  c6.x(0); c6.y(0); c6.y(0); c6.z(0)

  c7 = qiskit.QuantumCircuit(1)
  c7.x(0); c7.z(0); c7.y(0); c7.y(0)

  c8 = qiskit.QuantumCircuit(1)
  c8.y(0); c8.id(0); c8.y(0); c8.x(0)

  c9 = qiskit.QuantumCircuit(1)
  c9.y(0); c9.x(0); c9.y(0); c9.id(0)

  c10 = qiskit.QuantumCircuit(1)
  c10.y(0); c10.y(0); c10.x(0); c10.z(0)

  c11 = qiskit.QuantumCircuit(1)
  c11.y(0); c11.z(0); c11.x(0); c11.y(0)

  c12 = qiskit.QuantumCircuit(1)
  c12.z(0); c12.id(0); c12.z(0); c12.id(0)

  c13 = qiskit.QuantumCircuit(1)
  c13.z(0); c13.x(0); c13.z(0); c13.x(0)

  c14 = qiskit.QuantumCircuit(1)
  c14.z(0); c14.y(0); c14.id(0); c14.y(0)

  c15 = qiskit.QuantumCircuit(1)
  c15.z(0); c15.z(0); c15.id(0); c15.z(0)

  return [c0, c1, c2, c3, c4, c5, c6, c7, c8,
          c9, c10, c11, c12, c13, c14, c15]


def _rc_random4(basis_gates=None, optimization_level=1):
  """Randomly select four RC gates.

  Returns:
    (list) Four random gates from the set of 16 possible RC gates.
  """
  rc_gates = rc(basis_gates, optimization_level)
  random4 = []
  while len(random4) < 4:
    r = random.randint(0, len(rc_gates) - 1)
    if r not in random4:
      random4.append(r)
  return [rc_gates[g] for g in random4]


def _apply_rc(qcrc, reg):
  """Apply randomized compiling gates to a single cx gate in-situ.

  Args:
    qcrc (QuantumCircuit): Circuit to add RC to.
    reg ([quantumregister.Qubit]): Register with target and control qubits.
  Returns:
    (QuantumCircuit) A circuit with RC applied.
  """
  random_rc = _rc_random4()
  ctl = reg[0].index
  tgt = reg[1].index
  qcrc.append(random_rc[0], [ctl])
  qcrc.append(random_rc[1], [tgt])
  qcrc.cx(ctl, tgt)
  qcrc.append(random_rc[2], [ctl])
  qcrc.append(random_rc[3], [tgt])
