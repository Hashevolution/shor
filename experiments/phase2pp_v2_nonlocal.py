"""
Phase 2'' v2: NON-LOCAL magic on the RT surface (the decisive test).

v1 result: a single-qubit phase gate on one bond injects boundary magic but does
NOT make the area operator state-dependent -- because single-qubit bond magic is
LOCAL magic, removable, not the gravity-relevant kind. Exactly the theorem's
subtlety (Cao et al.): you need NON-LOCAL magic.

Here the two perfect tensors are joined by TWO bonds (the RT surface = 2 bonds),
and we dress that surface with a 2-qubit gate G4:
  baseline  I (x) I   -- stabilizer net
  Clifford  CZ        -- adds cross-surface entanglement, NO magic
  magic     CS / cT   -- non-local (2-qubit) non-Clifford gate = non-local magic
Question: does Var over bulk states of S_A stay ~0 for baseline+CZ (trivial area
operator) but become >0 for the non-Clifford gates (state-dependent area op)?
"""

import numpy as np
from phase2pp_holographic_magic import (
    five_qubit_code, perfect_tensor, entropy_region, m2_full, bloch)


def network_state_2bond(T, phi1, phi2, G4):
    """two perfect tensors joined by 2 bonds dressed with 4x4 gate G4.
    boundary (6 qubits): [A: T1.p1,p2,p3] + [B: T2.p3,p4,p5]."""
    T1 = np.tensordot(T, phi1, axes=([5], [0]))   # p1..p5
    T2 = np.tensordot(T, phi2, axes=([5], [0]))   # p1..p5
    g = G4.reshape(2, 2, 2, 2)                     # [a,b,a',b']
    psi = np.einsum('pqrab,abcd,cdxyz->pqrxyz', T1, g, T2)
    psi = psi.reshape(-1)
    return psi / np.linalg.norm(psi)


def var_SA(T, G4, bulk_states, phi2, region=(0, 1, 2)):
    vals = [entropy_region(network_state_2bond(T, p1, phi2, G4), 6, list(region))
            for p1 in bulk_states]
    return float(np.var(vals)), float(np.mean(vals)), float(min(vals)), float(max(vals))


def main():
    zeroL, oneL = five_qubit_code()
    T = perfect_tensor(zeroL, oneL)
    phi2 = bloch(0, 0)
    bulk_states = [bloch(t, p) for t in np.linspace(0, np.pi, 4)
                   for p in np.linspace(0, 2 * np.pi, 4, endpoint=False)]

    eye4 = np.eye(4, dtype=complex)
    CZ = np.diag([1, 1, 1, -1]).astype(complex)          # Clifford, no magic
    CS = np.diag([1, 1, 1, 1j]).astype(complex)          # controlled-S, magic
    cT = np.diag([1, 1, 1, np.exp(1j * np.pi / 4)])      # controlled-T, magic

    gates = [("baseline  I x I ", eye4, "stabilizer"),
             ("Clifford  CZ     ", CZ, "Clifford (no magic)"),
             ("magic     CS     ", CS, "non-Clifford (magic)"),
             ("magic     cT     ", cT, "non-Clifford (magic)")]

    print("2-bond holographic network; RT surface dressed by 2-qubit gate G4")
    print("region A = T1.{p1,p2,p3};  bulk1 swept over Bloch sphere (behind A)\n")
    print(f"{'gate':>18} {'Var_bulk(SA)':>13} {'meanSA':>8} {'SA range':>16} "
          f"{'M2_bdy':>8} {'kind':>22}")
    for name, G, kind in gates:
        v, m, lo, hi = var_SA(T, G, bulk_states, phi2)
        mb = m2_full(network_state_2bond(T, bloch(0.7, 0.3), phi2, G), 6)
        print(f"{name:>18} {v:13.3e} {m:8.4f} [{lo:6.4f},{hi:6.4f}] "
              f"{mb:8.4f} {kind:>22}")

    print("\n--- decisive diagnostic ---")
    print("trivial area op  -> Var_bulk(S_A) ~ 0  (baseline, and ideally CZ)")
    print("state-dependent  -> Var_bulk(S_A) > 0  (CS, cT)")
    print("If CZ stays ~0 but CS/cT lift Var>0: NON-LOCAL MAGIC (not entanglement)")
    print("makes the area operator state-dependent = the gravity-relevant effect.")


if __name__ == "__main__":
    main()
