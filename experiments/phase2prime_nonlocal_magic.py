"""
Phase 2' (work-order-gap-hunt.md): NON-LOCAL magic vs curvature.

The literature check (phase2_RESULTS.md) showed Phase 2 measured the WRONG
quantity: gravity needs NON-LOCAL magic (magic not removable by local unitaries;
Cao et al. arXiv:2306.14996), not total stabilizer Renyi entropy. Here we redo
the curvature scan with non-local magic.

Primary measure -- single-qubit-irreducible magic (a non-local / genuinely
multipartite magic proxy that vanishes for stabilizer states AND for any product
of single-qubit magic states):
    M2_nl(psi) = min over V = (x) V_i of  M2( V psi ),  V_i in SU(2) per site.
It strips gravity-irrelevant single-site magic, leaving magic that single-qubit
rotations cannot localise.

Secondary (optimisation-free cross-check) -- cut-crossing Pauli weight:
    w_X = sum of Xi_P over Paulis P that act non-trivially on BOTH halves,
    Xi_P = <psi|P|psi>^2 / 2^n. (coarse: contaminated by stabilizer entanglement,
    but trend-informative.)

Same inhomogeneous TFIM and curvature dial h as Phase 2.
"""

import numpy as np
from scipy.optimize import minimize
from phase2_magic_curvature import ground_state, half_chain_entropy

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
B_TENSOR = np.stack([I2, X, Y, Z], axis=0)  # B[a,r,c] = (P_a)[r,c]


def fast_pauli_coeffs(psi, n):
    """all 4^n Pauli expectations c_P = <psi|P|psi> as a (4,)*n real tensor."""
    psi = np.asarray(psi, dtype=complex).reshape(-1)
    M = np.outer(psi.conj(), psi)            # M[r,c] = psi*_r psi_c = rho_{c,r}
    T = M.reshape([2] * (2 * n))             # [r_0..r_{n-1}, c_0..c_{n-1}]
    for k in range(n):
        rem = n - k                          # r_k at axis 0, c_k at axis rem
        T = np.tensordot(T, B_TENSOR, axes=([0, rem], [1, 2]))
    return T.real                            # C[a_0..a_{n-1}] = c_P


def m2_from_coeffs(C, n):
    s4 = np.sum(C ** 4)
    return float(-np.log2(s4 / (2 ** n)))


def su2(a, b, g):
    """single-qubit unitary Rz(g) Ry(b) Rz(a)."""
    rz_a = np.array([[np.exp(-1j * a / 2), 0], [0, np.exp(1j * a / 2)]])
    ry_b = np.array([[np.cos(b / 2), -np.sin(b / 2)],
                     [np.sin(b / 2), np.cos(b / 2)]], dtype=complex)
    rz_g = np.array([[np.exp(-1j * g / 2), 0], [0, np.exp(1j * g / 2)]])
    return rz_g @ ry_b @ rz_a


def apply_product_unitary(psi, angles, n):
    t = psi.reshape([2] * n)
    for i in range(n):
        Vi = su2(*angles[3 * i:3 * i + 3])
        t = np.moveaxis(np.tensordot(Vi, t, axes=([1], [i])), 0, i)
    return t.reshape(-1)


def nonlocal_magic(psi, n, restarts=4, seed=0):
    """single-qubit-irreducible magic = min over product unitaries of M2."""
    rng = np.random.default_rng(seed)

    def cost(angles):
        return m2_from_coeffs(fast_pauli_coeffs(
            apply_product_unitary(psi, angles, n), n), n)

    best = cost(np.zeros(3 * n))             # identity start
    for _ in range(restarts):
        x0 = rng.uniform(0, 2 * np.pi, 3 * n)
        res = minimize(cost, x0, method="L-BFGS-B",
                       options={"maxiter": 60, "ftol": 1e-7})
        best = min(best, float(res.fun))
    return max(best, 0.0)


def crossing_weight(C, n):
    """optimisation-free: fraction of Pauli weight straddling the center cut."""
    A = list(range(n // 2))
    Bb = list(range(n // 2, n))
    idx = np.indices((4,) * n)               # (n, 4,4,...)
    nz = idx > 0
    maskA = nz[A].any(axis=0)
    maskB = nz[Bb].any(axis=0)
    Xi = C ** 2 / (2 ** n)
    return float(Xi[maskA & maskB].sum())


def self_check():
    # product of T magic states -> single-qubit-reducible -> M2_nl ~ 0
    t1 = np.array([1, np.exp(1j * np.pi / 4)]) / np.sqrt(2)
    prod = t1
    for _ in range(2):
        prod = np.kron(prod, t1)
    n = 3
    full = m2_from_coeffs(fast_pauli_coeffs(prod, n), n)
    nl = nonlocal_magic(prod, n, restarts=6, seed=1)
    assert full > 0.3, full
    assert nl < 1e-2, f"product magic should have ~0 non-local magic, got {nl}"
    # stabilizer |000> -> 0
    z = np.zeros(2 ** n, dtype=complex); z[0] = 1
    assert nonlocal_magic(z, n, restarts=2) < 1e-9
    print(f"[self-check] product-T: full M2={full:.3f}  non-local M2={nl:.4f}"
          f"  (local magic correctly removed)  OK")


def run(n=8, g=1.0, hs=None):
    if hs is None:
        hs = np.round(np.linspace(-0.30, 0.30, 7), 4)
    print(f"\ninhomogeneous TFIM  n={n}  g={g}   NON-LOCAL magic vs curvature\n")
    print(f"{'h':>7} {'S(half)':>9} {'M2_full':>9} {'M2_nonlocal':>12} {'w_cross':>9}")
    rows = []
    for h in hs:
        psi = ground_state(n, h, g)
        C = fast_pauli_coeffs(psi, n)
        S = half_chain_entropy(psi, n)
        m_full = m2_from_coeffs(C, n)
        m_nl = nonlocal_magic(psi, n, restarts=4)
        wX = crossing_weight(C, n)
        rows.append((h, S, m_full, m_nl, wX))
        print(f"{h:7.3f} {S:9.5f} {m_full:9.5f} {m_nl:12.6f} {wX:9.5f}")
    r = np.array(rows)
    h, S, mf, mnl, wX = r.T
    print("\n--- gate 2' diagnostics ---")
    print(f"corr(S, M2_nonlocal)     = {np.corrcoef(S, mnl)[0,1]:+.4f}")
    print(f"corr(M2_full, M2_nonloc) = {np.corrcoef(mf, mnl)[0,1]:+.4f}")
    print(f"S         peak at h = {h[np.argmax(S)]:+.3f}  "
          f"(monotonic {np.all(np.diff(S)>0) or np.all(np.diff(S)<0)})")
    print(f"M2_full   peak at h = {h[np.argmax(mf)]:+.3f}")
    print(f"M2_nonloc peak at h = {h[np.argmax(mnl)]:+.3f}  "
          f"(monotonic {np.all(np.diff(mnl)>0) or np.all(np.diff(mnl)<0)})")
    print(f"w_cross   peak at h = {h[np.argmax(wX)]:+.3f}")
    print(f"M2_nonlocal range   = {mnl.min():.4f} .. {mnl.max():.4f}")
    return r


if __name__ == "__main__":
    self_check()
    run(n=8, g=1.0)
