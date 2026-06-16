"""
Phase 2 (work-order-gap-hunt.md): magic axis on an inhomogeneous spin chain.

Question (gate 2): on a curved (inhomogeneous-coupling) spin chain, does the
*magic* (stabilizer Renyi entropy) respond to curvature with a profile that is
NOT reducible to the entanglement profile? The handover note's spin-chain work
(sections 9.8-9.15) only ever measured one axis -- entanglement. The 2026 field
consensus is that gravity/curvature needs entanglement AND magic. If magic shows
an independent curvature response here, that is the closest "unstepped" spot.

Model: open inhomogeneous transverse-field Ising chain
    H = - sum_i J_i Z_i Z_{i+1} - g sum_i X_i
Rainbow-like curvature dial h: central bonds strengthened (h>0, negative
curvature / AdS side) vs outer bonds strengthened (h<0, positive curvature /
de Sitter side). h scan therefore also probes the *sign of curvature* = the
Phase-1/A (de Sitter) question.

Measured per h:
  S  = half-chain von Neumann entanglement entropy
  M2 = stabilizer 2-Renyi entropy (magic), M2 = -log2( (1/2^n) sum_P <P>^4 )
"""

import itertools
import numpy as np
from numpy.linalg import eigh

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, X, Y, Z]


def kron_op(op, k, n):
    """single-site operator op at site k of n, as full 2^n matrix."""
    mats = [op if i == k else I2 for i in range(n)]
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def kron_pair(op, k, n):
    """two-site ZZ-style operator op acting on sites k,k+1."""
    mats = [I2] * n
    mats[k] = Z
    mats[k + 1] = Z
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def couplings(n, h):
    """rainbow profile: bond strength vs distance from chain center.
    h>0 -> central bonds strong (rainbow / negative curvature R=-h^2).
    h<0 -> outer bonds strong (inverse rainbow / positive-curvature analog)."""
    bonds = n - 1
    center = (bonds - 1) / 2.0
    d = np.array([abs(b - center) for b in range(bonds)])
    return np.exp(-h * d)


def ground_state(n, h, g):
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    J = couplings(n, h)
    for k in range(n - 1):
        H -= J[k] * kron_pair(None, k, n)
    for k in range(n):
        H -= g * kron_op(X, k, n)
    w, v = eigh(H)
    return v[:, 0]


def half_chain_entropy(psi, n):
    a = n // 2
    M = psi.reshape(2 ** a, 2 ** (n - a))
    s = np.linalg.svd(M, compute_uv=False)
    p = s ** 2
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log(p)))


def stabilizer_renyi_entropy(psi, n):
    """M2 = -log2( (1/2^n) sum_P <psi|P|psi>^4 ), sum over all 4^n Paulis."""
    psi_t = psi.reshape([2] * n)
    s4 = 0.0
    for combo in itertools.product(range(4), repeat=n):
        phi = psi_t
        for k in range(n):
            if combo[k] != 0:
                phi = np.moveaxis(
                    np.tensordot(PAULIS[combo[k]], phi, axes=([1], [k])), 0, k)
        c = np.vdot(psi_t, phi).real
        s4 += c ** 4
    val = s4 / (2 ** n)
    return float(-np.log2(val))


def self_check():
    """stabilizer state |0..0> must have M2 = 0."""
    for n in (2, 3, 4):
        psi = np.zeros(2 ** n, dtype=complex)
        psi[0] = 1.0
        m = stabilizer_renyi_entropy(psi, n)
        assert abs(m) < 1e-9, f"stabilizer M2 not 0 at n={n}: {m}"
    # a T-magic state on 1 qubit: should be > 0
    t = np.array([1.0, np.exp(1j * np.pi / 4)]) / np.sqrt(2)
    assert stabilizer_renyi_entropy(t, 1) > 1e-6
    print("[self-check] stabilizer states M2=0, T-state M2>0  OK")


def run(n=8, g=1.0, hs=None):
    if hs is None:
        hs = np.linspace(-0.30, 0.30, 7)
    print(f"\ninhomogeneous TFIM  n={n}  g={g}  (h<0: +curv/deSitter side, "
          f"h>0: -curv/AdS side)\n")
    print(f"{'h':>7} {'S(half)':>10} {'M2(magic)':>11}")
    rows = []
    for h in hs:
        psi = ground_state(n, h, g)
        S = half_chain_entropy(psi, n)
        M = stabilizer_renyi_entropy(psi, n)
        rows.append((h, S, M))
        print(f"{h:7.3f} {S:10.5f} {M:11.6f}")
    rows = np.array(rows)
    h, S, M = rows[:, 0], rows[:, 1], rows[:, 2]
    # gate-2 diagnostics
    cc = np.corrcoef(S, M)[0, 1]
    s_argpk = h[int(np.argmax(S))]
    m_argpk = h[int(np.argmax(M))]
    s_mono = np.all(np.diff(S) > 0) or np.all(np.diff(S) < 0)
    m_mono = np.all(np.diff(M) > 0) or np.all(np.diff(M) < 0)
    print("\n--- gate 2 diagnostics ---")
    print(f"corr(S, M)         = {cc:+.4f}")
    print(f"S peak at h        = {s_argpk:+.3f}   (monotonic: {s_mono})")
    print(f"M2 peak at h       = {m_argpk:+.3f}   (monotonic: {m_mono})")
    print(f"S range            = {S.min():.4f} .. {S.max():.4f}")
    print(f"M2 range           = {M.min():.4f} .. {M.max():.4f}")
    indep = (abs(cc) < 0.9) or (s_argpk != m_argpk) or (s_mono != m_mono)
    print(f"\nGATE 2 (magic responds to curvature differently than entanglement"
          f" -> candidate gap): {'PASS (b: residual)' if indep else 'fail (a)'}")
    return rows


if __name__ == "__main__":
    self_check()
    run(n=8, g=1.0)
