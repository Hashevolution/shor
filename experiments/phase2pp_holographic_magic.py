"""
Phase 2'' (work-order-gap-hunt.md): the RIGHT arena -- a holographic code.

Phase 2/2' showed 1D spin chains reduce magic to criticality/entanglement. The
literature (Cao et al. arXiv:2306.14996) locates the magic->gravity effect in
holographic CODES: stabilizer codes have a TRIVIAL (state-independent) area
operator; NON-LOCAL magic is required for a NON-TRIVIAL (state-dependent) one.

Minimal faithful test: two [[5,1,3]] perfect tensors (the canonical perfect
tensor) joined by one bond = the RT surface for boundary region A. We:
  (baseline) check that with a stabilizer network, the boundary entropy S_A is
             INDEPENDENT of the bulk state behind the surface (trivial area op),
  (test)     dress the RT-surface bond with a tunable phase gate and ask whether
             S_A becomes bulk-state-DEPENDENT, and whether that dependence peaks
             at NON-CLIFFORD (magic) angles and vanishes at Clifford angles.

A clean magic signature = Var over bulk states of S_A is ~0 at Clifford bond
angles (0, pi/2, pi) and >0 in between (max non-stabilizerness near pi/4).
"""

import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PD = {"I": I2, "X": X, "Y": Y, "Z": Z}


def pauli_string(s):
    out = np.array([[1]], dtype=complex)
    for ch in s:
        out = np.kron(out, PD[ch])
    return out


def five_qubit_code():
    """logical |0_L>,|1_L> of the [[5,1,3]] perfect code via stabilizer projector."""
    gens = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    ZL = pauli_string("ZZZZZ")
    XL = pauli_string("XXXXX")
    dim = 32
    P = np.eye(dim, dtype=complex)
    for g in gens:
        P = P @ ((np.eye(dim) + pauli_string(g)) / 2)
    P0 = P @ ((np.eye(dim) + ZL) / 2)
    # extract the (rank-1) code state |0_L>
    w, v = np.linalg.eigh(P0 @ P0.conj().T)
    zeroL = v[:, np.argmax(w)]
    zeroL = zeroL / np.linalg.norm(zeroL)
    oneL = XL @ zeroL
    oneL = oneL / np.linalg.norm(oneL)
    # verify
    for g in gens:
        assert np.allclose(pauli_string(g) @ zeroL, zeroL, atol=1e-8), g
    assert np.allclose(ZL @ zeroL, zeroL, atol=1e-8)
    assert np.allclose(ZL @ oneL, -oneL, atol=1e-8)
    assert abs(np.vdot(zeroL, oneL)) < 1e-8
    return zeroL, oneL


def perfect_tensor(zeroL, oneL):
    """6-leg perfect tensor T[p1,p2,p3,p4,p5, b] (5 physical, 1 bulk)."""
    V = np.stack([zeroL, oneL], axis=1)        # (32,2): col b -> |b_L>
    return V.reshape(2, 2, 2, 2, 2, 2)         # p1..p5, b


def m2_full(psi, n):
    """full stabilizer 2-Renyi entropy via fast Pauli transform."""
    B = np.stack([I2, X, Y, Z], axis=0)
    M = np.outer(psi.conj(), psi).reshape([2] * (2 * n))
    T = M
    for k in range(n):
        rem = n - k
        T = np.tensordot(T, B, axes=([0, rem], [1, 2]))
    s4 = np.sum(T.real ** 4)
    return float(-np.log2(s4 / 2 ** n))


def entropy_region(psi, n, region):
    region = list(region)
    rest = [i for i in range(n) if i not in region]
    t = psi.reshape([2] * n).transpose(region + rest).reshape(2 ** len(region), -1)
    s = np.linalg.svd(t, compute_uv=False)
    p = (s ** 2)
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log(p)))


def bloch(theta, phi):
    return np.array([np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)],
                    dtype=complex)


def network_state(T, phi1, phi2, Gbond):
    """contract two perfect tensors on one bond dressed by Gbond (2x2);
    bulk inputs phi1 (behind region A) and phi2. Returns 8-qubit boundary state.
    boundary order: [A: T1.p1..p4] + [B: T2.p2..p5]."""
    T1 = np.tensordot(T, phi1, axes=([5], [0]))   # p1..p5  (bulk1 contracted)
    T2 = np.tensordot(T, phi2, axes=([5], [0]))   # p1..p5  (bulk2 contracted)
    # bond: T1.p5 -- Gbond -- T2.p1
    psi = np.einsum('pqrse, ef, fwxyz -> pqrswxyz', T1, Gbond, T2)
    psi = psi.reshape(-1)
    return psi / np.linalg.norm(psi)


def var_SA_over_bulk(T, Gbond, bulk_states, region, phi2):
    vals = [entropy_region(network_state(T, p1, phi2, Gbond), 8, region)
            for p1 in bulk_states]
    return float(np.var(vals)), float(np.mean(vals)), float(min(vals)), float(max(vals))


def main():
    zeroL, oneL = five_qubit_code()
    assert m2_full(zeroL, 5) < 1e-6, "logical state must be stabilizer (M2=0)"
    # perfect-tensor check: any 2-qubit marginal of |0_L> is maximally mixed
    s2 = entropy_region(zeroL, 5, [0, 1])
    assert abs(s2 - 2 * np.log(2)) < 1e-6, s2
    print(f"[setup] [[5,1,3]] perfect code OK: M2(|0_L>)=0, "
          f"2-qubit marginal entropy = {s2:.4f} = 2 ln2 (perfect)")

    T = perfect_tensor(zeroL, oneL)
    A = [0, 1, 2, 3]                       # T1 exposed legs = region A
    phi2 = bloch(0, 0)                     # fix second bulk input |0>
    # bulk states behind A's RT surface, swept over Bloch sphere
    bulk_states = [bloch(t, p) for t in np.linspace(0, np.pi, 4)
                   for p in np.linspace(0, 2 * np.pi, 4, endpoint=False)]

    # baseline: stabilizer network (bond = identity)
    v0, m0, lo0, hi0 = var_SA_over_bulk(T, I2, bulk_states, A, phi2)
    print(f"\n[baseline  G=I, Clifford net]  S_A mean={m0:.4f}  "
          f"range[{lo0:.4f},{hi0:.4f}]  Var_bulk(S_A)={v0:.3e}")
    print("   -> trivial area operator if Var ~ 0 (S_A independent of bulk state)")

    # scan bond phase gate G(theta)=diag(1,e^{i theta}); Clifford at 0,pi/2,pi,3pi/2
    print("\n[bond-phase scan]  G=diag(1,e^{i th}); Clifford at th=0,pi/2,pi")
    print(f"{'theta/pi':>9} {'Var_bulk(SA)':>13} {'meanSA':>8} {'M2_bdy':>8} {'clifford?':>10}")
    thetas = np.linspace(0, 1.0, 9) * np.pi
    ref = network_state(T, bloch(0.7, 0.3), phi2, I2)
    for th in thetas:
        G = np.array([[1, 0], [0, np.exp(1j * th)]], dtype=complex)
        v, m, lo, hi = var_SA_over_bulk(T, G, bulk_states, A, phi2)
        mb = m2_full(network_state(T, bloch(0.7, 0.3), phi2, G), 8)
        is_cliff = any(abs((th - c) % (2 * np.pi)) < 1e-9 or
                       abs((th - c) % (2 * np.pi) - 2 * np.pi) < 1e-9
                       for c in (0, np.pi / 2, np.pi, 3 * np.pi / 2))
        print(f"{th/np.pi:9.3f} {v:13.3e} {m:8.4f} {mb:8.4f} {str(is_cliff):>10}")

    print("\n--- gate 2'' diagnostics ---")
    print("magic signature = Var_bulk(S_A) ~0 at Clifford angles (0, .5pi, 1.0pi)")
    print("                  and >0 in between (peak near 0.25pi). If so: bulk MAGIC")
    print("                  on the RT surface makes the area operator STATE-DEPENDENT.")


if __name__ == "__main__":
    main()
