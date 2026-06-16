"""
Phase 2b: self-falsification of the Phase-2 magic-vs-curvature signal.
(work-order-gap-hunt.md, section 9.15 discipline: change the measurement, see
if the pattern survives. The note died 3x this way -- do the same here.)

Phase 2 found: magic M2 has an INTERIOR peak in curvature h while entanglement S
is monotonic. Before promoting to a "gap" we must rule out the boring reading:
"magic just tracks criticality" -- a homogeneous TFIM at g=1 is critical and is
known to have near-maximal magic, so a peak near h=0 could be nothing new.

Tests:
  T1 finer h grid: is the interior magic peak robust (not a coarse-grid artifact
     like the note's 9.11 two-point artifact)?
  T2 g-scan: does the magic peak in h TRACK the critical point (move with g)?
     If it always sits where local coupling = critical -> criticality re-discovery (a).
     If it stays put / moves differently -> curvature-specific residual (b).
  T3 finite size: peak location consistent across n (n=6 vs n=8)?
"""

import numpy as np
from phase2_magic_curvature import (
    ground_state, half_chain_entropy, stabilizer_renyi_entropy)


def scan(n, g, hs):
    S, M = [], []
    for h in hs:
        psi = ground_state(n, h, g)
        S.append(half_chain_entropy(psi, n))
        M.append(stabilizer_renyi_entropy(psi, n))
    return np.array(S), np.array(M)


def peak_h(hs, vals):
    # quadratic refine around argmax for a sub-grid peak estimate
    i = int(np.argmax(vals))
    if 0 < i < len(vals) - 1:
        x0, x1, x2 = hs[i - 1], hs[i], hs[i + 1]
        y0, y1, y2 = vals[i - 1], vals[i], vals[i + 1]
        denom = (y0 - 2 * y1 + y2)
        if abs(denom) > 1e-12:
            return x1 - 0.5 * (x2 - x0) * (y2 - y0) / (2 * denom) * 0.5 + 0.0 \
                if False else x1 - 0.5 * (y2 - y0) / denom * (x1 - x0)
    return hs[i]


def main():
    # T1 ---- finer h grid at n=6 (fast) ----
    print("T1  finer h grid (n=6, g=1.0)")
    hs = np.round(np.linspace(-0.30, 0.30, 13), 4)
    S, M = scan(6, 1.0, hs)
    print(f"{'h':>7} {'S':>9} {'M2':>10}")
    for h, s, m in zip(hs, S, M):
        print(f"{h:7.3f} {s:9.5f} {m:10.5f}")
    print(f"  S monotonic: {np.all(np.diff(S)>0) or np.all(np.diff(S)<0)}")
    print(f"  M2 peak (refined) h = {peak_h(hs, M):+.4f}  "
          f"M2 monotonic: {np.all(np.diff(M)>0) or np.all(np.diff(M)<0)}")

    # T2 ---- does the magic peak track criticality? scan g ----
    print("\nT2  magic-peak-in-h vs transverse field g (n=6)")
    print("    if peak moves with g the way criticality does -> (a) re-discovery")
    print(f"{'g':>6} {'h*(M2 peak)':>12} {'h*(S behaviour)':>18}")
    for g in (0.6, 0.8, 1.0, 1.2, 1.5):
        S, M = scan(6, g, hs)
        s_mono = np.all(np.diff(S) > 0) or np.all(np.diff(S) < 0)
        print(f"{g:6.2f} {peak_h(hs, M):12.4f} "
              f"{'monotonic' if s_mono else 'NON-monotonic':>18}")

    # T3 ---- finite size: n=6 vs n=8 peak location (coarser grid for n=8) ----
    print("\nT3  finite-size check of magic peak (g=1.0)")
    hsc = np.round(np.linspace(-0.30, 0.30, 7), 4)
    for n in (6, 8):
        S, M = scan(n, 1.0, hsc)
        print(f"  n={n}:  M2 peak h = {peak_h(hsc, M):+.4f}   "
              f"S monotonic = {np.all(np.diff(S)>0) or np.all(np.diff(S)<0)}")

    print("\n--- reading guide ---")
    print("(b) residual/gap if: M2 interior peak is robust (T1), and does NOT")
    print("    simply track the critical coupling as g varies (T2), and is")
    print("    finite-size stable (T3).")
    print("(a) re-discovery if: M2 peak sits at/moves to wherever the chain is")
    print("    locally critical -> 'critical TFIM has max magic', already known.")


if __name__ == "__main__":
    main()
