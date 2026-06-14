"""
shor_n_scaling.py — Shor closed-form (p_0, ρ) 의 N-scaling 측정.

목적:
  Cryptographic regime 거동 추정. 본 framework 의 SR effect bound 가
    |Δp_max| = |p_0 - ρ|
  이므로 N 이 커지면서 어떻게 변하는지 측정.

설계:
  - N ∈ {437, 1147, 2491, 4087, 8009}  (paper §3.6 N 들 + 큰 N)
  - QPE (no b-trick) success criterion 사용 (Shor pure 보다 robust + faster).
  - 각 N 에서 3-5 random (a, r_a) setups, p_0 와 ρ MC 측정.
  - K(σ=0.1) prediction 만 추가 (single σ for time budget).

Reproduction:
  python -u -m experiments.shor_n_scaling
"""
from __future__ import annotations
import math
import statistics
import time
from pathlib import Path

import numpy as np

from classical import classical_order
from multi_base import convergent_denominators, minimize_order
from shor import simulate_period_finding


N_VALUES = [437, 1147, 2491, 4087, 8009]
N_SETUPS_PER_N = 3
N_P0_SAMPLES = 1000
N_RHO_SAMPLES = 1000
MAX_RUNS = 20

RESULTS_FILE = Path("experiments/shor_n_scaling_results.txt")


def qpe_success(k: int, a: int, N_: int, Q: int, r_a: int) -> bool:
    cands = set(convergent_denominators(k, Q, N_ - 1))
    valid = [d_ for d_ in cands if d_ > 0 and pow(a, d_, N_) == 1]
    if not valid:
        return False
    r = minimize_order(a, N_, min(valid))
    return r > 0 and r == r_a


def find_setups(N_: int, count: int) -> list[tuple[int, int]]:
    """Pick (a, r_a) pairs with diverse r."""
    found = {}
    for a in range(2, N_):
        if math.gcd(a, N_) != 1:
            continue
        r = classical_order(a, N_)
        if r < 3:
            continue
        if r not in found:
            found[r] = a
        if len(found) >= 30:
            break
    rs = sorted(found.keys())
    if not rs:
        return []
    step = max(1, len(rs) // count)
    chosen = []
    for i in range(0, len(rs), step):
        chosen.append((found[rs[i]], rs[i]))
        if len(chosen) >= count:
            break
    return chosen[:count]


def measure_p0(a, N_, r_a, n, rng):
    success = 0
    for _ in range(n):
        m = simulate_period_finding(a, N_, rng=rng)
        if qpe_success(m.k, a, N_, m.Q, r_a):
            success += 1
    return success / n


def measure_rho(a, N_, r_a, Q, n, rng):
    success = 0
    for _ in range(n):
        k = int(rng.integers(0, Q))
        if qpe_success(k, a, N_, Q, r_a):
            success += 1
    return success / n


def main():
    t0 = time.time()
    lines = []
    header = (
        f"# Shor cryptographic N-scaling (closed-form parameters)\n"
        f"# Model: p(σ) = ρ + (p_0 - ρ) · exp(-σ²)\n"
        f"# Effect bound: |Δp_max| = |p_0 - ρ|\n"
        f"# N values: {N_VALUES}\n"
        f"# {N_SETUPS_PER_N} setups per N, MC samples: p_0={N_P0_SAMPLES}, ρ={N_RHO_SAMPLES}\n\n"
    )
    print(header)
    lines.append(header)

    rows = []  # (N, a, r_a, p_0, rho, delta)

    for N_ in N_VALUES:
        Q = 1 << (2 * max(1, (N_ - 1).bit_length()))
        setups = find_setups(N_, N_SETUPS_PER_N)
        section = f"## N = {N_} (Q = {Q}, {len(setups)} setups)\n"
        print(section)
        lines.append(section)

        for idx, (a, r_a) in enumerate(setups, 1):
            rng_p0 = np.random.default_rng(N_ * 991 + idx * 17)
            rng_rho = np.random.default_rng(N_ * 991 + idx * 23)
            t_mc = time.time()
            p0 = measure_p0(a, N_, r_a, N_P0_SAMPLES, rng_p0)
            rho = measure_rho(a, N_, r_a, Q, N_RHO_SAMPLES, rng_rho)
            mct = time.time() - t_mc
            delta = abs(p0 - rho)
            row = (
                f"  setup {idx}: a={a}, r_a={r_a}, p_0={p0:.4f}, ρ={rho:.4f}, "
                f"|Δ|={delta:+.4f}  ({mct:.0f}s)\n"
            )
            print(row, end="")
            lines.append(row)
            rows.append((N_, a, r_a, p0, rho, delta))

    # Summary table
    summary = "\n## Summary: |Δ| = |p_0 - ρ| scaling with N\n"
    print(summary)
    lines.append(summary)
    tbl_hdr = "| N | mean p_0 | mean ρ | mean \\|Δ\\| | max \\|Δ\\| | min \\|Δ\\| |\n"
    tbl_hdr += "|---:|---:|---:|---:|---:|---:|\n"
    print(tbl_hdr, end="")
    lines.append(tbl_hdr)
    for N_ in N_VALUES:
        rs = [r for r in rows if r[0] == N_]
        p0s = [r[3] for r in rs]
        rhos = [r[4] for r in rs]
        deltas = [r[5] for r in rs]
        if not p0s:
            continue
        row = (
            f"| {N_} | {statistics.mean(p0s):.4f} | {statistics.mean(rhos):.4f} | "
            f"{statistics.mean(deltas):.4f} | {max(deltas):.4f} | {min(deltas):.4f} |\n"
        )
        print(row, end="")
        lines.append(row)

    elapsed = time.time() - t0
    footer = f"\n# Elapsed: {elapsed:.1f}s ({elapsed/60:.1f} min)\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
