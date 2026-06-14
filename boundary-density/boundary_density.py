"""
boundary_density.py — K-distribution geometry analyzer (boundary density rho_b).

Context
-------
The §3.6 noise-as-resource (SR) program established that the *net* noise effect is
null across every algorithm structure (|mean SR| <= 1.5%, no significance, direction
base-set-stochastic), while the *per-seed magnitude* distribution is algorithm-
dependent. This script formalizes the structural quantity behind that magnitude:
the **boundary density** rho_b, the fraction of trials sitting on a flippable K-bin
boundary.

It also makes precise an observation already implicit in the paper ("Pure Shor shows
the smallest magnitudes ... dead-trial dominance"): the mean K (`K_baseline`) is a
*misleading* proxy for boundary density. Decompose it:

    mean K  =  (dead / ceiling mass at >= MAX_RUNS)
             + (floor mass at K = 1)
             + (live boundary structure, 2 <= K < MAX_RUNS)

Only the last term — rho_b — is the population noise can actually flip. A large mean K
driven by dead mass (e.g. Shor's ~48% b-trick failures) carries *no* boundary density.

This is pure stdlib (no numpy/sympy): it reads the committed K-histograms in
../experiments/*_histograms.txt. Regenerate those via the engine as documented in
paper.md §3.6 (`python -m experiments.pure_shor_sr`, etc.).

Histogram file format (tab-separated, '#' comments):
    seed    sigma   K   count

Usage:
    python boundary-density/boundary_density.py                 # scan ../experiments
    python boundary-density/boundary_density.py FILE [FILE ...] # specific files
"""

from __future__ import annotations

import collections
import glob
import os
import sys
from dataclasses import dataclass

# A trial is "dead" (ceiling / max_runs failure) when K >= MAX_RUNS.
# NOTE: all current histograms use max_runs = 20. When the algorithm-generalization
# experiments land with a different cap, make this per-file (auto-detect the top bin).
MAX_RUNS = 20
FLOOR_K = 1  # a single successful run


@dataclass
class CellStats:
    """Boundary-density decomposition of one K-distribution."""

    n: int
    mean_K: float          # full mean (the headline "K_baseline")
    mean_K_live: float     # mean over live (K < MAX_RUNS) trials only
    f_dead: float          # P(K >= MAX_RUNS)      -- unflippable ceiling mass
    f_floor: float         # P(K == 1)             -- floor reservoir
    rho_b: float           # P(2 <= K < MAX_RUNS)  -- *boundary density*

    def row(self, label: str) -> str:
        return (
            f"{label:<30} n={self.n:>5}  "
            f"meanK={self.mean_K:7.3f}  meanK_live={self.mean_K_live:6.3f}  "
            f"dead={self.f_dead*100:5.1f}%  floor={self.f_floor*100:5.1f}%  "
            f"rho_b={self.rho_b*100:5.1f}%"
        )


def parse_histograms(path: str):
    """Return {(seed, sigma): {K: count}} from a histogram file."""
    cells: dict = collections.defaultdict(lambda: collections.defaultdict(int))
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) != 4 or parts[0] == "seed":
                continue
            seed, sigma, K, count = parts
            cells[(int(seed), float(sigma))][int(K)] += int(count)
    return cells


def decompose(hist: dict) -> CellStats:
    """Compute the boundary-density decomposition of one K-histogram {K: count}."""
    n = sum(hist.values())
    if n == 0:
        return CellStats(0, 0.0, 0.0, 0.0, 0.0, 0.0)
    total_K = sum(K * c for K, c in hist.items())
    dead = sum(c for K, c in hist.items() if K >= MAX_RUNS)
    floor = sum(c for K, c in hist.items() if K == FLOOR_K)
    boundary = sum(c for K, c in hist.items() if FLOOR_K < K < MAX_RUNS)
    live = n - dead
    live_K = sum(K * c for K, c in hist.items() if K < MAX_RUNS)
    return CellStats(
        n=n,
        mean_K=total_K / n,
        mean_K_live=(live_K / live) if live else 0.0,
        f_dead=dead / n,
        f_floor=floor / n,
        rho_b=boundary / n,
    )


def merge(cells: dict, sigma: float | None = None) -> dict:
    """Sum histograms across seeds (optionally filtered to one sigma)."""
    merged: dict = collections.defaultdict(int)
    for (seed, sig), hist in cells.items():
        if sigma is not None and abs(sig - sigma) > 1e-12:
            continue
        for K, c in hist.items():
            merged[K] += c
    return merged


def analyze_file(path: str) -> CellStats:
    """Print per-seed and aggregate decomposition for one file; return aggregate."""
    name = os.path.basename(path)
    cells = parse_histograms(path)
    sigmas = sorted({sig for (_seed, sig) in cells})
    base_sigma = sigmas[0] if sigmas else 0.0  # noise-free baseline

    print(f"\n=== {name}  (baseline sigma={base_sigma}) ===")
    for seed in sorted({seed for (seed, _sig) in cells}):
        print("  " + decompose(cells[(seed, base_sigma)]).row(f"seed {seed}"))

    agg = decompose(merge(cells, base_sigma))
    print("  " + "-" * 96)
    print("  " + agg.row("AGGREGATE (sigma=baseline)"))
    return agg


def find_data_files() -> list[str]:
    """Locate *_histograms.txt: prefer ../experiments, then ./experiments, ./data."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, "..", "experiments"),
        os.path.join(here, "experiments"),
        os.path.join(here, "data"),
        os.path.join(os.getcwd(), "experiments"),
    ]
    for d in candidates:
        hits = sorted(glob.glob(os.path.join(d, "*_histograms.txt")))
        if hits:
            return hits
    return []


def main(argv: list[str]) -> int:
    paths = argv if argv else find_data_files()
    if not paths:
        print("no *_histograms.txt found (looked in ../experiments)", file=sys.stderr)
        return 1

    print("boundary-density decomposition  (MAX_RUNS=%d, FLOOR_K=%d)" % (MAX_RUNS, FLOOR_K))
    print("rho_b = P(2 <= K < MAX_RUNS) is the flippable boundary population.")

    aggs = {}
    for path in paths:
        aggs[os.path.basename(path)] = analyze_file(path)

    print("\n=== cross-structure summary (noise-free baseline) ===")
    print(f"{'histogram file':<40}{'meanK':>9}{'dead%':>8}{'floor%':>8}{'rho_b%':>9}")
    for name, st in aggs.items():
        short = name.replace("_histograms.txt", "")
        print(f"{short:<40}{st.mean_K:>9.3f}{st.f_dead*100:>8.1f}"
              f"{st.f_floor*100:>8.1f}{st.rho_b*100:>9.1f}")
    print(
        "\nReading: a large meanK with small rho_b (e.g. pure_shor) means the high\n"
        "K_baseline is a dead-mass artifact, not boundary structure. rho_b is the real\n"
        "driver of SR magnitude — and it is set by algorithm structure, not by noise."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
