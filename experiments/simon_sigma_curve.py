"""
simon_sigma_curve.py — Simon's algorithm σ-curve closed-form 검증.

목적:
  Simon's algorithm 의 phase noise σ-curve 가 우리 universal closed form
  p(σ) = p_∞ + (p_0 - p_∞) · exp(-α σ²) 에 적용되는지 검증.

Theoretical derivation:
  Simon's circuit 의 output measurement 분포 (s = hidden string):
    Noise-free: P_0(y) = 1/2^{n-1} if y · s = 0 (mod 2), else 0.
    Phase noise σ: E[|FFT(a · e^{iε})_y|²] = (1-e^{-σ²})/2^n + e^{-σ²} · P_0(y).

  Orthogonality probability:
    p_o(σ) = Σ_{y: y·s=0} P_σ(y)
           = 2^{n-1} · (1-e^{-σ²})/2^n + e^{-σ²} · 1
           = (1-e^{-σ²})/2 + e^{-σ²}
           = 1/2 + (1/2) · e^{-σ²}
           = p_∞ + (p_0 - p_∞) · exp(-σ²)
    with p_0 = 1, p_∞ = 1/2, α = 1.

K-distribution:
  Simon 의 success criterion 은 n-1 개 linearly-independent orthogonal y 수집.
  매 measurement 가 prob p_o(σ) 로 orthogonal, conditional on orthogonal 그 y 가
  prior 들과 independent 일 prob ≈ 1 - 2^{-(n-1-k_collected)} (k_collected = 현재 모은
  independent count). 작은 n 에서 근사 (n-1)/p_o(σ) 사용.

설계:
  - n = 5 (32-element state space)
  - 100 random hidden strings s
  - σ ∈ {0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500}
  - K trials per (s, σ): 200

Reproduction:
  python -u -m experiments.simon_sigma_curve
"""
from __future__ import annotations
import math
import statistics
import time
from pathlib import Path

import numpy as np


N_BITS = 5
DIM = 1 << N_BITS  # 2^n
SIGMAS = [0.000, 0.025, 0.050, 0.100, 0.150, 0.200, 0.300, 0.500]
N_HIDDEN_STRINGS = 30
N_K_TRIALS = 200
MAX_RUNS = 30  # cap on K (n-1=4 needs at least 4 measurements)

RESULTS_FILE = Path("experiments/simon_sigma_curve_results.txt")


def simon_state_amplitudes(s: int, x0: int, n_bits: int) -> np.ndarray:
    """
    Simon's circuit post-oracle, pre-final-Hadamard state amplitudes.
    Input register has equal superposition over {x0, x0 ⊕ s}.
    """
    dim = 1 << n_bits
    amps = np.zeros(dim, dtype=np.complex128)
    amps[x0] = 1.0 / math.sqrt(2)
    amps[x0 ^ s] = 1.0 / math.sqrt(2)
    return amps


def _make_hadamard(n_bits: int) -> np.ndarray:
    """Hadamard^n matrix: H_{y,x} = (-1)^{popcount(x AND y)} / sqrt(2^n)."""
    dim = 1 << n_bits
    H = np.empty((dim, dim), dtype=np.float64)
    for y in range(dim):
        for x in range(dim):
            H[y, x] = -1.0 if bin(x & y).count("1") % 2 else 1.0
    H /= math.sqrt(dim)
    return H


_HADAMARD = _make_hadamard(N_BITS)


def measure_simon(s: int, sigma: float, rng: np.random.Generator) -> int:
    """
    Single Simon's measurement: simulate post-oracle state + phase noise + final
    Hadamard^n + measurement.  Returns y ∈ {0, ..., 2^n - 1}.

    Hadamard transform preserves XOR structure: orthogonal y values have
    nonzero amplitude (1/sqrt(2^{n-1})), non-orthogonal y have zero amplitude.
    """
    x0 = int(rng.integers(0, DIM))
    amps = simon_state_amplitudes(s, x0, N_BITS)
    if sigma > 0:
        phases = rng.normal(0.0, sigma, size=DIM)
        amps = amps * np.exp(1j * phases)
    out = _HADAMARD @ amps
    probs = np.abs(out) ** 2
    total = probs.sum()
    if total > 0:
        probs = probs / total
    else:
        probs = np.ones(DIM) / DIM
    y = int(rng.choice(DIM, p=probs))
    return y


def is_orthogonal(y: int, s: int) -> bool:
    """y · s = 0 (mod 2)."""
    return bin(y & s).count("1") % 2 == 0


def gf2_rank(vectors: list[int], n_bits: int) -> int:
    """Compute GF(2) rank of bit vectors (Gauss elimination)."""
    rows = list(vectors)
    rank = 0
    for col in range(n_bits - 1, -1, -1):
        # Find row with leading 1 in this column
        pivot = -1
        for i in range(rank, len(rows)):
            if (rows[i] >> col) & 1:
                pivot = i
                break
        if pivot == -1:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        # Eliminate column in other rows
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> col) & 1):
                rows[i] ^= rows[rank]
        rank += 1
        if rank == n_bits:
            break
    return rank


def simon_one_trial(s: int, sigma: float, rng: np.random.Generator) -> int:
    """K-loop: collect measurements until n-1 independent orthogonal y's gathered."""
    collected = []  # nonzero orthogonal y's
    for K in range(1, MAX_RUNS + 1):
        y = measure_simon(s, sigma, rng)
        if y == 0 or not is_orthogonal(y, s):
            continue
        collected.append(y)
        rank = gf2_rank(collected, N_BITS)
        if rank >= N_BITS - 1:
            return K
    return MAX_RUNS


def measure_orthogonal_prob(s: int, sigma: float, n: int, rng: np.random.Generator) -> float:
    """Monte Carlo estimate of P(y · s = 0)."""
    count = 0
    for _ in range(n):
        y = measure_simon(s, sigma, rng)
        if is_orthogonal(y, s):
            count += 1
    return count / n


def predicted_K_simon(p_o: float, n_bits: int, max_runs: int = MAX_RUNS) -> float:
    """
    Approximate E[K] for Simon's: need n-1 independent orthogonal y's.
    Each measurement: prob p_o orthogonal. Conditional on orthogonal, prob of
    being independent given k collected ≈ 1 - 2^{(k - (n-1))} (for k < n-1).

    Simplified: E[K] ≈ Σ_{k=0}^{n-2} 1 / (p_o · (1 - 2^{k - (n-1)}))
    """
    if p_o <= 0:
        return float(max_runs)
    expected = 0.0
    for k in range(n_bits - 1):
        prob_indep = 1.0 - 2.0 ** (k - (n_bits - 1))
        if prob_indep <= 0:
            return float(max_runs)
        expected += 1.0 / (p_o * prob_indep)
    return min(expected, float(max_runs))


def closed_form_p_o(sigma: float) -> float:
    """p_o(σ) = 1/2 + (1/2) · exp(-σ²)."""
    return 0.5 + 0.5 * math.exp(-sigma * sigma)


def main():
    t0 = time.time()
    lines = []
    header = (
        f"# Simon σ-curve closed-form 검증\n"
        f"# Model: p_o(σ) = 1/2 + (1/2)·exp(-σ²)  (p_0=1, p_∞=1/2, α=1)\n"
        f"# K-distribution: Σ_{{k=0}}^{{n-2}} 1 / (p_o · (1 - 2^{{k-(n-1)}}))\n"
        f"# n = {N_BITS}, dim = {DIM}, σ ∈ {SIGMAS}\n"
        f"# {N_HIDDEN_STRINGS} random hidden strings × {N_K_TRIALS} trials\n"
        f"# max_runs = {MAX_RUNS}\n\n"
    )
    print(header)
    lines.append(header)

    # First: verify p_o(σ) closed form across multiple s
    print("## p_o(σ) closed form vs measured\n")
    lines.append("## p_o(σ) closed form vs measured\n")
    tbl_hdr = "| σ | p_o predicted | p_o measured (avg over s) | diff |\n|---:|---:|---:|---:|\n"
    print(tbl_hdr, end="")
    lines.append(tbl_hdr)

    rng_s = np.random.default_rng(20260614)
    s_list = []
    while len(s_list) < N_HIDDEN_STRINGS:
        cand = int(rng_s.integers(1, DIM))  # s != 0
        s_list.append(cand)

    rng_po = np.random.default_rng(123)
    for sigma in SIGMAS:
        pred = closed_form_p_o(sigma)
        meas_list = [measure_orthogonal_prob(s, sigma, 500, rng_po) for s in s_list[:10]]
        meas = statistics.mean(meas_list)
        row = f"| {sigma:.3f} | {pred:.4f} | {meas:.4f} | {meas - pred:+.4f} |\n"
        print(row, end="")
        lines.append(row)

    # K-distribution measurements
    print("\n## K(σ) closed form vs measured (Simon's)\n")
    lines.append("\n## K(σ) closed form vs measured (Simon's)\n")

    all_pred = []
    all_meas = []

    for sigma in SIGMAS:
        p_o = closed_form_p_o(sigma)
        K_pred = predicted_K_simon(p_o, N_BITS)
        rng_k = np.random.default_rng(int(sigma * 1e6) * 31 + 7)
        Ks = []
        t_s = time.time()
        for s in s_list:
            for _ in range(N_K_TRIALS):
                K = simon_one_trial(s, sigma, rng_k)
                Ks.append(K)
        K_mean = statistics.mean(Ks)
        K_sd = statistics.stdev(Ks)
        ts = time.time() - t_s
        row = (
            f"σ={sigma:.3f}: p_o_pred={p_o:.4f}  K_pred={K_pred:.3f}  "
            f"K_meas={K_mean:.3f}±{K_sd:.2f}  diff={K_mean - K_pred:+.3f}  ({ts:.0f}s)\n"
        )
        print(row, end="")
        lines.append(row)
        all_pred.append(K_pred)
        all_meas.append(K_mean)

    if len(all_pred) >= 2:
        mm = statistics.mean(all_meas)
        ss_res = sum((m - p) ** 2 for m, p in zip(all_meas, all_pred))
        ss_tot = sum((m - mm) ** 2 for m in all_meas)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        rmse = math.sqrt(ss_res / len(all_meas))
        agg = (
            f"\n# Aggregate fit: R² = {r2:+.4f}, RMSE = {rmse:.3f}, "
            f"n = {len(all_meas)} σ values\n"
        )
        print(agg, end="")
        lines.append(agg)

    elapsed = time.time() - t0
    footer = f"# Elapsed: {elapsed:.1f}s\n"
    print(footer, end="")
    lines.append(footer)

    RESULTS_FILE.write_text("".join(lines), encoding="utf-8")
    print(f"# Results saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
