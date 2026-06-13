# Shor's Algorithm — Multi-Base Order Finding and Noise-as-Resource Mechanism

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20679807.svg)](https://doi.org/10.5281/zenodo.20679807)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)

A numpy-only implementation and analysis of Shor's quantum factoring algorithm with:

1. **Five theorems** on noise-invariant determinism, logarithmic coverage time, exact noise scaling, conditional Regev compatibility, and a hybrid (C) + Regev *b*-trick factoring algorithm.
2. **Universal trial-level boundary-flip mechanism** observation: 13 base sets × 12 σ values × 200 trials = 31,200 trial-measurements at `(N, d) = (437, 4)`.
3. **Engineered amplification** demonstration: removing the (C) augmentation amplifies the per-seed |SR| by ~5× — confirming the borderline-trial population is the bottleneck.
4. **Algorithm-structure regime map** for noise-as-resource susceptibility in multi-base quantum factoring.

> Main paper: **[paper.md](paper.md)** (Markdown) / **[paper.tex](paper.tex)** (LaTeX).
> A Korean theoretical companion is at [README.ko.md](README.ko.md).

## Headline results

### Theorems 1–5 (paper §3.1–3.5)

- **Theorem 1 — Noise-invariant determinism**: maintaining an accumulated exponent candidate `L = lcm(r_a₁, ..., r_aₖ)` and augmenting standard continued-fraction post-processing with a divisor search over `L` yields a procedure `(C)` that recovers the order `r_a` of any new base `a` deterministically whenever `r_a | L`, *independent of the measurement distribution*.
- **Theorem 2 — Logarithmic coverage time**: for a semiprime `N = pq`, `E[K_λ] ≤ 1 + Σ_{ℓ | λ(N)} 1/(ℓ^{s_ℓ} - 1)`, refining a corollary of Pomerance et al. 2017.
- **Theorem 3 — Exact noise scaling**: under destructive noise, `E[K_λ^{alg}(η)] = E[K_λ^{ideal}] / g_M(η)`.
- **Theorem 4 — Conditional Regev compatibility**: `(C)` applies coordinate-wise to Regev's multi-base measurements under a marginal-distribution assumption.
- **Theorem 5 — Hybrid (C) + Regev b-trick factoring**: empirically requires ~5× fewer runs than standalone Regev at `N ∈ {437, 1147, 2491, 4087}`, with **Lemma 5.1** giving a closed-form per-`b` nontrivial-sqrt probability bound.

Verified across **17,700 measurements at six composite sizes**, zero theorem violations.

### §3.6 — Multi-boundary mechanism observation (this release)

A high-statistics measurement at `(N, d) = (437, 4)` (13 seeds × 12 σ values × 200 trials = 31,200 trial-measurements) shows:

- **13/13 seeds exhibit a boundary-flip mechanism** — universal at the base-set level.
- **K-bin boundary distribution**: 76.9% K=1/K=2, 15.4% K=2/K=3, 7.7% K=3↔K=1 long-jump.
- **σ-curve direction asymmetry**: positive-direction seeds saturate + decline; negative-direction seeds monotonically worsen (a consequence of K-distribution skew).
- **Direction independence from K_baseline**: seeds with identical K_baseline = 1.720 can show opposite SR directions, confirming direction is determined by base-set-specific K-distribution structure, not by aggregate K_mean.
- **Net SR at σ = 0.050**: mean = +0.144%, t = 0.51, p = 0.31 — *not statistically significant*. Mechanism is universal; net direction is base-set-stochastic.

### Cross-cell verification at `(N, d) = (1147, 2)` (★ new)

A 5-seed σ-scan at K_baseline ≈ 2.92 reveals a *richer* mechanism than the K=1/K=2 boundary at (437, 4):

- **High-K rescue** (new pattern): trials at K = 8, 11, 15, 20 move to moderate K = 4, 5 under noise.
- **Boundary diversification**: only 40% of seeds show K=1/K=2 as the dominant transition (vs 77% at (437, 4)).
- Cross-seed mean SR = +3.35%, p ≈ 0.12 (t-distribution) — marginal, dominated by high-K-rescue seeds.

### Engineered amplification (★ new)

| Variant | K_baseline | Per-seed \|SR\| | Mechanism |
|---|---|---|---|
| Full hybrid (with (C)) | 2.08 | 0–1.16% | (C) buffers borderline trials into K=1 |
| Over-thinned (smallest convergent only) | 19.87 (9.55×) | **0.00%** | no borderline population left to flip |
| Mild-thinned (no (C) augmentation) | 2.92 (1.4×) | **4.03–4.44% (~5× amplified)** | (C) buffer removed; mechanism exposed |

The (C) augmentation acts as a *noise buffer*. Its removal amplifies per-seed |SR| ~5× but leaves direction stochastic.

### Algorithm-structure regime map (testable conjecture)

| Algorithm structure | SR magnitude | Source |
|---|---|---|
| Single-base Shor (1994) | small (≤1%) | **predicted** |
| Multi-base Regev (LLL post-processing) | negative | **predicted** |
| Hybrid full (this work) | +0.14% at (437,4) | **measured** |
| Hybrid mild-thinned | ~5× amplified | **measured** |
| Hybrid over-thinned | zero | **measured** |

The regime map predicts that noise-as-resource in quantum factoring is naturally maximized in a *multi-base + per-coordinate independent recovery + no buffer* variant — neither standalone Shor nor Regev. Verification at standalone Shor and Regev is left to follow-up work.

## What this paper does *not* claim

- **No cryptographic implication**: the effect is too small (1–7 trial flips per 200 trials, ±0.3–2% per seed) and there is no statistically significant net direction. RSA factoring estimates are unchanged.
- **No paradigm shift**: this is a mechanism observation paper, not a breakthrough algorithm. The contribution is *conceptual* — a clean bridge between quantum factoring and the noise-as-resource literature (Benzi et al. 1981 on classical SR, Plenio–Huelga 2008 on ENAQT).

## Repository layout

| File | Role |
|------|------|
| `paper.md` / `paper.tex` | **The formal paper (English): Theorems 1–5 + Lemma 5.1 + §3.6 mechanism observation** |
| `README.md` | This file (English) |
| `README.ko.md` | Korean theoretical companion (algorithm walkthrough, basic Shor mathematics) |
| `release_notes_v0.2.0.md` | Detailed release notes for v0.2.0 |
| `compile-notes.md` | LaTeX compilation + arXiv/Zenodo upload guide |
| `classical.py` | Classical reduction + classical order computation (baseline) |
| `shor.py` | Quantum order finding (numpy state-vector simulation) |
| `multi_base.py` | Multi-base λ(N) accumulation with `(C)` post-processing |
| `noise.py` | 6 noise models + 5 hybrid-b-trick corruption models |
| `demo.py` | Factoring demo + recovery rate comparison |
| `experiments/` | Experiments reproducing Theorems 2–5 + §3.6 |
| `experiments/sigma_scan_437.py` | §3.6 baseline σ-scan at (437, 4) (3 seeds × 12 σ × 200 trials) |
| `experiments/sigma_scan_437_extend.py` | §3.6 extended seeds 4–13 + K-histogram backfill |
| `experiments/sigma_scan_general.py` | σ-scan at arbitrary (N, d) cell |
| `experiments/analyze_histograms.py` | Per-seed K-bin flip identification |
| `experiments/sr_amplification.py` | Over-thinned amplification test (null result) |
| `experiments/sr_mild_amplification.py` | Mild-thinned amplification test (5× amplification) |
| `experiments/*_results.txt`, `*_histograms.txt` | Raw measurement data |
| `summary.md` / `frontier.md` / `hypotheses.md` / `roadmap.md` | Research notes in Korean (working diary) |

## Quick start

```powershell
pip install numpy

python demo.py                              # default factoring demo (15, 21, 35)
python demo.py 21                           # factor N=21
python demo.py --multi 91                   # multi-base mode
python demo.py --compare 33 35 77           # measurement count comparison
python classical.py 21                      # classical baseline

# Reproduce §3.6 results
python -m experiments.sigma_scan_437                # baseline σ-scan at (437, 4)
python -m experiments.sigma_scan_437_extend         # extended seeds + histograms
python -m experiments.analyze_histograms            # per-seed flip identification
python -m experiments.sigma_scan_general 1147 2 5 100 compact   # cross-cell verification
python -m experiments.sr_amplification              # over-thinned (null result)
python -m experiments.sr_mild_amplification         # mild thinned (5× amplification)
```

All experimental scripts are *resumable* (immediate per-cell save, skip-existing on re-run). Raw measurement data is committed under `experiments/*.txt`.

## Simulation limits

- Qubits ≈ 3 · ⌈log₂ N⌉. State-vector dimension 2^{t+n} → memory explodes as N grows.
- We use a two-stage trick for efficiency:
  1. Measure the work register first to collapse the counting register to a partial state.
  2. Apply numpy FFT directly to the counting register (inverse QFT = normalized DFT).
- This trick yields the same measurement statistics as the full quantum circuit.
- Largest N tested: 4087 = 61 · 67 (24-qubit counting register, 256 MB at complex128).

## Citation

```
Hashevolution. (2026). A Noise-Invariant Determinism Theorem for Multi-Base
Post-Processing in Shor's Order Finding (Version 0.2.0). Zenodo.
https://doi.org/10.5281/zenodo.20679807
```

BibTeX:
```bibtex
@software{hashevolution_shor_2026,
  author       = {Hashevolution},
  title        = {{A Noise-Invariant Determinism Theorem for
                   Multi-Base Post-Processing in Shor's Order
                   Finding}},
  month        = jun,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v0.2.0},
  doi          = {10.5281/zenodo.20679807},
  url          = {https://doi.org/10.5281/zenodo.20679807}
}
```

See `CITATION.cff` for the GitHub-rendered citation widget (BibTeX, APA, etc. auto-generated).

## License

[MIT](LICENSE)

## References

- Shor, P. W. (1994). *Algorithms for quantum computation*. FOCS.
- Nielsen, M. A. & Chuang, I. L. (2000). *Quantum Computation and Quantum Information*. Cambridge.
- Knill, E. (1995). *On Shor's quantum factor finding algorithm*. LANL tech report.
- McAnally, D. (2001). *A Refinement of Shor's Algorithm*. arXiv:quant-ph/0112055.
- Bach, E. & Shallit, J. (1996). *Algorithmic Number Theory*, Vol. 1. MIT Press.
- Pomerance, C. et al. (2017). *Expected number of generators of finite groups*. arXiv:1707.07193.
- Regev, O. (2023). *An efficient quantum factoring algorithm*. arXiv:2308.06572. JACM 2025.
- Ragavan, S. & Vaikuntanathan, V. (2023). *Space-Efficient and Noise-Robust Quantum Factoring*. arXiv:2310.00899.
- Ekerå, M. (2024). *On the success probability of quantum order finding*. arXiv:2201.07791.
- Benzi, R., Sutera, A., & Vulpiani, A. (1981). *The mechanism of stochastic resonance*. J. Phys. A.
- Plenio, M. B. & Huelga, S. F. (2008). *Dephasing-assisted transport: quantum networks and biomolecules*. New J. Phys.
- Wellens, T., Shatokhin, V., & Buchleitner, A. (2004). *Stochastic resonance*. Rep. Prog. Phys.

A full bibliography is in `paper.md` / `paper.tex` §References.

## Acknowledgements

This work benefited from extensive iteration. Initial single-seed observations gave way to multi-seed validations that overturned several over-claims (six explicit retractions are documented in paper §3.6). The published narrative reflects the corrections, not the original mistakes.
