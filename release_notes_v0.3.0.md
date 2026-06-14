# Release notes — v0.3.0 (2026-06-14)

## Summary

**v0.3.0** is an *errata + extension* release of paper v0.2.1 (DOI
10.5281/zenodo.20681847).

The §3.6 boundary-flip mechanism described in v0.2.1 is replaced by a single
analytic closed form derived from standard phase-noise dephasing of the FFT.
We verify the closed form across **five quantum algorithm classes** (R² ∈
[0.88, 0.99]) and position the work as an *analytical complement* to
Yang-Markidis (arXiv:2605.16074, ICS Workshops '26).

## What changed

### Self-correction (paper §3.6.bis, new)

The "boundary flip + universal direction stochasticity + plateau/overload"
framework of v0.2.1 §3.6 is the **finite-trial expression of a single closed
form**:

```
p(σ) = ρ + (p_0 - ρ) · exp(-σ²)
E[K(σ)] = (1 - (1-p)^M) / p
```

derived directly from `E[|FFT(a·e^{iε})_k|²] = (1-e^{-σ²})/Q + e^{-σ²}·P_0(k)`.

### Cross-algorithm verification (new)

| Algorithm | R² | n |
|---|---:|---:|
| Grover | +0.88 | 40 |
| Shor pure (b-trick) | +0.95 | 40 |
| QPE isolated (no b-trick) | +0.96 | 40 |
| Simon | +0.99 | 8 |
| Hybrid (C)+b-trick (paper §3.6 setup) | +0.91 | 40 |

The same closed form fits all five — a small *organizing* contribution.

### Retracted (within §3.6)

1. *"Boundary flip" lexicon* as a distinct mechanism — it is the K-binning of a
   smooth p(σ) shift under finite trials.
2. *"Deterministic flip set within plateau"* reading — actually statistical.
3. *"Universal direction stochasticity"* as unexplainable — `sign(p_0 - ρ)`
   explains it per (a, b) setup.

### Retained (unchanged)

- All raw measurement data (31,200 K-measurements).
- The five-cell regime map predictions (5/5 measured).
- The conclusion that SR-based factoring acceleration is precluded.
- Theorems 1–5 of v0.2.1 (independent of §3.6 framework).

## New artifacts

- `sr_sigma_curve_model.md` — unified closed-form framework (§1-§7).
- `sr_generalization.md` — SR generalization scoping + Grover/QPE results.
- `experiments/grover_sr.py`, `grover_sr_focused.py`, `grover_sr_ceiling.py`,
  `grover_sigma_curve_model.py` — Grover suite.
- `experiments/shor_sigma_curve_model.py` — Shor closed-form fit.
- `experiments/qpe_isolated_sigma.py` — QPE internal consistency.
- `experiments/simon_sigma_curve.py` — Simon multi-coordinate.
- `experiments/hybrid_sigma_curve.py` — paper §3.6 setup direct fit.
- `experiments/shor_n_scaling.py` — cryptographic regime |Δ| scaling.
- `paper.md` §3.6.bis — self-correction section.

## Positioning

- **Analytical complement** to Yang-Markidis (arXiv:2605.16074, May 2026):
  their empirical noise propagation model `(1-ε)·P_s + ε·distractors` has the
  same structural form as our closed form; we provide `ε = 1 - exp(-σ²)`
  analytically.
- **Self-correction** of v0.2.1 §3.6 SR framework — boundary-flip lexicon
  retracted, closed-form replacement.
- **No new mechanism claim**. The exp(-σ²) form is the standard dephasing
  result; our contribution is its specific application to Shor-class SR claims
  and its cross-algorithm universality demonstration.

## What this is NOT

- Not a new quantum advantage.
- Not a new mechanism.
- Not a refutation of v0.2.1 — the data is retained; only the §3.6 *interpretation*
  is corrected.

## What this IS

- An honest scientific cycle: claim → verification → self-correction.
- A small bridge between quantum-algorithm SR claims and standard dephasing
  literature.
- A cross-algorithm organizing result.

## Related work

- Yang-Markidis (arXiv:2605.16074, ICS Workshops '26) — empirical recoverability
  via ML features; same noise propagation structure, complementary approach.
- Tight success bounds for period finding (arXiv:2506.20527, quant-ph/0607148) —
  noise-free bounds; our closed form is the noisy extension.
- Coherence/decoherence in noisy Shor (arXiv:2508.11962) — lower bounds; our
  approach gives the explicit functional form.
- Quantum SR for channels (quant-ph/9903062, 1109.4147) — general quantum SR;
  we narrow to Shor-class algorithms.

## Citing

If you use this work, cite both v0.2.1 (for theorems 1-5 + raw data) and v0.3.0
(for the §3.6 self-correction and closed-form framework):

```
@misc{shor_v0_3_0,
  title  = {A noise-invariant determinism theorem for multi-base post-processing
            in Shor's order finding (v0.3.0, §3.6 errata)},
  author = {{Hashevolution}},
  year   = {2026},
  doi    = {TBD},
  note   = {v0.3.0 supersedes §3.6 interpretation of v0.2.1
            (10.5281/zenodo.20681847)}
}
```
