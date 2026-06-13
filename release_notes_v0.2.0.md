# v0.2.0 — Multi-boundary mechanism observation

**Date**: 2026-06-13

This release upgrades the §3.6 stochastic-resonance observation from a single-cell "Goldilocks" framing to a **universal trial-level mechanism observation** based on 13-seed × 12-σ high-statistics measurement at `(N, d) = (437, 4)`.

The main theoretical contribution — Theorems 1–5 on noise-invariant determinism, logarithmic coverage time, exact noise scaling, conditional Regev compatibility, and a hybrid (C) + Regev b-trick factoring algorithm — is unchanged. The mechanism observation has been substantially deepened.

## What's new

### §3.6 — Universal trial-level boundary-flip mechanism

- **13 independent base sets × 200 trials × 12 σ values = 31,200 trial-measurements** at (N, d) = (437, 4)
- **13/13 seeds exhibit boundary-flip mechanism** — universal at the base-set level
- **K-bin boundary distribution**: 76.9% K=1/K=2, 15.4% K=2/K=3, 7.7% K=3↔K=1 long-jump
- **σ-curve direction asymmetry** identified: positive-direction seeds saturate + decline, negative-direction seeds monotonically worsen (consequence of K-distribution skew at this cell)
- **Direction independence from K_baseline**: seeds with identical K_baseline = 1.720 can show opposite SR directions, confirming that direction is determined by base-set-specific K-distribution structure, not by aggregate K_mean
- **Mechanism follows the classical Benzi–Buchleitner stochastic resonance shape**: sub-threshold + saturation plateau + overload decline
- **Cross-cell verification**: ceiling cells (4) and noise-floor cell (1, multi-seed) and active-boundary cell (1, multi-seed) are all consistent with regime predictions

### ★ Mechanism diversification at higher K_baseline (new — (1147, 2))

A 5-seed σ-scan at (N, d) = (1147, 2) with K_baseline ≈ 2.92 reveals a *richer* mechanism than the K = 1 / K = 2 boundary flip at (437, 4):

- **High-K rescue** ★ : trials at very high K bins (K = 8, 11, 15, even K = 20) move to *moderate* K bins (K = 4, 5) under noise (seed 3: K = 15 → K = 5, SR = +9.44%; seed 4: 3 × (K = 8 → K = 4), SR = +8.56%)
- **Per-seed |SR| amplifies** with K_baseline: max +9.44% at (1147, 2) vs +1.93% at (437, 4)
- **Boundary distribution shifts**: only 40% K = 1 / K = 2 at (1147, 2) vs 77% at (437, 4); the remainder includes K = 8 → K = 4, K = 15 → K = 5, etc.
- **Cross-seed mean**: +3.35% (sd 5.37, SE 2.40, t = 1.39, p ≈ 0.12 with t-distribution / 0.082 normal approx) — borderline at 5-seed sample, dominated by high-K-rescue seeds

### ★ Engineered amplification — (C) augmentation as a noise buffer

Two thinned variants of the hybrid test the structural prediction that SR magnitude is bounded by the *borderline-trial population*:

- *Over-thinned* (smallest convergent only + no (C)): K_baseline = 19.87 (9.55× sub-functional), **SR = 0.00%** — sub-functional alone does not amplify (no borderlines left)
- *Mild thinned* (all convergents + no (C) augmentation): K_baseline = 2.92 (1.4× sub-functional), **per-seed |SR| amplifies to 4.03–4.44% (~5× larger than full hybrid)**

The (C) augmentation acts as a *noise buffer*: its removal exposes the underlying boundary-flip mechanism in amplified form, while leaving the direction stochastic.

### ★ Algorithm-structure regime map for noise-as-resource (testable conjecture)

| Algorithm structure | SR magnitude | Source |
|---|---|---|
| Single-base Shor (1994) | small (≤ 1%) | **predicted** |
| Multi-base Regev (LLL) | negative SR | **predicted** |
| Hybrid full (this work) | +0.14% at (437, 4) | **measured** |
| Hybrid mild-thinned | ~5× amplified | **measured** |
| Hybrid over-thinned | zero | **measured** |

The regime map predicts that noise-as-resource in quantum factoring is naturally maximized in a *multi-base + per-coordinate independent recovery + no buffer* variant. Standalone Shor and Regev predictions remain testable with the same protocol; we leave their verification to follow-up work.

### Statistical caveat

- **Net SR direction across 13 seeds at (437, 4): mean = +0.144%, t = 0.51, p = 0.31** — *not* statistically significant
- Sign test: 8/13 positive (p = 0.29)
- The mechanism is universal; the net direction is base-set-stochastic at our sample size
- (1147, 2) 5-seed mean +3.35% (p ≈ 0.12 t-dist / 0.082 normal approx) is marginal; high-K-rescue seeds dominate

### Retractions

Six earlier claims have been retracted with explicit footnotes:

1. *17.86% peak* at (1147, 2, σ=0.01) — single-seed direction fluctuation
2. *Polynomial scaling SR ∝ N^α* — rejected when N=2491 cells gave near-zero or negative SR
3. *σ_opt ∝ N^α* ("small lock, small wiggle" intuition) — rejected by σ-scan finding σ_opt ≈ 0.010 independent of N
4. *Anti-Optimization Principle (d=1 universal positive SR)* — undermined by multi-seed re-measurement at (1147, 1) giving −0.53% ± 4.28%
5. *V3 sign-test p = 0.03 as significance* — invalid because σ values within the saturation plateau are perfectly correlated (same boundary trials flip in all)
6. *Goldilocks single-cell robust* — refined: K_baseline ≈ 2 marks the regime where mechanism is detectable, but direction is stochastic, not systematically positive

### Reproducibility

New scripts to reproduce §3.6 results:

```
python -m experiments.sigma_scan_437                # (437, 4) baseline σ-scan (3 seeds × 12 σ × 200 trials)
python -m experiments.sigma_scan_437_extend         # (437, 4) extended seeds 4-13 + K-histogram backfill
python -m experiments.analyze_histograms            # per-seed K-bin flip identification
python -m experiments.sigma_scan_general 1147 2 5 100 compact  # (1147, 2) cross-cell σ-scan
python -m experiments.sr_amplification              # over-thinned (null amplification result)
python -m experiments.sr_mild_amplification         # mild thinned (5× amplification)
```

Raw measurement data is committed:
- `experiments/sigma_scan_437_d4_results.txt` — K_mean at (437, 4)
- `experiments/sigma_scan_437_d4_extended.txt` — extended seeds K_mean
- `experiments/sigma_scan_437_d4_histograms.txt` — per-seed K-histograms at (437, 4)
- `experiments/sigma_scan_N1147_d2_compact_results.txt` — (1147, 2) K_mean
- `experiments/sigma_scan_N1147_d2_compact_histograms.txt` — (1147, 2) K-histograms
- `experiments/sr_mild_amplification_results.txt` — mild thinned amplification data

All scripts are resumable (immediate per-cell save, skip-existing on re-run).

## Unchanged from v0.1.0

- Theorems 1–5 (paper §3.1–3.5) with proofs and empirical verification
- Lemma 5.1 (per-`b` nontrivial-sqrt probability ≥ 1/2 for any semiprime)
- 17,700 measurements across 6 composite sizes verifying Theorem 1
- Cross-cell verification of Theorem 5 hybrid algorithm at N ∈ {437, 1147, 2491, 4087}
- Hardware-calibrated noise simulation (Appendix E)

## Connection to noise-as-resource literature

The mechanism we identify is the discrete analog of stochastic resonance in continuous systems. It provides the first explicit, mechanism-level bridge between **integer factoring quantum algorithms** and the **noise-as-resource** paradigm (Benzi et al. 1981 on classical SR; Wellens–Buchleitner 2004 on quantum SR; Plenio–Huelga 2008 on ENAQT).

The effect magnitude (1–7 trial flips per 200 trials, ±0.3–2% K change per seed, no statistically significant cross-seed net direction) is too small to enable cryptographic advantage; the contribution is conceptual.

## Open questions / future work (v0.3 candidates)

1. **Direction bias at larger sample sizes** (30+ seeds) — would a small net positive bias emerge?
2. **Active-boundary determinants** at the base-set level — what predicts which K-boundary flips?
3. **Extended (1147, 2) verification** — script `sigma_scan_N1147_d2_extended.py` (3 seeds × 7 σ × 100 trials, ~2 hours) tests high-K rescue universality at seeds 2, 3, 6
4. **Algorithm-structure regime map verification** — measure SR at standalone Shor and Regev with LLL post-processing; confirms or refines the predictions in §3.6
5. **Universality across N** — verification at (1147, 3), (4087, 4) at V3-style scale
6. **Long-jump events** — true noise-induced trajectory divergence or cascade through adjacent boundaries?
7. **Hardware verification** — survives structured (non-iid) phase noise?

## How to cite

See `CITATION.cff` (GitHub will auto-format BibTeX from this file).

```
Hashevolution. (2026). A Noise-Invariant Determinism Theorem for Multi-Base
Post-Processing in Shor's Order Finding (Version 0.2.0). Zenodo.
https://doi.org/[TO BE ADDED]
```

## License

MIT (see LICENSE)

## Companion code

~1,000 lines of numpy. No quantum libraries required. Verified on Python 3.13.

## Acknowledgements

This work benefited from extensive iteration: initial single-seed observations gave way to multi-seed validations that overturned several over-claims (see Retractions above). The published narrative reflects the corrections, not the original mistakes.

---

🤖 Generated with [Claude Code](https://claude.com/claude-code)
