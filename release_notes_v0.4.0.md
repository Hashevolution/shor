# Release notes — v0.4.0 (2026-06-14)

**DOI**: minted by Zenodo on release (backfilled after archival).

## Summary

**v0.4.0** is a *formalization + honesty* release on top of v0.3.0 (DOI
[10.5281/zenodo.20685015](https://doi.org/10.5281/zenodo.20685015)).

It adds **Theorem 6**, a sign-agnostic no-go result that closes the
stochastic-resonance (SR) factoring question raised in §3.6, backs it with a
cryptographic-regime N-scaling measurement, and downgrades the Yang-Markidis
positioning from "analytical complement / gap-filling" to the honest
"numerical verification + boundary mapping."

## What changed

### Theorem 6 — SR no-go (paper §3.3.ter, new)

For any noise in the coherence-loss class (★), `g(ε) = (1−ε)·g_0 + ε·g_∞`
with `0 < g_0, g_∞ ≤ 1`, **without assuming the sign of `g_0 − g_∞`**:

1. **No tuned resonance.** `E[K(ε)]` is monotone in `ε`; the optimum sits at an
   endpoint, never an interior `ε*`. There is no noise level to "tune into."
2. **Closed-form swing.** `|ΔK| = E[K_λ^ideal] · |1/g_∞ − 1/g_0|`, fixed by the
   two endpoint probabilities alone, independent of the σ-profile.
3. **No asymptotic advantage.** If the reciprocal gap is `O(1)`, the swing is
   `O(log log N)` — a constant multiplicative factor, never a speedup.

The earlier draft wrongly assumed `g_∞ ≤ g_0` ("noise always hurts"). Live
N-scaling data falsified that: the smallest-order setup at each N is a genuine
positive-SR cell (`g_∞ > g_0`). Theorem 6 is therefore stated sign-agnostically.

### N-scaling measurement (paper §3.6.bis, new table)

`|Δ| = |p_0 − ρ|` over 12 setups (3 per N, 300–1,000 MC samples,
`shor_n_scaling.py`):

| N | mean p_0 | mean ρ | mean \|Δ\| | min \|Δ\| | max \|Δ\| |
|---:|---:|---:|---:|---:|---:|
| 437 | 0.741 | 0.410 | 0.446 | 0.173 | 0.607 |
| 1147 | 0.443 | 0.402 | 0.211 | 0.118 | 0.260 |
| 2491 | 0.486 | 0.287 | 0.372 | 0.260 | 0.483 |
| 4087 | 0.476 | 0.390 | 0.259 | 0.220 | 0.297 |

Across a 9× range in N the gap stays in `[0.12, 0.61]` with no monotone growth
toward 1 — the empirical regularity behind Theorem 6's `O(1)` reciprocal-gap
assumption. The result is reproducible to the digit under fixed seeds.

### Honest reframing of Yang-Markidis positioning

- Their two-stage model is their **Eq. (3), §5**; the mixing weight `ε` is a
  *conceptual* parameter ("total weight transferred out of the intended
  family"), left unspecified — neither fitted nor given as `ε(σ)`. Full text
  (incl. appendix) verified 2026-06-14.
- The `exp(−σ²)` decay is the standard dephasing result (Nielsen–Chuang §8.3).
  **We claim no new mechanism.** Our contribution is narrow and verificational:
  (i) verify their qualitative weight equals the standard-dephasing factor
  (R²=0.95 at N=437); (ii) map the boundary (breaks for amplitude damping,
  R²=0.03); (iii) confirm across five algorithm classes.
- README / arxiv_draft / sr_sigma_curve_model wording corrected accordingly.
- We **do not** position this as filling an analytical gap they could not close.

### Canonical source declaration

- `paper.md` is now the **single canonical source of truth**.
- `paper.tex` is **deprecated / frozen at v0.2.1** — to be regenerated wholesale
  from `paper.md` only at arXiv-submission time (banner added; see
  `compile-notes.md`).

## Retained (unchanged)

- Theorems 1–5 and all raw measurement data.
- The v0.3.0 closed-form σ-curve and five-algorithm verification.
- The conclusion that SR-based factoring acceleration is precluded — now
  *formalized* as Theorem 6 rather than asserted.

## What this is NOT

- Not a new quantum advantage, not a new mechanism.
- Not an unconditional impossibility result — Theorem 6 is a no-go *within the
  framework and under its measured regularities* (see its "Scope and honest
  caveats").

## Citing

If you use this work, cite v0.4.0 for Theorem 6 + N-scaling + honest framing,
alongside v0.3.0 (closed-form framework) and v0.2.1 (theorems 1–5 + raw data):

```
@misc{shor_v0_4_0,
  title  = {A noise-invariant determinism theorem for multi-base post-processing
            in Shor's order finding (v0.4.0, Theorem 6 SR no-go)},
  author = {{Hashevolution}},
  year   = {2026},
  doi    = {10.5281/zenodo.XXXXXXXX},
  url    = {https://doi.org/10.5281/zenodo.XXXXXXXX},
  note   = {v0.4.0 adds Theorem 6 and N-scaling to v0.3.0
            (10.5281/zenodo.20685015)}
}
```

(DOI backfilled once Zenodo archives the release.)
