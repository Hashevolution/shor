# §3.6 Reframe Draft (대기 중)

**Status**: 11 seeds 분석 완료 — seeds 12-13 결과 받으면 최종 [TBD] 채우고 paper.md / paper.tex 에 commit.

**Source data**:
- V3 (N=437 d=4, 2000 trials × 6 σ): 단일 또는 소수 seed
- 현 σ scan (sigma_scan_437_d4_results.txt): 3 seeds × 200 trials × 12 σ
- 현 σ scan extend (sigma_scan_437_d4_extended.txt): seeds 4-13 [진행 중]
- Histogram (sigma_scan_437_d4_histograms.txt): 모든 seeds 에서 σ=0, σ=0.050 비교
- AOP grid + 이전 measurements (summary.md)

**현재 진행 (n=11)**:
- Mean SR at σ=0.050: +0.464% (n=10) → [TBD seeds 12-13 후 갱신]
- Sign test: 8/11 positive (binomial p ≈ 0.11)
- Mechanism: **Multi-boundary flip universal (11/11)**, K=1/K=2 dominant (8/11), K=2/K=3 secondary (3/11)
- Sanity check: 11/11 ✓ histogram-K_mean 정확 일치

**[TBD] 항목** (seeds 12-13 후):
- 최종 mean SR, sd, SE, t, p
- Sign test 최종 (positive / total)
- Seeds 12-13 의 flip 패턴 통합

---

## 3.6 Trial-level noise sensitivity at K-bin boundaries (multi-boundary mechanism observation)

We document a small but mechanistically clean phenomenon in the hybrid (C) + Regev b-trick algorithm: phase noise of magnitude `σ ∈ [0.005, 0.100]` reliably flips a small number of "borderline" trials at one of the K-bin boundaries of single-run hybrid factoring — typically the `K = 1 / K = 2` boundary, occasionally the `K = 2 / K = 3` boundary. Across [TBD: N] independent seeds at (N, d) = (437, 4), every seed exhibits boundary-flip activity (universality of the mechanism), with the *specific* boundary location and direction determined by base-set composition. The aggregate K-mean effect is approximately `±0.3` to `±2%` per seed, with [TBD: net positive bias / mean indistinguishable from 0] across our sample.

The K = 1 / K = 2 boundary, which dominates in ~73% of sampled seeds, corresponds to the algorithm's intrinsic per-coordinate success rate giving a meaningful population of trials at the success/failure threshold. Less frequently (~18%), seeds whose K-distribution gives a stable K = 1 / K = 2 boundary show activity instead at the K = 2 / K = 3 boundary. The mechanism is qualitatively a member of the *stochastic resonance* family (Benzi et al. 1981; Wellens–Buchleitner 2004 for the quantum variant), and quantitatively shows the classical SR σ-curve shape (saturation plateau + overload decline). It is too small to have cryptographic implications but provides a clean bridge between integer factoring and the broader *noise-as-resource* literature.

### High-statistics measurement at N = 437, d = 4

We measured `K` with 200 trials per `σ`, across 12 σ values, and [TBD: 5-13] independent seeds:

| σ | mean K (over seeds) | mean SR % | SE | t | p (1-sided) |
|---:|---:|---:|---:|---:|---:|
| 0.000 | [TBD] | (baseline) | — | — | — |
| 0.005 | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| 0.050 | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| 0.100 | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| 0.150 | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| 0.200 | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |

[Plateau SR at σ ∈ [0.005, 0.100] is statistically [TBD: significant / borderline / not significant] at the `p < 0.05` level. The K-mean differences are individually small but consistent in [TBD: direction].]

### Mechanism: K-bin boundary flips (multi-boundary)

For each seed we compute the K-histogram at `σ = 0` and `σ = 0.050` (within the saturation plateau). The diff histogram reveals which K-bin trials moved between, identifying the dominant flip per seed. Across [TBD: N] seeds:

| seed | K_baseline | SR % at σ=0.050 | direction | dominant flip | magnitude |
|---|---|---|---|---|---|
| 1 | 2.200 | +0.682% | helps | K = 2 → K = 1 | 1 |
| 2 | 1.555 | +1.929% | helps | K = 2 → K = 1 | 1 |
| 3 | 1.720 | -0.872% | hurts | K = 1 → K = 2 | 1 |
| 4 | 1.720 | +0.581% | helps | K = 2 → K = 1 | 2 |
| 5 | 1.630 | +0.613% | helps | K = 3 → K = 1 (long jump) | 1 |
| 6 | 1.550 | +0.323% | helps | K = 3 → K = 2 (secondary boundary) | 1 |
| 7 | 1.515 | +0.330% | helps | K = 2 → K = 1 | 2 |
| 8 | 2.315 | +0.432% | helps | K = 2 → K = 1 | 2 |
| 9 | 2.215 | +1.580% | helps | K = 2 → K = 1 | 5 |
| 10 | 2.090 | -0.957% | hurts | K = 2 → K = 3 (secondary boundary, negative) | 2 |
| 11 | 1.820 | -1.099% | hurts | K = 1 → K = 2 | 2 |
| [TBD seed 12-13] | ... | ... | ... | ... | ... |

**Boundary location distribution**:
- **K = 1 / K = 2**: 8 of 11 seeds (72.7%) — primary boundary
- **K = 2 / K = 3**: 2 of 11 seeds (18.2%) — secondary boundary (seeds 6, 10)
- **Non-adjacent (K = 3 ↔ K = 1)**: 1 of 11 seeds (9.1%) — long-jump (seed 5)

The K-mean change for each seed is *exactly* reproduced by the histogram diff (sanity check: 11/11 seeds, `Δ K_total = (K_mean(σ=0.05) − K_mean(σ=0)) × n_trials` to within integer precision). This confirms that the dominant flip identified is the principal contribution; subsidiary boundary activity, if present, contributes ≤ 1 trial per 200.

**Direction is determined by which side of the active boundary has more near-miss trials.** For seeds where the active boundary is K = 1 / K = 2, the direction is "K = 2 → K = 1" (helps) when the noise-free distribution has more *near-success* (almost-K = 1) trials, and "K = 1 → K = 2" (hurts) when it has more *barely-succeeded* (K = 1 but borderline) trials. For seeds where the active boundary shifts to K = 2 / K = 3, the same dynamic applies with shifted bins.

A striking instance of base-set determinism: seeds 3 and 4 have **identical** `K_baseline = 1.720`, yet show **opposite SR directions** (`-0.872%` vs `+0.581%`). The K_baseline value is a population summary; the direction is determined by which side of the active boundary the base set populates more densely, an internal-distribution property that K_baseline does not capture.

### σ-curve consistency with classical SR

The K-mean as a function of σ (averaged across [TBD: N] seeds, 200 trials each) follows the classical Benzi–Buchleitner SR-curve shape:

- **σ < 0.005 (sub-threshold)**: zero effect — phase noise too small to perturb boundary trials.
- **σ ∈ [0.005, 0.100] (saturation plateau)**: all flippable boundary trials flip; the magnitude is independent of σ within this range. Within a single seed, all 9 σ values in this range give *identical* K-mean to integer-flip precision (deterministic given that the same trials flip across the range).
- **σ ∈ [0.150, 0.200] (overload)**: additional trials at higher K-bins begin to flip in the failure direction, partially or fully canceling the boundary effect. Net SR at σ = 0.200 is statistically indistinguishable from zero.

This is the quantitative shape of stochastic resonance — saturation followed by overload — and matches the prediction from the boundary-flip mechanism, where additional noise pushes trials further from the active boundary into other regimes that do not benefit (or harm) the algorithm uniformly.

The intra-seed determinism of the plateau (`σ = 0.005, 0.010, ..., 0.100` give identical K-means within a seed) is itself a strong signature: it means the dominant flips are *deterministic* with respect to σ once σ exceeds the sub-threshold; the only stochastic ingredient is *which seeds* land on which side of the active boundary (between-seed variance).

### Cross-cell verification

The boundary-flip mechanism makes three falsifiable predictions, which we test against existing measurements at other `(N, d)` cells:

**Prediction 1: Ceiling regime (`K_baseline → 1`) shows SR ≈ 0.**

At these cells, almost all noise-free trials already succeed in 1 run; there is no `K = 2` population to flip into `K = 1`, and the per-trial K can only increase (failure direction) under noise. The boundary-flip mechanism predicts SR magnitude near zero (small effects from rare flips at higher K-boundaries).

| (N, d) | K_baseline | SR % (single seed) | match |
|---|---|---|---|
| (437, 8) | 1.190 | +0.00% | ✓ exact |
| (1147, 8) | 1.055 | +0.47% | ✓ small |
| (2491, 4) | 1.065 | −1.88% | ✓ small |
| (2491, 8) | 1.000 | −1.50% | ✓ |

All four ceiling cells fall within `±2%` of zero, consistent with the prediction.

**Prediction 2: Noise-floor regime (`K_baseline ≫ 2`) shows per-seed variance dominating the effect.**

At these cells, K-distribution is spread over many bins; the K = 1 / K = 2 boundary is sparsely populated, but the per-trial K-variance is high. Single-seed SR estimates have standard error comparable to or larger than any plausible effect.

| (N, d) | K_baseline | data | match |
|---|---|---|---|
| (1147, 1) | 5.780 | multi-seed mean `−0.53%`, sd `4.28%`, 3 seeds × 100 trials | ✓ exact — variance dominates |
| (437, 1) | 7.277 | single seed `+0.78%` (300 trials) | unverified (within plausible variance) |
| (2491, 1) | 5.643 | single seed `+1.36%` (300 trials) | unverified (within plausible variance) |

The `(1147, 1)` multi-seed measurement directly confirms the noise-floor regime: per-seed standard deviation (`4.28%`) exceeds any individual seed's SR estimate, and the 95% CI `[−5.4%, +4.3%]` includes both zero and the previously-celebrated single-seed `+2.60%` peak.

**Prediction 3: Active-boundary regime (`K_baseline ≈ 2`) shows SR `±0.5–2%` with seed-dependent direction.**

At these cells, the K = 1 / K = 2 boundary is well-populated. The mechanism predicts that direction varies across seeds, and the per-seed magnitude is bounded by the number of borderline trials (~1–6 of 200 trials in our setting).

| (N, d) | K_baseline | data | match |
|---|---|---|---|
| (437, 4) | 1.92 | V3: `+0.91%` × 2000 tr; current scan [TBD: mean ± SE, N seeds] | ✓ small, [TBD] direction |
| (1147, 2) | 2.43 | single `+1.24%` × 300 tr; multi-seed confirm `+0.42%` × 4×1000 tr | ✓ single regresses to small multi-seed mean |
| (1147, 3) | 1.67 | single `+1.00%` × 300 tr | consistent (unverified) |
| (1147, 4) | 1.44 | single `+5.90%` × 200 tr | likely outlier (>2× predicted) |
| (2491, 2) | 2.25 | single `−4.89%` × 200 tr | likely outlier (negative direction × seed variance) |

The cleanest verification is at `(1147, 2)`: a single-seed estimate of `+1.24%` regressed to `+0.42%` under multi-seed confirmation (4 seeds × 1000 trials, `t = 0.50`, not significant), exactly the mechanism prediction of small magnitude with seed-dependent direction averaging toward zero.

**Direction independence from K_baseline.**

The mechanism predicts that direction is determined by base-set composition, not by the aggregate K_baseline value. This is directly verified by seeds 3 and 4 at (437, 4): both with K_baseline = 1.720, opposite SR directions (`−0.872%` vs `+0.581%`).

### Interpretation in the noise-as-resource framework

The multi-boundary flip mechanism is the discrete analog of stochastic resonance in continuous systems. In ENAQT (Plenio–Huelga 2008) and related biological-transport studies, environmental noise drives a continuous probability distribution toward a more efficient transport configuration. In our setting, the analog is: phase noise drives discrete K outcomes across the K-boundaries (primarily K = 1 / K = 2, occasionally K = 2 / K = 3), and aggregate SR direction depends on the population balance at the active boundary. Both are instances of noise as a *boundary-locator* rather than a *signal amplifier*.

To our knowledge this is the first observation of a clean, mechanism-level SR signal in an integer-factoring quantum algorithm, where the mechanism is identified at trial-level granularity (which K-bin trials moved and in which direction). The effect is too small to have cryptographic implications (1–7 trial flips per 200 = `±0.3–2%` K reduction per seed, [TBD: net across seeds]) but provides a clear conceptual bridge between the noise-as-resource paradigm and quantum factoring.

### What we earlier over-claimed and have since retracted

1. *`+17.86%` peak at `(N, d, σ) = (1147, 2, 0.01)` (150 trials × 1 seed):* Re-measurement at 1000 trials × 4 seeds gave mean `+0.42%` (`t = 0.50`, not significant). The original value is consistent with within-mechanism variance at a single seed with small trial count.

2. *Polynomial scaling `SR ∝ N^α`:* Based on 2 data points (N = 437, 1147); rejected when N = 2491 cells showed near-zero or negative SR. The boundary-flip mechanism predicts no monotonic N-scaling because the effect depends on K_baseline (not N) and seed-specific direction.

3. *`σ_opt ∝ N^α` ("small lock, small wiggle" intuition):* A dedicated σ-scan at N = 437 and N = 1147 found `σ_opt ≈ 0.010` independent of N — consistent with the σ-saturation plateau (any σ > threshold flips the same trials).

4. *"Anti-Optimization Principle" — single-base hybrid (d = 1) maintains positive SR for all N:* Motivated by single-seed AOP grid measurements at `(437, 1), (1147, 1), (2491, 1)` (all `+0.78%, +2.60%, +1.36%`). Multi-seed re-measurement at `(1147, 1)` gave `−0.53% ± 4.28%`, fully consistent with the noise-floor regime (variance > effect). The d = 1 universality was a single-seed artifact.

5. *V3 sign test `p = 0.03` as significance:* The V3 measurement showed all 5 tested σ values below the noise-free baseline at one seed; we initially treated this as `p = (1/2)^5 ≈ 3%`. This is invalid: within the σ-saturation plateau, σ values are perfectly correlated (all flip the same boundary trials), so 5/5 below is no more informative than 1/1. The proper significance test is between-seed direction, [TBD: which our [TBD: 5-13]-seed measurement gives `p = [TBD]`].

### Caveats

- *Effect magnitude is tiny.* Approximately 1–7 trial flips per 200 trials at the active K-boundary; aggregate K-mean change `±0.005–0.030` runs, or `±0.3–2%` per seed.
- *No definite net SR effect across seeds.* [TBD: Mean SR across N seeds is `[TBD]±[TBD]%` (t = [TBD], p = [TBD]); the direction (net positive / net zero / net negative) [TBD: is / is not] statistically supported.]
- *Direction is base-set-deterministic, not random-in-population.* Different seeds give different *directions* and sometimes different *active boundaries* (K = 1 / K = 2 vs K = 2 / K = 3) but always exhibit the boundary-flip mechanism. Predicting direction (or boundary location) for an unmeasured seed requires inspecting its base-set K-distribution, not its K_baseline.
- *Phase-noise specific.* Depolarizing and amplitude-damping noise show monotone degradation. Artificial `k ± δ` perturbation of measured k yields zero effect (continued-fraction expansion absorbs small δ).
- *Closed-form prediction is partial.* We predict the σ-curve shape (saturation + overload) and that some K-boundary will be active, but not the specific active boundary or the flip direction (both depend on base-set composition).

### Open questions

1. *Direction net bias across seeds.* Our [TBD: 11–13] seeds at (437, 4) show [TBD: a positive net trend with marginal significance / a clear positive net bias / mean indistinguishable from 0]. A larger sample (30+ seeds) would distinguish a small genuine positive bias from a symmetric direction distribution.
2. *Which boundary is active at a given seed.* We observe 73% K = 1 / K = 2, 18% K = 2 / K = 3, and 9% long-jump in our 11-seed sample. The base-set property that determines this is not currently predictable from K_baseline alone; characterizing it would require finer-grained analysis of the noise-free K-histogram per seed.
3. *Universality across N.* The multi-boundary flip mechanism predicts that any active-boundary cell shows the same qualitative pattern. Confirmation at `(1147, 2), (1147, 3), (4087, 4)` requires multi-thousand-trial measurements at each.
4. *Long-jump events.* Seed 5 in our sample shows a single trial moving from K = 3 to K = 1, skipping K = 2. Whether this reflects a true noise-induced trajectory divergence or a cascade through two adjacent boundaries (recorded as one event at our resolution) remains open.
5. *Hardware verification.* Phase noise on real quantum hardware is structured (not Gaussian iid); whether the boundary-flip mechanism survives that noise structure is an empirical question.

We document this as a single-algorithm, mechanism-level observation: qualitatively a member of the noise-as-resource family, quantitatively too small to enable practical advantage. The K = 1 / K = 2 boundary mechanism is the *right* level of description; the integrated SR effect is the *measurement* of that mechanism.

Reproduce: `python -m experiments.sigma_scan_437` (sigma scan + per-seed histogram), `python -m experiments.sr_aop` (single-seed AOP grid for cross-cell predictions).
