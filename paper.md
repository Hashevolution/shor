# A Noise-Invariant Determinism Theorem for Multi-Base Post-Processing in Shor's Order Finding

*Draft, 2026-06-11 (v0.2.1); §3.6.bis self-correction added 2026-06-14 (v0.3.0).*
*Companion code: [github.com/Hashevolution/shor](https://github.com/Hashevolution/shor).*

**Current release: v0.3.0**
[![DOI v0.3.0](https://zenodo.org/badge/DOI/10.5281/zenodo.20685015.svg)](https://doi.org/10.5281/zenodo.20685015) (v0.3.0, §3.6 self-correction + closed-form σ-curve framework)
[![DOI v0.2.1](https://zenodo.org/badge/DOI/10.5281/zenodo.20681847.svg)](https://doi.org/10.5281/zenodo.20681847) (v0.2.1, theorems 1–5 + §3.6 original framework)

**Reading note**: §3.6 of this paper describes the original "boundary-flip mechanism" framework (v0.2.1). §3.6.bis (v0.3.0) replaces that interpretation with a single closed-form σ-curve `p(σ) = ρ + (p_0 − ρ)·exp(−σ²)` verified across five algorithm classes (Grover, Shor, QPE, Simon, Hybrid; R² ∈ [0.88, 0.99]). The §3.6 measurement data is retained unchanged; the *interpretation* is corrected in §3.6.bis. See [release_notes_v0.3.0.md](release_notes_v0.3.0.md) and [sr_sigma_curve_model.md](sr_sigma_curve_model.md).

## Abstract

We give five theorems, a v0.3.0 generalization (Theorem 3') with a no-go corollary (Theorem 6), and one conditional compatibility observation about classical post-processing of multiple-base order-finding measurements in Shor's algorithm. (1) **Noise-invariant determinism**: maintaining an accumulated exponent candidate `L` — the least common multiple of orders recovered from previous bases — and augmenting standard continued-fraction post-processing with a divisor search over `L` yields a procedure (which we call (C)) that recovers the order `r_a` of any new base `a` deterministically whenever `r_a | L`, *independent of the measurement distribution*. (2) **Logarithmic coverage time (ideal)**: for a semiprime `N = pq`, the expected number of independent uniform bases `K_λ` required for `L = λ(N)` satisfies `E[K_λ] ≤ 1 + Σ_{ℓ | λ(N)} 1/(ℓ^{s_ℓ} - 1)`, where `s_ℓ ∈ {1, 2}` records the ℓ-Sylow overlap of `(Z/p)*` and `(Z/q)*`. (3) **Noise scaling**: under a class of "destructive" noise models (depolarizing, bias, modexp), the actual algorithm K_λ scales exactly as `E[K_λ^{alg}(η)] = E[K_λ^{ideal}] / g_M(η)` where `g_M(η)` is the per-base extraction probability. Together: a noise-adjusted logarithmic number of measurements suffice to enter the noise-immune regime of (1). We verify (1) across 17,700 measurements (zero violations), (2) across 17,000 trials (mean within bound), and (3) across 9 noise setups on N=437 (mean error ~11% for the destructive class). Theorem 4 demonstrates conditional compatibility with Regev's 2023 multi-base framework (numpy simulation, 200 trials × 4 N), provided each measurement coordinate's marginal is Shor-like. Theorem 5 introduces the **(C) + Regev b-trick hybrid**: applying (C) per coordinate while immediately checking each newly-recovered `ord(a_i)` for a nontrivial square root via `b_i^{ord(a_i)}`. Empirically, the hybrid factors N=437 in 1.03-1.70 runs (100% success across 5 noise/corruption conditions × 30 trials), versus 6.7-7.1 (70%) for (C) alone and 3.1-4.0 (90%) for pure b-trick. We formalize this with **Lemma 5.1** (per-`b` nontrivial-sqrt probability ≥ 1/2 for any semiprime) and a closed-form bound `E[K] ≤ 1/(1 − (1 − g(η)·c)^d)` matching the N=437 empirical mean within 4% (the bound is correspondingly loose at larger N, where the per-coordinate recovery rate `g` drops below 1; the empirical mean stays in `[1.1, 1.5]` across `N ∈ {1147, 2491, 4087}`). We discuss the theorem's relation to prior work — Knill's lcm trick (1995), McAnally's larger-Q convergent enumeration (2001), the Bach-Shallit textbook treatment of `r | λ(N)`, and the quantum algorithm for computing the Carmichael function (2021) — and conclude that while every individual ingredient is folklore, the *clean statement of the noise-invariance corollary* appears not to be made explicit in surveyed literature. We additionally show that two natural attempts to improve measurement count beyond this scheme — adaptive base selection and lattice-based joint post-processing — yield no further reduction, suggesting the theorem captures the essential structure of the problem at this scale. Finally, we document a universal trial-level mechanism observation at `(N, d) = (437, 4)`: across 13 independent base sets (200 trials × 12 σ values each), phase noise σ ∈ [0.005, 0.100] reliably flips a small number of trials at one of the K-bin boundaries (primarily K=1/K=2, occasionally K=2/K=3) of single-run factoring — but the direction (success ↔ failure) is base-set-stochastic, and the cross-seed mean SR (+0.14%, t=0.51, p=0.31 at σ=0.050) shows no statistically significant net bias. The mechanism is qualitatively a member of the stochastic-resonance family and provides a conceptual bridge to the noise-as-resource literature (ENAQT) but is too small in absolute terms to enable cryptographic advantage. **Theorem 6** formalizes this: for any noise in the coherence-loss class (★) — *regardless of whether noise helps or hurts* (some cells have `g_∞ > g_0`, a genuine positive-SR effect) — the expected run count is monotone in the noise weight, so no tuned (interior-optimal) resonance exists, and the total swing is the closed-form constant `|1/g_∞ − 1/g_0|`; under a bounded reciprocal-gap regularity (measured across `N ∈ {437, …, 4087}`) this rules out any asymptotic SR-based factoring acceleration within the framework. Finally, a **companion direction** (§9, v0.5.0+) measures the *nonstabilizerness ("magic")* spent across the speedup ladder (Simon $0$ / Grover finite, peak exactly 3 bits / Shor growing), develops a coding-theory specialization of marker-set magic, and identifies oracle-hiding with non-Clifford ($T$) gate cost — crediting the flat-state closed form to Tarabunga–Castelnovo and Shor's magic↔period law to Paviglianiti et al.

## 1. Background and notation

Let `N` be a composite integer, not a prime power, with `N` odd. Define `Q = 2^t` where `t = 2⌈log₂ N⌉` (so `Q ≥ N²`). For `a ∈ (Z/NZ)*`, let `r_a` denote the multiplicative order of `a` modulo `N` (the smallest positive integer with `a^{r_a} ≡ 1 mod N`). Shor's order-finding measures `k ∈ {0, 1, ..., Q-1}` from a distribution concentrated near `j · Q / r_a` for integer `j`. Standard classical post-processing uses the continued-fraction expansion of `k/Q` to recover `r_a`.

Two classical facts are central:

(F1) **Carmichael bound.** For all `a ∈ (Z/NZ)*`, `r_a | λ(N)` where `λ` is the Carmichael function — the exponent of the group `(Z/NZ)*`.

(F2) **Order as minimum.** `r_a` is by definition the smallest positive integer `d` with `a^d ≡ 1 mod N`. Any `d` satisfying this is a positive multiple of `r_a`.

## 2. The (C) procedure

Let `L ∈ ℕ` be a quantity maintained by the algorithm — initially `L = 1`, updated as orders are recovered (see Algorithm 1 below). For a measurement `k ∈ {0, ..., Q-1}` and current `L`, define the candidate pool:

```
candidates(k, L) := convergents(k/Q) ∪ divisors(L)
```

where `convergents(k/Q)` is the set of denominators of continued-fraction convergents of `k/Q` up to `N-1`, and `divisors(L)` is the set of positive divisors of `L`.

The (C) post-processing is:

```
def C(a, N, k, Q, L):
    cands = convergents(k/Q, max_denom=N-1) ∪ divisors(L)
    valid = {d ∈ cands : a^d ≡ 1 mod N}
    if valid == ∅: return 0
    return minimize_order(a, N, min(valid))
```

where `minimize_order(a, N, d)` reduces `d` to its smallest positive divisor still satisfying `a^? ≡ 1` (by removing small prime factors). (See `multi_base.py` for the reference implementation.)

**Algorithm 1 (multi-base accumulation).**

```
L ← 1
for iteration in 1..max_iter:
    pick a ∈ {2, ..., N-1} (random; if gcd(a, N) > 1, return factor)
    if L > 1 and a^L ≡ 1 mod N:
        # fast path: recover r_a classically
        r ← min{d | L : a^d ≡ 1 mod N}
    else:
        # slow path: measure k and apply (C)
        k ← measure(a, N)
        r ← C(a, N, k, Q, L)
    if r == 0: continue
    L ← lcm(L, r)
    attempt factoring with (a, r) via standard Shor reduction or factor_from_exponent(L)
```

## 3. Main theorems

This section states two theorems. Theorem 1 ((C)-determinism) characterizes the noise-invariant region of (C). Theorem 2 (number of bases) quantifies how fast the algorithm enters that region.

### 3.1 Theorem 1: (C)-determinism

**Theorem 1 ((C)-determinism).** If `r_a | L`, then `C(a, N, k, Q, L) = r_a` for every `k ∈ {0, ..., Q-1}`.

**Proof.** Assume `r_a | L`. Then `r_a ∈ divisors(L) ⊆ candidates(k, L)`. By construction, `a^{r_a} ≡ 1 mod N`, so `r_a ∈ valid`. Hence `valid ≠ ∅`. By (F2), every `d ∈ valid` is a positive multiple of `r_a`, so `r_a ≤ d`. Combining with `r_a ∈ valid` gives `min(valid) = r_a`. Finally, `minimize_order(a, N, r_a) = r_a` since `r_a` is already the order. ∎

**Corollary 1 (Noise invariance).** Theorem 1 makes no assumption on the distribution from which `k` is drawn. Hence: under any noise model that affects only the measurement step — including arbitrary distortions of the QFT output distribution — `(C)` returns `r_a` whenever `r_a | L`.

**Corollary 2 (L integrity).** When Algorithm 1 updates `L ← lcm(L, r)`, the value `r` returned by `(C)` is always either `r_a` (success) or `0` (failure). In the success case, `r_a | λ(N)` by (F1), so `L` always remains a divisor of `λ(N)`. Under noise, `L` cannot be corrupted to a value outside the divisors of `λ(N)`.

**Corollary 3 (Failure region).** `(C)` can fail only when `r_a ∤ L`. In that case, success depends on `convergents(k/Q)` containing some positive multiple of `r_a` — the classical "(B) path."

### 3.2 Theorem 2: Number of bases to reach λ(N)

Theorem 1 holds *given* `r_a | L`. We now bound how many independent bases suffice to reach the global condition `L = λ(N)`, after which Theorem 1 applies to every subsequent base.

Let `N = pq` be a semiprime with `p, q` distinct odd primes. Let `a_1, ..., a_K` be independent uniform samples from `(Z/N)*`, with orders `r_1, ..., r_K`, and `L_K := lcm(r_1, ..., r_K)`. Let `K_λ := min{K : L_K = λ(N)}`. For each prime `ℓ | λ(N)` define

```
s_ℓ := |{ξ ∈ {p, q} : v_ℓ(ξ-1) = v_ℓ(λ(N))}|  ∈ {1, 2}
```

where `v_ℓ` is the `ℓ`-adic valuation.

**Theorem 2 (K_λ distribution).** For `K ≥ 1`,

```
P[L_K < λ(N)]  ≤  Σ_{ℓ | λ(N)} ℓ^{-K · s_ℓ}  ≤  ω(λ(N)) · 2^{-K}    (tail bound)
E[K_λ]         ≤  1 + Σ_{ℓ | λ(N)} 1/(ℓ^{s_ℓ} - 1)                  (expectation)
P[K_λ > ⌈log₂(ω(λ(N))/ε)⌉]  ≤  ε    for any  ε ∈ (0, 1)              (high-probability)
```

**Proof.** By CRT, `(Z/N)* ≅ C_{p-1} × C_{q-1}`, and a uniform `a ∈ (Z/N)*` corresponds to an independent uniform pair `(x, y)`.

*Step 1 (cyclic-group valuation distribution).* In cyclic `C_n` with `v := v_ℓ(n)`, the set `{x : ord(x) | m}` equals the unique subgroup of order `m`, for `m | n`. Hence `{x : v_ℓ(ord(x)) ≤ k-1}` is the subgroup of order `max{m | n : v_ℓ(m) ≤ k-1} = n · ℓ^{-(v - k + 1)}`. Therefore for uniform `x`:

```
P[v_ℓ(ord(x)) ≥ k] = 1 - 1/ℓ^{v - k + 1}    (1 ≤ k ≤ v)
```

*Step 2 (per-base miss probability).* Let `v_p := v_ℓ(p-1)`, `v_q := v_ℓ(q-1)`, `v := v_ℓ(λ(N)) = max(v_p, v_q)`. Then `v_ℓ(ord(a)) = max(v_ℓ(ord(x)), v_ℓ(ord(y)))`, and using independence of `x, y`:

```
P[v_ℓ(ord(a)) < v] = P[v_ℓ(ord(x)) < v] · P[v_ℓ(ord(y)) < v] = (1/ℓ)^{s_ℓ}
```

(if `v_p < v`, the x-factor is 1; otherwise 1/ℓ by Step 1 with `k = v`. Likewise for `y`. `s_ℓ` counts how many of {p, q} achieve `v_p = v` or `v_q = v`.)

*Step 3 (tail bound).* By independence across bases, `P[ℓ-component of L_K not covered] = ℓ^{-K · s_ℓ}`. Union bound over primes `ℓ | λ(N)` gives the first inequality. `ℓ^{-K · s_ℓ} ≤ 2^{-K}` (since `ℓ ≥ 2, s_ℓ ≥ 1`) gives the second.

*Step 4 (expectation).* `E[K_λ] = Σ_{K ≥ 0} P[L_K < λ(N)]`. The `K = 0` term is `1` (since `L_0 = 1 < λ(N)`). For `K ≥ 1`, apply the tail bound and interchange sums:

```
E[K_λ]  ≤  1 + Σ_ℓ Σ_{K ≥ 1} ℓ^{-K · s_ℓ}  =  1 + Σ_ℓ 1/(ℓ^{s_ℓ} - 1).
```

*Step 5 (high-probability).* Solve `ω · 2^{-K} ≤ ε` for `K`. ∎

**Asymptotic.** By the Hardy-Ramanujan theorem, `ω(p-1) ≤ (1 + o(1)) log log p` for almost all primes `p`, so for typical semiprimes `ω(λ(N)) = O(log log N)` and `E[K_λ] = O(log log log N)`.

**Empirical verification.** We measure the empirical `E[K_λ]` over 1,000 independent trials per `N`, across 17 semiprimes from `N = 15` to `N = 4087`. In every case, the empirical mean stays within 0.4 of the bound (c), and is often within 0.05–0.1. See `experiments/k_lambda_dist.py` and Appendix D.

**Remark.** Pomerance et al. (2017) prove `e(G) ≤ d + 2.752` for any finite abelian `G`, where `e(G)` is the expected number of uniform samples needed to *generate* `G` and `d` is the maximum Sylow generator count. For `G = (Z/N)*` with `N = pq`, `d ≤ 2` so `e(G) ≤ 4.752`. Since reaching `λ(N)` via lcm is weaker than generating, `E[K_λ] ≤ e((Z/N)*) ≤ 4.752` follows as a corollary. Theorem 2(c) refines this in two ways: (i) by using the explicit ℓ-Sylow structure of `(Z/N)*` rather than a generic abelian bound; (ii) by tailoring to "reach exponent" rather than "generate." For e.g. `N = 4087`, Theorem 2(c) gives `E[K_λ] ≤ 2.475` versus Pomerance et al.'s `≤ 4.752` and the empirical mean of `2.26`.

### 3.3 Theorem 3: Algorithm K_λ under destructive noise

Theorem 2 bounds the *ideal* `K_λ` — the number of independent uniform bases needed for `L = lcm(r_{a_1}, \ldots, r_{a_K}) = λ(N)`, assuming the orders `r_{a_i}` are extracted noise-freely. The full algorithm performs (C) post-processing on noisy measurements, which can fail to recover `r_{a_i}` and waste the corresponding base.

We isolate a class of noise models for which this overhead admits a clean closed form. A noise model `M` is **destructive at rate η** if each measurement is, independently, replaced with probability `η` by a draw from an `a`-independent distribution `D_M`, and otherwise is noise-free. Examples: depolarizing (`D_M` uniform), bias-to-zero (`D_M = δ_{k=0}`). Counter-examples: phase decoherence, amplitude damping (both are functional transformations of the ideal distribution, not independent replacements).

Let `g_0 := P[(C) recovers r_a | noise-free]` and `g_{unif}^M := P[(C) recovers r_a | k drawn from D_M]`. Define `g_M(η) := (1-η) g_0 + η · g_{unif}^M`.

**Theorem 3 (algorithm K_λ for destructive noise).** Assume `M` is destructive at rate `η`. Then `K_λ^{alg}(η)`, the number of bases drawn by Algorithm 1 before `L = λ(N)`, satisfies:

```
E[K_λ^{alg}(η)] = E[K_λ^{ideal}] / g_M(η)
```

In particular, `E[K_λ^{alg}(η)] / E[K_λ^{alg}(0)] = g_0 / g_M(η)`.

**Proof.** Markov chain on states `s ∈ divisors(λ(N))`. At state `s`, per trial:

- With prob `ε_s := P[r_a | s]` (covered): fast path of Algorithm 1 returns `r_a` deterministically (Theorem 1); state unchanged.
- With prob `1 - ε_s` (`r_a ∤ s`): slow path measures and applies (C). Under destructive `M`, recovery succeeds with prob `g_M(η)` (state-independent: the replacement distribution `D_M` is `a`-independent, so its convergent enumeration is independent of `s`). On success, state advances to `lcm(s, r_a)`.

The nontrivial transition rate from `s` is `(1 - ε_s) · g_M(η)`, giving `E[\#\,trials\,in\,s] = 1/[(1-ε_s) g_M(η)]`. By linearity over the path:

```
E[K_λ^{alg}(η)]
  = Σ_s 1/[(1-ε_s) g_M(η)]
  = (1/g_M(η)) · Σ_s 1/(1-ε_s)
  = E[K_λ^{ideal}] / g_M(η).      ∎
```

**Empirical verification (N = 437).** Using `g_M(η)` measured by a separate per-base experiment (`experiments/g_eta.py`) and `E[K_λ^{ideal}] = 1.83` (measured):

| Noise | g_M(η) | Thm 3 prediction | Empirical K_λ^{alg} | Error |
|---|---:|---:|---:|---:|
| noise-free          | 0.380 |  4.82 |  4.97 | +3% |
| depol p=0.1         | 0.340 |  5.38 |  6.09 | +13% |
| depol p=0.3         | 0.262 |  6.98 |  7.40 | +6% |
| depol p=0.5         | 0.226 |  8.10 |  9.19 | +13% |
| depol p=0.7         | 0.156 | 11.73 | 14.11 | +20% |
| bias_zero p=0.5     | 0.188 |  9.73 | 10.80 | +11% |

Theorem 3 predicts K_λ^{alg} within ~11% mean error across 6 destructive setups. The residual error reflects state-dependence in `g_M(s, η)` for `(C)` post-processing (a minor effect for destructive noise but not absent).

**Out of scope: structural noise.** For `phase_sigma` and `amplitude_damp` (which are not destructive in the above sense), the per-base recovery probability `g_M(s, η)` depends strongly on the state `s` (specifically, recovery for larger `r_a` is harder under peak smearing). The analogous prediction `E[K_λ^{ideal}] / g_M(L=1, η)` is a *lower* bound on `E[K_λ^{alg}]`, but underestimates the truth by a factor 1.4–2.4 at moderate noise and >5 at extreme noise. Deriving a closed-form `g_M(s, η)` for structural noises is left as future work — *see §3.3.bis (v0.3.0 extension) for the phase-noise case, which we resolve in closed form*.

### 3.3.bis Universal coherence-loss form (v0.3.0 extension)

Theorem 3 fits the *destructive at rate η* class via `g_M(η) = (1 − η)·g_0 + η·g_{unif}^M`. We now show that this is one specialization of a more general form, derived in v0.3.0 (§3.6.bis) for phase noise via the FFT noise-averaging identity:

```
E[ |FFT(a · e^{iε})_k|² ] = (1 − exp(−σ²)) / Q  +  exp(−σ²) · P_0(k),
```

which gives the per-run success probability under phase noise of magnitude σ as `p(σ) = ρ + (p_0 − ρ)·exp(−σ²)`, where `ρ = |S|/Q` is the uniform success density and `p_0 = Σ_k P_0(k)·I(k)` is the noise-free success probability.

Comparing to Theorem 3: the structure is identical with the substitutions

```
g_0  ↔  p_0,         g_{unif}^M  ↔  ρ,       η  ↔  ε,
```

provided we identify `ε = 1 − exp(−σ²)` for phase noise. This motivates the following generalization:

**Theorem 3' (universal coherence-loss form).** For any noise model `M` for which there exist `ε(M) ∈ [0, 1]` and `g_∞(M) ∈ [0, 1]` such that the noise-averaged per-base success probability factors as

```
g_M = (1 − ε(M)) · g_0  +  ε(M) · g_∞(M),                              (★)
```

Theorem 3's Markov-chain argument yields `E[K_λ^{alg}] = E[K_λ^{ideal}] / g_M`. The Theorem 3 case is `ε = η`, `g_∞ = g_{unif}^M` (destructive at rate η). The phase-noise case is `ε = 1 − exp(−σ²)`, `g_∞ = ρ`. We classify five noise models of paper §3.3 / Appendix E by whether they admit a (`ε`, `g_∞`) pair satisfying (★):

| Noise model        | ε(M)                         | g_∞(M)                       | Form (★) applies | R² (PR #5, 3 setups × 8 levels × 200 trials) |
|--------------------|------------------------------|------------------------------|------------------|----------------------------------------------|
| Phase σ            | `1 − exp(−σ²)`               | `ρ` (uniform success density)| **yes**          | +0.9519 (closed form derivation, §3.6.bis)   |
| Depolarizing p     | `p`                          | `ρ` (uniform mix)            | **yes**          | **+0.9953**                                  |
| Bias zero p        | `p`                          | `I(k = 0)` (often 0)         | **yes**          | **+0.9963**                                  |
| Modexp error q     | `q`                          | success rate on modified S   | **yes** (Theorem 3 case) | (qualitative, in §3.3 table)           |
| Amplitude damping γ| no clean ε                   | no clean g_∞                 | **no**           | +0.0328 (structural; see below)              |

**Phase noise's promotion.** Theorem 3's footnote labels phase σ as "structural" (out of scope). The v0.3.0 derivation in §3.6.bis shows that phase noise admits the form (★) with explicit `ε = 1 − exp(−σ²)` and `g_∞ = ρ`, and the resulting universal form fits R² = +0.95 at N = 437 across five setups × 8 σ × 500 trials. This *moves phase noise from the structural class into the destructive-equivalent class*, and validates Theorem 3's bound on K_λ^{alg} for phase noise (not just the destructive subset originally treated).

**Amplitude damping's boundary status.** The same verification (R² = +0.03) confirms amplitude damping is genuinely *outside* form (★) — the distortion `a_x \mapsto a_x · \exp(−γ x)` is asymmetric across the amplitude index and is not a linear mixture of `P_0` and any single distractor distribution. Deriving an analogous closed form for amplitude damping requires modeling the index-dependent attenuation directly and remains open.

**Operational consequence (algorithmic SR upper bound).** Whenever (★) holds, the maximum K-change as the noise sweeps `ε ∈ [0, 1]` is

```
|ΔK_max|  =  |1/g_∞ − 1/g_0|,        (capped at the max_runs cutoff M)
```

a per-setup constant determined only by the noise-free and uniform-mix recovery probabilities. This is a *closed-form upper bound on any algorithmic SR effect* for noises in class (★). The sign of `g_0 − g_∞` is per-setup: most Shor-class phase setups have `g_0 > g_∞` (noise increases K), but some have `g_∞ > g_0` (e.g. `a = 254`, `N = 437`, where noise genuinely lowers K) — a real positive-SR cell. Either way the swing is bounded by the closed form and never tunable to an interior optimum (Theorem 6, §3.3.ter). The cliff edge of algorithmic optimization is thereby quantified; the path to noisy-regime advantage lies in hardware-side noise reduction (QEC) or noise-aware post-processing (e.g., Yang–Markidis [arXiv:2605.16074], whose two-stage noise propagation model has the same structure as (★); the identification `ε = 1 − exp(−σ²)` is our verification, not a formula stated in their paper).

Reproduce: `python -m experiments.multi_noise_closed_form` (verifies form (★) for depol, bias_zero, amp_damp at N = 437), `python -m experiments.shor_sigma_curve_model` (verifies phase σ case, R² = 0.95).

### 3.3.ter Theorem 6: No-go for SR-based factoring acceleration

The closed-form swing above promotes to a no-go statement for the specific algorithm-level stochastic-resonance (SR) mechanism, within class (★). Crucially, we do **not** assume noise always hurts: some setups satisfy `g_∞ > g_0` (e.g. `a = 254`, `N = 437`: `ρ = 0.86 > p_0 = 0.69`), where increasing noise genuinely *lowers* the run count — a real positive-SR cell. Theorem 6 bounds the effect regardless of sign.

**Theorem 6 (No-go for SR-based factoring acceleration).** Fix a noise model `M` in class (★) with per-base extraction probability `g(ε) = (1 − ε)·g_0 + ε·g_∞`, `ε ∈ [0, 1]`, `0 < g_0, g_∞ ≤ 1`. Let `E[K(ε)] = E[K_λ^{ideal}] / g(ε)` be the expected run count of the (C)-accumulation algorithm (Theorem 3'). Then, *without any assumption on the sign of `g_0 − g_∞`*:

1. **(No tuned resonance.)** `E[K(ε)]` is monotone on `[0, 1]` (non-decreasing if `g_∞ ≤ g_0`, non-increasing if `g_∞ > g_0`). Hence its optimum is always attained at an endpoint `ε ∈ {0, 1}`; no interior optimum exists. Classical stochastic resonance — an *optimal nonzero* noise level — is therefore impossible: the best operating point is always either zero noise or the maximal-noise endpoint, never a tuned intermediate.
2. **(Closed-form swing.)** The total run-count change over the full noise range is
   ```
   |ΔK| = E[K_λ^{ideal}] · |1/g_∞ − 1/g_0|,
   ```
   fixed by the two endpoint probabilities alone and independent of the σ-profile. For the single-run factoring observable of §3.6 the corresponding bound is `|ΔK_max| = |1/g_∞ − 1/g_0|`, capped at the cutoff `M`. When `g_∞ > g_0` this `|ΔK|` is exactly the maximal *benefit* noise can confer.
3. **(No asymptotic advantage.)** If the reciprocal gap is bounded, `|1/g_∞ − 1/g_0| = O(1)` as `N → ∞`, then `|ΔK| = O(E[K_λ^{ideal}]) = O(log log N)` (Theorem 2) — noise (helpful or harmful) shifts the run count by at most a constant multiplicative factor, never an asymptotic speedup.

*Proof.* (1) `g(ε)` is affine in `ε` with slope `g_∞ − g_0`, hence monotone; `t ↦ 1/t` is monotone on `(0, 1]`, so `1/g(ε)` is monotone, as is `E[K(ε)]` after scaling by the positive constant `E[K_λ^{ideal}]`. A monotone function on a compact interval attains its extrema at the endpoints. (2) Evaluate at `g(0) = g_0`, `g(1) = g_∞` and take absolute value. (3) Substitute the `O(1)` reciprocal-gap hypothesis into (2) and use `E[K_λ^{ideal}] = O(log log N)`. ∎

**Scope and honest caveats.** Theorem 6 rules out the *specific* algorithm-level SR-factoring mechanism for noise in class (★): even in genuine positive-SR cells (`g_∞ > g_0`), the benefit is a closed-form constant `|1/g_∞ − 1/g_0|`, not a tunable resonance, and does not scale with `N`. It does **not** bound (a) noise outside (★) — e.g. amplitude damping, which violates (★) (§3.3.bis); (b) hardware-level error-mitigation or QEC gains; or (c) SR in observables other than run count. The `O(1)` reciprocal gap is an empirical regularity (`|p_0 − ρ| ∈ [0.12, 0.61]` over 12 setups, with no growth toward 1 across `N ∈ {437, 1147, 2491, 4087}`, `shor_n_scaling.py`; see §3.6.bis table), **not** a proven property of Shor-class distributions at cryptographic scale. Theorem 6 is therefore a no-go *within the framework and under its measured regularities*, not an unconditional impossibility result.

Reproduce: `python -m experiments.shor_n_scaling` (endpoint gap `|p_0 − ρ|` vs N), `python -m experiments.multi_noise_closed_form` (class-(★) membership per noise model).

### 3.4 Theorem 4: Compatibility with Regev's multi-base measurement (conditional)

Regev's 2023 algorithm (arXiv:2308.06572) uses `d ≈ √(log N)` bases per quantum circuit and recovers factorization via lattice reduction on `√n + 4` independent measurement vectors. We observe that (C) post-processing applies *coordinate-wise* to Regev's measurements, conditional on a marginal-distribution assumption.

**Assumption (Regev marginal — *partial validity*).** Each measurement of Regev's circuit produces a vector `(k_1, ..., k_d) ∈ {0, ..., Q-1}^d`. For the *marginal* of each coordinate `k_i` to be Shor-like — i.e., `k_i ≈ j_i · Q / r_{a_i}` so that continued fractions of `k_i / Q` recover `r_{a_i}` with positive probability — is sufficient for the (C) framework to apply coordinate-wise.

**Caveat on the assumption.** Public summaries of Regev (2023) indicate the measurements jointly satisfy a *linear constraint* of the form `Σ b_i k_i ≈ 0 (mod r)`, where the `b_i` encode quadratic-character structure. The constraint introduces correlations across coordinates, but does not prohibit Shor-like marginals — it adds a low-dimensional structure that Regev's LLL post-processing exploits for *additional* efficiency beyond what coordinate-wise recovery would provide. Whether each `k_i` is *individually* Shor-like (allowing (C) per coordinate) is a property of the marginal that should be verified directly against Regev's distribution. We have not done so; the empirical numbers below assume independent Shor-like marginals as a working model.

**Theorem 4 (Regev-(C) compatibility, conditional).** Under the Regev marginal assumption, applying (C) post-processing independently to each coordinate of each Regev measurement yields:

- **(a) Per-coordinate Theorem 1.** Whenever `r_{a_i} | L_{before}` for any (run, coordinate) pair, recovery is deterministic regardless of how the marginal distribution of `k_i` is distorted by noise.
- **(b) Run-level K_λ bound.** Let `K_λ^{Regev-(C)}` be the number of Regev runs before `L = λ(N)`. Then `E[K_λ^{Regev-(C)}] ≤ E[K_λ^{ideal}] / d` (under noise-free, ignoring intra-run early termination).
- **(c) Noise invariance.** Corollary 1 of Theorem 1 applies coordinate-wise — any measurement-layer noise on Regev's output preserves the recovery of `r_{a_i}` whenever the corresponding component is covered.

**Proof.** (a) and (c) are direct applications of Theorem 1 and its Corollary 1 to each coordinate independently (which is valid because (C) acts only on the marginal of `k_i`, by the assumption). (b) follows from Theorem 2(b) applied to `K · d` independent samples: `P[L_{K \cdot d} < λ(N)] ≤ ω(λ(N)) · 2^{-K · d}`, so `E[K · d] ≤ log_2 ω(λ(N)) + O(1)`, giving `E[K] ≤ (log_2 ω(λ(N)) + O(1)) / d`. ∎

**Empirical (numpy simulation with assumption).** We simulate Regev runs as `d` independent Shor measurements per run and apply (C) coordinate-wise. Trials = 200, noise-free:

| N    | d | mean K runs | bases ≤ K·d | p99 (runs) | max (runs) |
|------|--:|-----------:|------------:|----------:|----------:|
| 77   | 3 | 2.41       |  7.23       | 11        | 12        |
| 143  | 3 | 3.02       |  9.07       | 11        | 11        |
| 437  | 4 | 1.75       |  7.02       |  6        |  7        |
| 1147 | 4 | 2.50       |  9.98       |  9        |  9        |

For `N = 437`, Regev-(C) requires `~ 1.75` runs in expectation, compared to Regev's `√n + 4 ≈ 4` runs claimed for LLL post-processing. (C) coordinate-wise may be more measurement-efficient under the marginal assumption, but ignores joint correlations that LLL exploits — a trade-off, not a strict improvement.

**Empirical evidence that the marginal assumption is robust.** To test whether the joint constraint actually hurts (C) coordinate-wise, we simulate a *joint-constrained* model: sample `(k_1, ..., k_d)` from independent Shor distributions, then project to the affine subspace `Σ b_i k_i ≡ 0 (mod λ(N))` for random `b_i ∈ Z/N` (a stand-in for Regev's quadratic-character coefficients). We measure `K_λ^{Regev-(C)}` under this model and compare to the independent-marginal model:

| N    | d | indep model | joint model | ratio (joint/indep) |
|------|--:|------------:|------------:|---------------:|
| 77   | 3 | 2.41        | 2.23        | 0.93 |
| 143  | 3 | 3.02        | 3.12        | 1.03 |
| 437  | 4 | 1.75        | 1.66        | 0.95 |
| 1147 | 4 | 2.50        | 2.38        | 0.96 |

Across 4 semiprimes (200 trials each), the joint constraint changes `K_λ^{Regev-(C)}` by at most 7%, and the joint model is *slightly better* in 3 of 4 cases. This is strong empirical evidence that the marginal Shor-likeness is preserved under linear joint constraints — the constraint is "soft" relative to the convergent structure that (C) uses.

**Noise robustness (N=437, d=4, 100 trials).** We additionally measure `K_λ^{Regev-(C)}` on the joint-constrained model under several noise channels and compare to single-base (C) (Phase 2):

| Noise | K_λ^{Regev-(C)} runs | overhead | K_λ^{single-(C)} bases | overhead |
|---|---:|---:|---:|---:|
| noise-free  | 1.59 | 1.00× | 4.97 | 1.00× |
| depol p=0.5 | 2.58 | 1.62× | 9.19 | 1.85× |
| modexp q=0.3| 3.26 | 2.05× | 11.50| 2.31× |
| phase σ=1.0 | 3.36 | 2.11× | 13.06| 2.63× |

**Regev-(C) has *lower* noise overhead than single-base (C) across all three noise channels** — the `d` parallel bases per run amortize per-measurement noise. This is a practical consequence of Theorem 4(c) and suggests Regev-(C) coordinate-wise post-processing may be more hardware-robust than single-base (C), even before accounting for the run-count savings of (b).

**Caveat.** Our joint model uses uniform random `b_i ∈ Z/N` as a stand-in for Regev's quadratic characters; the qualitative result (joint constraint does not destroy marginals, noise overhead is reduced) is robust, but exact efficiency depends on Regev's specific circuit. The noise-invariance of Theorem 1 carries over to *any* setting where per-coordinate (C) succeeds, regardless of joint structure. Reproduce: `python -m experiments.regev_joint` and `python -m experiments.regev_joint --noise 437`.

**Partial head-to-head with RV's filter-then-LLL (N=437, d=4, 50 trials).** We implement a simplified version of RV's Algorithm 6.1 (sympy LLL on the lattice `[S·I_d  S·W; 0  I_|E|]` with `S = 100`) and measure its corruption-filter precision against the same noise model used for our (C) coordinate-wise approach. Corrupted runs are simulated as uniform-random samples (RV's "overwrite" model with corruption probability `p`):

| corruption p | (C) runs to λ(N) | (C) overhead | RV filter precision |
|---:|---:|---:|---:|
| 0.00 | 1.72 | 1.00× | 100.0% |
| 0.10 | 1.86 | 1.08× |  88.5% |
| 0.20 | 2.08 | 1.21× |  78.0% |
| 0.30 | 2.10 | 1.22× |  66.5% |

**Observation.** (C) coordinate-wise degrades gracefully — at 30% corruption rate, it needs only 22% more runs than noise-free. RV's filter precision drops from 100% to 66.5% over the same range (i.e., one-third of "uncorrupted" samples are actually corrupt), at which point RV must rely on Regev's downstream LLL tolerating these stragglers. The (C) framework avoids the filter entirely — per-coordinate verification (`a^d ≡ 1 mod N`) is itself a built-in corruption check.

**Caveat (full comparison requires Regev LLL).** Our RV implementation uses a simplified scaling `S = 100` rather than RV's `S = 2^{An/d}`. We also implement a skeleton of Regev's downstream Algorithm B.1 (`regev_algorithm_b1` in the companion code: builds the lattice `[I_d  ε^{-1} W; 0  I_k]` and LLL-reduces it). The skeleton's output (LLL-reduced basis) does not currently translate to a nontrivial factor of `N` without the full lattice `L` construction over Regev's quadratic characters `b_i` (with `b_i^2 ≡ a_i mod N`) — only with that construction can short vectors give rise to nontrivial square roots and hence factors. The (C) framework bypasses this step entirely: `λ(N)` recovery suffices via Miller-Rabin reduction. A full head-to-head — Regev's filter-then-LLL pipeline (with proper `b_i` and `S = 2^{An/d}`) vs (C) coordinate-wise, both ending in a factor of `N` — is left to future work. Reproduce: `python -m experiments.rv_filter_lll compare`.

### 3.5 Theorem 5: Hybrid (C) + Regev b-trick (factoring)

The empirical comparison in §3.4 reveals an unexpected consequence of Regev's setup. In Regev's framework, the bases used by the quantum circuit are `a_i = b_i² mod N`, where `b_i` are chosen *first* (random in `(Z/N)*`) and known to the algorithm. This is the key structural feature that lets Regev's classical post-processing extract a nontrivial square root of `1` from a short vector `z ∈ L = {z : ∏ a_i^z_i ≡ 1 mod N}`: with `z` in hand, one computes `b = ∏ b_i^z_i mod N`, and `b² ≡ 1 mod N` automatically; if also `b ≢ ±1`, then `gcd(b ± 1, N)` gives a factor.

A purely (C)-based extraction misses this: in Regev's setup, the orders `ord(a_i) = ord(b_i² mod N)` are the *odd part* of `ord(b_i)` — they never contain the factor of 2 that `λ(N)` does for typical semiprimes. Hence pure (C) accumulates `L = lcm(\mathrm{ord}(a_i)) \leq \lambda(N)/2`, and the standard Miller–Rabin reduction (`factor_from_exponent` on `L`) fails because `L` is odd. Pure Regev's b-trick succeeds but requires recovering individual orders `ord(a_i)` independently (no multi-base accumulation gain).

Combining the two: apply (C) coordinate-wise for noise-tolerant per-base order recovery, and on each newly-recovered `ord(a_i)`, immediately compute `b_i^{ord(a_i)} mod N` and check for a nontrivial square root. We call this `(C) + b-trick`.

**Algorithm 2 (Hybrid (C) + Regev b-trick).**

```
For each Regev run:
  For each coordinate (a_i, b_i, k_i):
    Apply (C) post-processing to (a_i, k_i, L)
    If a new r = ord(a_i) is recovered:
      Set L ← lcm(L, r)
      Compute b ← b_i^r mod N
      If b ∉ {1, N-1} and b² ≡ 1 (mod N):
        Return gcd(b ± 1, N) as a factor.
```

The algorithm inherits both: (a) Theorem 1's noise-invariance per coordinate, and (b) Regev's b-trick's ability to factor without requiring `L = λ(N)`. Empirically (N = 437, d = 4, 30 trials each):

| Condition          | (C) lcm only   | Regev b-trick | (C)+b-trick hybrid |
|--------------------|---------------:|--------------:|-------------------:|
| noise-free         | 6.70 (21/30)   | 3.13 (27/30)  | **1.03 (30/30)**   |
| corruption p=0.2   | 6.77 (21/30)   | 3.27 (27/30)  | **1.10 (30/30)**   |
| corruption p=0.3   | 6.87 (21/30)   | 3.50 (27/30)  | **1.27 (30/30)**   |
| depol p=0.3        | 6.73 (21/30)   | 3.37 (27/30)  | **1.10 (30/30)**   |
| phase σ=1.0        | 7.10 (21/30)   | 4.00 (27/30)  | **1.70 (30/30)**   |

Values are mean `K` (runs to factor) and success / 30 trials. The hybrid achieves 100% success in 1–2 runs across all noise/corruption conditions tested.

**Interpretation.** The hybrid is *strictly stronger* than either ingredient alone in this setting: (C) alone fails on the Regev setup (odd `L`), pure b-trick recovers each order independently (no multi-base speedup), and the hybrid combines (C)'s multi-base efficiency with the b-trick's direct factor route. This is a concrete operational benefit of the (C) framework specifically inside Regev's circuit.

**Formal analysis.** We sharpen the empirical statement above with two analytical bounds. Let `v_p := v_2(p - 1)`, `v_q := v_2(q - 1)` (the 2-adic valuations of `p-1`, `q-1`; both `≥ 1` since `p`, `q` are odd). Decompose `(Z/N)* ≅ C_{p-1} × C_{q-1}` via CRT, and for `b ∈ (Z/N)*` let `α_p := v_2(ord(b mod p))`, `α_q := v_2(ord(b mod q))`.

**Lemma 5.1 (per-`b_i` nontrivial sqrt probability).** For uniform random `b ∈ (Z/N)*`:

```
P[b^{ord(b²)} is a nontrivial square root of 1 mod N] = P[α_p ≠ α_q]
                                                      = 1 − 2^{−(v_p + v_q)} · (4^{min(v_p, v_q)} + 2)/3
                                                      ≥ 1/2
```

with equality at `v_p = v_q = 1`, and approaching `2/3` as `min(v_p, v_q) → ∞`.

**Proof.** Let `a = b² mod N`, so `ord(a) = ord(b) / gcd(2, ord(b))`. If `ord(b)` is odd (equivalently `α_p = α_q = 0`), then `ord(a) = ord(b)` and `b^{ord(a)} = 1` — trivial.

Otherwise (at least one of `α_p, α_q ≥ 1`), `ord(a) = ord(b)/2`. Working coordinate-wise in `C_{p-1}` and `C_{q-1}`:

- If `α_p < α_q`: `ord(b mod p) = 2^{α_p} m_p` divides `ord(a) = 2^{α_q - 1} lcm(m_p, m_q)` (since `α_p ≤ α_q - 1` and `m_p | lcm`), so `b^{ord(a)} ≡ 1 mod p`. Meanwhile in `C_{q-1}`, `b^{ord(a)}` has order 2 (since `2 · ord(a) = ord(b)` is a multiple of `ord(b mod q)`, while `ord(a)` itself is not), so `b^{ord(a)} ≡ -1 mod q`. By CRT, `b^{ord(a)} ≡ (1, -1)` — a *nontrivial* square root of 1 mod `N`.
- If `α_p > α_q`: symmetric — `b^{ord(a)} ≡ (-1, 1)`, also nontrivial.
- If `α_p = α_q ≥ 1`: both coordinates yield `-1`, so `b^{ord(a)} ≡ -1 mod N` — trivial.

Hence `P[nontrivial] = P[α_p ≠ α_q]`. The marginals of `α_p, α_q` are computed by standard cyclic-group analysis (in `C_n` with `n = 2^v · m`, the elements with `v_2(ord) ≤ j` form the unique subgroup of order `2^j · m`):

```
P[α_p = 0] = 2^{-v_p}, P[α_p = j] = 2^{j - 1 - v_p} for 1 ≤ j ≤ v_p.
```

`α_p` and `α_q` are independent (CRT). Computing `P[α_p = α_q]` as a sum and using `Σ_{j=1}^{m} 4^j = 4(4^m - 1)/3` gives the closed form. ∎

**Theorem 5 (Hybrid factoring — quantitative bound).** Apply Algorithm 2 with `d` uniform random `b_1, ..., b_d` ∈ (Z/N)*, `a_i = b_i² mod N`. Let `g(η)` denote the per-coordinate (C) recovery probability at `L = 1` under noise model `η` (the `g_M(η)` of §3.3). Let `c := P[nontrivial sqrt]` from Lemma 5.1, so `c ≥ 1/2`.

**(a) Single-run success.** The probability that the hybrid algorithm factors `N` in **one** Regev run is at least:

```
P[1-run success] ≥ 1 − (1 − g(η) · c)^d.
```

**(b) K-run success (with fixed bases).** Letting `X_i` be the indicator that `b_i^{ord(a_i)}` is a nontrivial square root of 1, and conditioning on the fixed `b_i`:

```
P[K-run success | b_1, ..., b_d] = 1 − ∏_{i=1}^{d} (1 − X_i (1 − (1 − g(η))^K)).
```

Taking expectation over `b_i` (independent, each with `E[X_i] ≥ c`):

```
P[K-run success] ≥ 1 − ((1 − c) + c · (1 − g(η))^K)^d.
```

In particular, as `K → ∞` with fixed bases, the success probability approaches `1 − (1 − c)^d`, *bounded by base randomness*: a fraction `≤ (1 − c)^d ≤ 2^{−d}` of base choices have all `X_i = 0` and the algorithm cannot succeed without restarting.

**(c) Restart-augmented expectation.** Modify the algorithm to draw fresh `b_1, ..., b_d` after `K_max` failed runs (a standard amplification trick). Then `E[\text{total runs}] ≤ K_max / (1 − (1 − c)^d) \cdot 1/(1 − \text{conditional fail prob})`. For `K_max = O(\log(1/\epsilon)/g)` and `d ≥ \log_2(1/\epsilon)`, the algorithm succeeds with probability `1 − \epsilon` in `O(K_max)` runs.

**(d) Asymptotic.** For `d = \lceil \log_2(1/\epsilon) \rceil` (so `(1 − c)^d ≤ \epsilon`), and `g(η) = \Omega(1)`, the expected total runs is `O(1)` as `N → ∞`.

At `N = 437`, `d = 4`, with empirical `g ≈ 1` for the small orders typical when `a_i = b_i²`: the 1-run bound gives `P[\text{success}] ≥ 1 − (1/2)^4 = 15/16 ≈ 0.94`, so `E[K] ≤ 16/15 ≈ 1.07`. Empirical mean: `1.03` and success `30/30` in our original §3.5 table; in a second independent run with the same parameters but `max_runs = 10` we observed `27/30` success — matching the `~ 1/16` expected fraction of "all bases bad" trials.

**Scaling to larger N.** We verify the asymptotic prediction `E[K] → 1` by running the hybrid at semiprimes `N ∈ {1147, 2491, 4087}` (same `d = 4`, 30 trials each, `max_runs = 20`, noise-free):

| N    | mean K (hybrid) | success | mean K (b-trick alone) | mean K ((C) lcm alone) |
|------|---:|---:|---:|---:|
| 1147 | **1.50** | 30/30 | 1.60 | 2.57 (29/30) |
| 2491 | **1.10** | 30/30 | 1.10 | 3.17 (27/30) |
| 4087 | **1.43** | 30/30 | 1.53 | 4.33 (27/30) |

The hybrid maintains `E[K] ∈ [1.1, 1.5]` across two orders of magnitude in N, matching the theorem's prediction. For comparison, pure (C) lcm fails to factor in ~3/30 trials (the odd-L issue) and uses 2.6–4.3 runs when it succeeds; pure b-trick succeeds always but needs 1.1–1.6 runs (no multi-base accumulation). The hybrid combines the strengths.

**Noise robustness at large N (N = 4087, d = 4, 20 trials).** We further confirm that the hybrid's noise tolerance carries over to the largest N tested:

| Noise condition | hybrid mean K | success |
|---|---:|---:|
| noise-free      | 1.55 | 20/20 |
| depol p=0.3     | 2.15 | 20/20 |
| phase σ=1.0     | 2.70 | 20/20 |

At `N = 4087`, hybrid still factors in `≤ 2.7` runs even under heavy `phase σ = 1.0` noise, with 100% success. This is consistent with Theorem 5's noise-invariance claim and the destructive-noise scaling of Theorem 3.

**Proof of (a) and (b).** For coordinate `i`, the hybrid succeeds in `K` runs iff `X_i = 1` AND `(C)` recovers `ord(a_i)` in at least one of the `K` runs. The second event has probability `1 − (1 − g(η))^K` (each run is an independent Bernoulli `g(η)`). The two are independent (b_i and k_i are independent). So `P[\text{coord } i \text{ succeeds in } K | b_i] = X_i (1 − (1 − g(η))^K)`, giving the conditional formula. Marginalizing over the d independent b_i with `E[X_i] = c`:

```
P[fail in K | b₁, …, b_d] = ∏ᵢ (1 − Xᵢ(1 − (1 − g)^K))
E[above]                  = ((1 − c) + c(1 − g)^K)^d
```

(by independence of the d coordinates in our model). Hence (b). Setting `K = 1` and using `(1 − c) + c(1 − g) = 1 − cg` recovers (a). ∎

Reproduce: `python -m experiments.rv_filter_lll factor`.

### 3.5.bis Empirical robustness: the hybrid is ε-dominant (no noise-adaptive selector needed)

A natural follow-on question is whether the *best* post-processing method changes with the noise level — i.e. whether there is a noise regime in which standalone (C)-lcm or standalone b-trick beats the hybrid, which would justify a **noise-adaptive selector** that switches method as a function of measured ε. We tested this directly with an ε×method bake-off (`experiments/method_bakeoff.py`): at each σ all three methods consume the *same* measurement sequence per trial, so the comparison is fair, and we report the always-hybrid *regret* = `E[K]_hybrid − min_method E[K]`.

| σ | ε = 1−e^(−σ²) | (C) lcm | Regev b-trick | hybrid | always-hybrid regret |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.000 | 2.55 | 1.30 | 1.30 | +0.00 |
| 0.20 | 0.039 | 2.75 | 1.60 | 1.50 | +0.00 |
| 0.40 | 0.148 | 2.85 | 1.65 | 1.50 | +0.00 |
| 0.60 | 0.302 | 3.05 | 1.70 | 1.55 | +0.00 |
| 0.80 | 0.473 | 4.60 | 2.25 | 2.10 | +0.00 |
| 1.00 | 0.632 | 5.05 | 2.90 | 2.75 | +0.00 |
| 1.50 | 0.895 | 7.45 | 5.55 | 5.25 | +0.00 |
| 2.00 | 0.982 | 8.80 | 4.90 | 4.60 | +0.00 |

(N = 1147, d = 4, 20 trials, max_runs = 20; success ≥ 80% gating; the σ = 0 tie 1.30 = 1.30 resolves to hybrid.) Across the full measured range `ε ∈ [0, 0.98]` the hybrid is the (weakly) best method at **every** noise level: max regret is `+0.00`, well inside one standard error. There is no crossover, so a noise-keyed selector would have **zero synergy** to capture.

We record this honestly as an *empirical observation at these parameters*, not a theorem: the hybrid is ε-dominant over the tested grid, so adaptive method-switching on ε is unnecessary. This is distinct from (i) the Ragavan–Vaikuntanathan fixed filter-then-LLL pipeline (a single non-adaptive recovery rule), (ii) the Yang–Markidis recoverability predictor (which predicts *whether* a run is usable, not *which* method to run), and (iii) windowed phase estimation (arXiv:2509.05010, a circuit-level change). It is also the constructive counterpart to Theorem 6: not only is there no tunable interior noise optimum, there is likewise no noise regime in which a different post-processing method should be selected. Reproduce: `python -m experiments.method_bakeoff --N 1147 --trials 20`.

### 3.6 Trial-level noise sensitivity at K-bin boundaries (universal mechanism observation)

We document a mechanistically clean phenomenon in the hybrid (C) + Regev b-trick algorithm at `(N, d) = (437, 4)`: under phase noise of magnitude `σ ∈ [0.005, 0.100]`, every sampled base set exhibits a small set of trials moving across one of the K-bin boundaries of single-run factoring — most often the `K = 1 / K = 2` boundary, occasionally the `K = 2 / K = 3` boundary. Across 13 independent seeds × 200 trials × 12 σ values (= 31,200 trial-measurements), the *boundary-flip mechanism is universal* (13/13 seeds), but the *direction* (success ↔ failure) is base-set-dependent and shows no statistically significant net bias: mean SR `= +0.144%`, sd `1.016%`, SE `0.282%`, `t = 0.51`, `p (1-sided) = 0.31` at σ = 0.050.

The mechanism is qualitatively a member of the stochastic resonance family (Benzi et al. 1981; Wellens–Buchleitner 2004 for the quantum variant) and the σ-curve quantitatively matches the classical Benzi shape (saturation plateau followed by overload decline). It provides a clean conceptual bridge between integer factoring and the broader noise-as-resource literature (ENAQT, Plenio–Huelga 2008) but is too small in absolute terms (1–7 trial flips per 200 = `±0.3–2%` per seed) to have cryptographic implications.

**13-seed σ-scan at `N = 437`, `d = 4`.** We measured `K` with 200 trials per σ across 12 σ values and 13 seeds:

| σ | mean SR % (over seeds) | sd | SE | t | p (1-sided) |
|---:|---:|---:|---:|---:|---:|
| 0.005 | +0.211% | 1.052 | 0.292 | +0.72 | 0.236 |
| 0.010-0.025 | +0.190% | 1.054 | 0.292 | +0.65 | 0.258 |
| 0.035-0.075 | +0.165% | 1.041 | 0.289 | +0.57 | 0.284 |
| 0.100 | +0.142% | 1.041 | 0.289 | +0.49 | 0.312 |
| 0.150 | +0.001% | 1.000 | 0.277 | +0.00 | 0.500 |
| 0.200 | −0.318% | 0.974 | 0.270 | −1.18 | 0.881 |

The plateau structure of the σ-curve is preserved: within each seed, σ ∈ [0.005, 0.100] gives identical K-mean (deterministic flip set). The decline at σ ≥ 0.150 is real but the net direction is small. At σ = 0.200 the mean turns slightly negative — see *direction asymmetry* below.

**Per-seed mechanism observation.** For each seed we computed the K-histogram at σ = 0 and σ = 0.050 and identified the dominant K-bin flip:

| seed | K_base | SR at σ=0.050 | direction | dominant flip | mag |
|---|---|---|---|---|---|
| 1 | 2.200 | +0.682% | helps | K=2→K=1 | 1 |
| 2 | 1.555 | +1.929% | helps | K=2→K=1 | 1 |
| 3 | 1.720 | −0.872% | hurts | K=1→K=2 | 1 |
| 4 | 1.720 | +0.581% | helps | K=2→K=1 | 2 |
| 5 | 1.630 | +0.613% | helps | K=3→K=1 (long jump) | 1 |
| 6 | 1.550 | +0.323% | helps | K=3→K=2 (secondary) | 1 |
| 7 | 1.515 | +0.330% | helps | K=2→K=1 | 2 |
| 8 | 2.315 | +0.432% | helps | K=2→K=1 | 2 |
| 9 | 2.215 | +1.580% | helps | K=2→K=1 | 5 |
| 10 | 2.090 | −0.957% | hurts | K=2→K=3 (secondary, neg) | 2 |
| 11 | 1.820 | −1.099% | hurts | K=1→K=2 | 2 |
| 12 | 1.815 | −1.377% | hurts | K=1→K=2 | 5 |
| 13 | 1.720 | −0.291% | hurts | K=1→K=2 | 3 |

**Boundary distribution**: K = 1 / K = 2 in 10/13 seeds (76.9%); K = 2 / K = 3 in 2/13 (15.4%); K = 3 ↔ K = 1 long-jump in 1/13 (7.7%). All 13 seeds show boundary-flip activity: *mechanism universality is empirically robust*. The K-mean change for each seed exactly matches the histogram-derived total (`Δ K_total / n_trials = mean SR / 100 × K_base`, integer precision, 13/13).

**Direction independence from K_baseline.** Seeds 3 and 4 share `K_baseline = 1.720` but show opposite directions (−0.872% vs +0.581%); seeds 4 and 13 also share this `K_baseline` but give +0.581% vs −0.291%. Direction is determined by the trial-level K-distribution of the specific base set drawn, not by the K-mean aggregate.

**σ-curve direction asymmetry.** Positive-SR seeds and negative-SR seeds follow qualitatively different σ-curves:

- *Positive seeds*: saturation plateau (σ ∈ [0.005, 0.100]), then decline back toward baseline as σ ≥ 0.150 (some flipped trials are pushed back). The effect is bounded above by the small K = 2 (or K = 3) population available to flip into K = 1.
- *Negative seeds*: monotonic worsening — K-mean increases with σ throughout the range, no decline. The K = 1 population (180+ trials of 200) is virtually unbounded relative to K = 2, so noise has a continuous reservoir of trials to push out of K = 1.

This asymmetry — a direct consequence of the heavily K = 1-skewed noise-free K-distribution at this cell — explains why the cross-seed mean SR at high σ approaches zero or slightly negative: positive seeds saturate while negative seeds amplify.

**Speculative mechanism.** Small phase noise broadens the QFT measurement peak around `j · Q / r` slightly. For each trial the broadening either exposes an alternative convergent candidate that succeeds (K decrease) or perturbs a previously-found one into failure (K increase); the outcome per base set is determined by which side of the active boundary holds more *near-miss* trials in the noise-free measurement. The plateau structure follows from the discreteness of the affected trial set: once σ exceeds the per-trial peak width that the borderline trials require to flip, additional σ has no marginal effect within the plateau. The overload decline at large σ reflects the same dynamic in reverse — flipped trials become unstable and revert.

**Why this effect is small compared to existing noise-as-resource literature.** Existing examples in the broader noise-as-resource paradigm — environment-assisted quantum transport in the FMO complex (Plenio–Huelga 2008, with efficiency boost from ~13% to ~47%), classical stochastic resonance in neural detection (Benzi et al. 1981; with detection-rate improvement from ~0% to ~50%), and bistable-potential SR (signal amplification 30%+) — all share a common structural feature: the *noise-free* baseline is *sub-functional* (stuck in a local minimum, trapped by destructive interference, or sub-threshold), and noise *restores function* by breaking the trap. Effect magnitudes are large because the system goes from non-functional to functional. Our system has no such structural feature: the noise-free hybrid algorithm already factors in ~K_baseline runs (mean K = 1.92 at (N, d) = (437, 4), meaning ~95% of trials succeed in 1–2 runs), and there is no destructive interference, bistability, or threshold trapping the noise-free baseline. Phase noise can therefore only flip a small number of *borderline* trials at K-bin boundaries (1–7 per 200 = ±0.3–2% per seed), not switch the algorithm from "broken" to "working." The mechanism we identify is the *same* as in classical SR (saturation plateau + overload decline), but the system is in the regime of "marginal improvement on a working system" rather than "function-restoration on a stuck system."

**Engineered amplification — borderline-population is the bottleneck (★ direct evidence).** We tested two thinned variants of the hybrid to probe the structural prediction that the SR magnitude is bounded by the *borderline-trial population* at the active boundary, not by the algorithm's overall success rate.

- *Over-thinned* (smallest convergent denominator only + no (C) augmentation): K_baseline = 19.87 (vs 2.08 in the full hybrid with 100 trials, 9.55× sub-functional). At this baseline almost all trials reach max_runs = 20, leaving *no borderline population* to flip. The measured SR is 0.00% at σ ∈ {0.050, 0.150} across 3 seeds — *sub-functional alone does not amplify the effect*.

- *Mild thinned* (all convergent denominators, (C) augmentation removed): K_baseline = 2.92 (1.4× sub-functional, a moderate increase). The per-seed |SR| at σ = 0.050 amplifies to 4.03–4.44% (vs 0–1.16% in the full hybrid, ~5× larger). The cross-seed mean shifts from +0.06% (full) to −1.43% (mild), with 1 positive / 2 negative direction split in the mild case. The mechanism is unchanged (boundary-flip + occasional high-K rescue), but the (C) augmentation evidently acts as a noise *buffer*: its removal exposes the underlying mechanism in a more amplified form at the per-seed magnitude level, without itself biasing the net direction.

These two thinning experiments — null result with over-thinning, positive amplification with mild thinning — jointly confirm the structural prediction. Engineered amplification of *net direction* (rather than per-seed magnitude) would require a separate mechanism (e.g., direction-correlated noise or base-set selection); we leave this to future work. The mild-thinned variant is intentionally worse than the full hybrid and has no practical advantage; it demonstrates that the same boundary-flip mechanism operates in the sub-functional regime, just with a larger borderline population available to noise.

**Algorithm-structure regime map for noise-as-resource (★ now empirically validated).** The amplification results suggested a regime map for noise-as-resource susceptibility in multi-base quantum factoring. We now empirically validate all five entries through direct measurement (5/5 measured):

| Algorithm structure | Per-seed \|SR\| | Mean SR | K_base | Source |
|---|---|---|---|---|
| Single-base Shor (1994) | 0–1.10% | −0.04% (3+/2−) | 10.38 | **measured** (5 seeds × 100 trials × 3 σ) |
| Multi-base Regev (LLL post-processing) | 1.27–3.95% | −0.31% (2+/3−) | 2.44 | **measured** (5 seeds × 50 trials × 3 σ, faithful LLL implementation) |
| Hybrid (C) + b-trick — full (this work) | 0–1.93% at (437, 4) | +0.14% (8+/5−) | 1.82 | **measured** (13 seeds × 200 trials × 12 σ) |
| Hybrid mild-thinned (no (C) augmentation) | 4.03–4.44% | −1.43% (1+/2−) | 2.92 | **measured** (3 seeds × 100 trials × 3 σ) |
| Hybrid over-thinned (smallest convergent only) | 0% | 0% | 19.87 | **measured** (3 seeds × 100 trials × 3 σ) |

**★ Universal direction stochasticity across algorithm structures.** Our 5-algorithm regime map measurement reveals an unexpected unifying pattern: cross-seed direction is *base-set-stochastic* (not algorithm-determined) in *every* tested algorithm structure. The mean SR magnitudes are all small (|mean SR| ≤ 1.5%) and none reach statistical significance at our sample sizes. The regime map's original directional predictions (Shor "small", Regev "negative") are *qualitatively* confirmed (Shor mean ≈ 0, Regev mean weakly negative −0.31%) but the magnitudes are smaller than the naive "LLL is fragile" picture suggests. The per-seed magnitude distribution *does* depend on algorithm structure: Pure Shor is most narrowly distributed (max |SR| = 1.10%), Pure Regev has a wider distribution (max = 3.95%), and our Hybrid at (1147, 2) reaches max +9.44% via high-K rescue. The dominant first-order observation is therefore *direction stochasticity is universal, magnitude distribution is algorithm-dependent* — not the "single algorithm is best" picture our regime map originally suggested.

**Pure Shor σ-scan verification.** Implementing single-base Shor (d = 1, no (C) accumulation, no fast path; b-trick factoring retained) at N = 437 with 5 seeds × 100 trials × 3 σ gives K_baseline = 10.38, dominated by ~50% trials reaching max_runs = 20 (where the chosen `b` does not satisfy the b-trick nontrivial-sqrt condition). The borderline-trial population is small (~15% of trials at K = 2–7), and noise can flip only this small population. Mean SR = −0.04% at σ = 0.050 (sd 0.85%, SE 0.38%, t = −0.11, p ≈ 0.46) confirms the prediction of small effect. Per-seed |SR| reaches a maximum of 1.10%. Sign test: 3 of 5 positive direction.

**Pure Regev σ-scan verification (faithful LLL implementation).** Implementing the b-trick factoring through LLL on a Regev-style lattice yields K_baseline = 2.44, matching the c ≥ 1/2 lower bound from Lemma 5.1. Our implementation uses self-contained LLL reduction (δ = 0.75), Regev-style (d+1)-dimensional embedding lattice with k_vec measurements, multi-S scaling (S ∈ {Q, Q/2, 2Q}), enumerated short vector search (basis rows + small linear combinations with coefficients ∈ [−2, 2]), and multi-measurement accumulation over K runs. Mean SR at σ = 0.050 is −0.31% (sd 2.91%, SE 1.30%, t = −0.24, p ≈ 0.60); per-seed |SR| range is 1.27–3.95%. Direction is 2 positive / 3 negative (weak negative lean as predicted). The result *qualitatively* confirms the "LLL is fragile" prediction in direction (weakly negative mean) but the *magnitude* is smaller than the naive picture suggested. Our LLL implementation does not use BKZ reduction or the exact lattice from Regev (2023/JACM 2025) Algorithm B.1; a full Regev implementation might show different per-seed magnitudes, but the qualitative direction-stochasticity finding should persist.

This regime-map perspective and the universal-stochasticity finding together suggest a refined understanding: *noise-as-resource in quantum factoring is a mechanism-level phenomenon (boundary flips) that operates across all multi-base algorithm structures with stochastic direction; the differences between algorithms are in per-seed magnitude distribution rather than systematic direction bias*. The hybrid (C) framework is unique in providing a *buffer* that absorbs borderline trials and thereby reduces per-seed magnitude — its removal (mild thinning) restores larger per-seed magnitudes. Pure Shor shows the smallest magnitudes (narrow K-distribution + dead-trial dominance). Pure Regev with LLL shows intermediate magnitudes with a *weak* negative direction lean (consistent with LLL fragility but smaller than naively predicted).

**Cross-cell verification of the mechanism.** The boundary-flip mechanism makes three regime-level predictions, which we test against existing measurements at other `(N, d)` cells:

| (N, d) | K_baseline | regime | observed SR % | match |
|---|---|---|---|---|
| (437, 8) | 1.19 | ceiling | +0.00% | ✓ exact |
| (1147, 8) | 1.06 | ceiling | +0.47% | ✓ small |
| (2491, 4) | 1.07 | ceiling | −1.88% | ✓ small |
| (2491, 8) | 1.00 | floor | −1.50% | ✓ noise |
| (1147, 1) | 5.78 | noise floor | −0.53% (multi-seed 3×100tr, sd 4.28%) | ✓ variance > effect |
| (437, 4) | ~1.9 (V3: 1.92, 13-seed mean: 1.82) | active boundary | +0.144% (13×200tr, p=0.31) | ✓ small, no net direction |
| (1147, 2) | 2.43–2.92 | active boundary | +1.24% single → +0.42% (4×1000tr) → +3.35% (5 seeds × 100tr, p≈0.12) | ✓ multi-seed reveals high-K rescue, regresses partially |
| (2491, 2) | 2.30 | active boundary | -4.89% single → +0.30% (5 seeds × 100tr, p=0.45) | ✓ single-seed outlier regresses to small mean |

The ceiling cells show near-zero SR (no `K ≥ 2` population to flip); the noise-floor cell shows variance dominating effect; the active-boundary cells show small magnitudes that regress under multi-seed measurement. The single-seed grid measurements at `(437, 1, 2, 3), (1147, 3, 4), (2491, 1, 2)` (mostly positive, one −4.89% at 2491, 2) are consistent with this picture: per-seed variance ≈ effect magnitude, mean direction not yet determined.

**Mechanism diversification at higher `K_baseline` (★ new observation).** A 5-seed σ-scan at `(N, d) = (1147, 2)` with `K_baseline ≈ 2.92` (range 2.24–3.39, 100 trials × 5 σ values per seed) reveals a *richer* mechanism than the K = 1 / K = 2 boundary flip observed at `(437, 4)`. While 2 of 5 seeds (seeds 1 and 5) show the same K = 1 / K = 2 boundary mechanism familiar from `(437, 4)` (±1–3% SR), the remaining 3 seeds exhibit a qualitatively different pattern we call *high-K rescue*: trials sitting at very high K bins (`K = 8, 11, 15`, and even `K = 20`, i.e., max_runs failures) move to *moderate* K bins (`K = 4, 5`) under noise. In seed 3, for instance, one trial moves from K = 15 to K = 5, one from K = 11 to K = 5, and one from K = 20 to K = 6, producing a net K reduction of 32 across 100 trials (SR = +9.44%). In seed 4, three trials move K = 8 → K = 4 (a 12-K-unit shift), producing SR = +8.56%. The cross-seed mean SR at σ = 0.050 is +3.35% (sd 5.37, SE 2.40, t = 1.39, p ≈ 0.12 with t-distribution (df = 4), or 0.082 under normal approximation), borderline at the 5-seed sample but dominated by the high-K-rescue seeds.

This high-K rescue mechanism is *absent* at `(437, 4)` simply because the K-distribution there has very few high-K trials to rescue (`K_baseline = 1.92` implies almost all trials succeed in K = 1 or K = 2). At cells with broader K-distributions (`K_baseline ≈ 3+`), high-K trials become available targets for noise-induced rescue, providing a *new* mechanism channel that amplifies per-seed |SR| substantially. The boundary distribution at `(1147, 2)` is correspondingly more diverse: only 40% of seeds show K = 1 / K = 2 as the dominant transition (vs 77% at `(437, 4)`), with K = 8 → K = 4, K = 15 → K = 5, K = 3 → K = 1 long-jump, and K = 1 → K = 2 each represented in different seeds.

Note that the *net cross-seed mean* at `(1147, 2)` remains marginal (+3.35%, p ≈ 0.12 with t-distribution (df = 4), or 0.082 under normal approximation) because the high-K-rescue seeds (3, 4) are partially offset by a K = 1 → K = 2 negative-direction seed (5, SR = −2.68%). The mechanism *magnitude* amplifies with K_baseline; the *net direction* does not.

**What we earlier over-claimed and have since retracted.**

1. *`+17.86%` peak at `(1147, 2, 0.01)`* (150 trials × 1 seed): re-measured at 1000 trials × 4 seeds gave mean `+0.42%` (`t = 0.50`, not significant). The original is a single-seed direction fluctuation, within mechanism variance.
2. *Polynomial scaling `SR ∝ N^α`*: rejected once `N = 2491` cells showed near-zero or negative SR. The boundary-flip mechanism predicts no N-scaling; effect depends on `K_baseline` and seed-specific direction.
3. *`σ_opt ∝ N^α`* ("small lock, small wiggle"): rejected by σ-scan finding `σ_opt ≈ 0.010` independent of N — consistent with the discrete saturation plateau (any σ > threshold flips the same trials).
4. *"Anti-Optimization Principle" — `d = 1` universal positive SR*: motivated by single-seed AOP measurements; multi-seed re-measurement at `(1147, 1)` gave `−0.53% ± 4.28%`, consistent with the noise-floor regime.
5. *V3 sign test `p = 0.03` as significance*: V3 showed all 5 σ values below baseline at a single seed; we treated this as `(1/2)^5 ≈ 3%`. This is invalid because σ values within the saturation plateau are perfectly correlated (same boundary trials flip in all). The proper significance test is between-seed direction, which our 13-seed measurement finds *not* significant (`p = 0.31`).
6. *"Goldilocks K_baseline ≈ 2 regime as a robust SR cell"*: refined. `K_baseline ≈ 2` *does* mark the regime where boundary-flip activity is detectable, but the *direction* is base-set-stochastic, not systematically positive. The boundary flip mechanism is the robust observation; the SR-as-positive-effect interpretation was a single-seed artifact.

**Caveats.**

- *Effect magnitude is tiny.* 1–7 trial flips per 200, ±0.3–2% K change per seed, no cryptographic implication.
- *No net SR direction.* Mean SR across 13 seeds (`+0.144%`, p = 0.31) is not statistically distinguishable from zero. Direction is base-set-stochastic.
- *Direction asymmetry of σ-curve* reflects K-distribution skew, not a feature of phase noise itself.
- *Phase-noise specific.* Depolarizing and amplitude damping show monotone degradation, not boundary flips.
- *Mechanism not formalized at base-set level.* We predict that some K-boundary is active and the σ-curve shape, but not the specific boundary location or flip direction (both depend on base-set-specific K-distribution structure).

**Open questions.**

1. *Direction bias at larger samples.* 30+ seeds may reveal a small systematic bias hidden in our 13-seed noise. The current point estimate (+0.14%) is consistent with both a true zero and a small positive effect.
2. *Active-boundary determinants.* Which K-boundary becomes active at a given seed is not predictable from `K_baseline` alone; finer-grained analysis of base-set composition is needed.
3. *Universality across N.* The mechanism likely extends to other active-boundary cells; concrete verification candidates: `(1147, 2), (1147, 3), (4087, 4)`, each at ~2000-trial multi-seed scale.
4. *Long-jump events* (e.g., K = 3 → K = 1 in seed 5) may reflect cascades or true noise-induced trajectory divergence.

We document this as a universal trial-level mechanism observation with stochastic per-seed direction: the multi-base hybrid Shor algorithm has K-bin boundaries that are sensitive to small phase noise, the mechanism is robust across all sampled base sets, but the cross-seed direction is not systematically biased toward improvement. The qualitative connection to the noise-as-resource paradigm is clear; the quantitative net effect is null at our sample size.

Reproduce: `python -m experiments.sigma_scan_437` (12 σ × 3 seeds × 200 trials baseline scan), `python -m experiments.sigma_scan_437_extend 4 13` (additional 10 seeds), `python -m experiments.analyze_histograms` (boundary-flip analysis).

### 3.6.bis Self-correction (v0.3.0, 2026-06-14): closed-form replaces boundary-flip framework

**The mechanism described above as "boundary flip" is, on closer analysis, the finite-trial expression of a single smooth closed form**:

```
p(σ) = ρ + (p_0 - ρ) · exp(-σ²)
E[K(σ)] = (1 - (1-p(σ))^M) / p(σ)
```

derived directly from the noise-averaged FFT under per-amplitude phase noise:
`E[|FFT(a·e^{iε})_k|²] = (1-e^{-σ²})/Q + e^{-σ²}·P_0(k)`.

**This form has been verified to fit five algorithm classes**:

| Algorithm | R² | Note |
|---|---|---|
| Grover | +0.88 | k-iter accumulation |
| Shor (pure + b-trick) | +0.95 | gating |
| QPE isolated (no b-trick) | +0.96 | clean QFT |
| Simon | +0.99 | Hadamard + XOR |
| Hybrid (C)+b-trick (this paper's §3.6 setup) | +0.91 | direct internal fit |

**Reinterpretation of §3.6 observations under the closed form**:

- *"Boundary flip" mechanism*: not a separate mechanism — it is the K-binning of a smoothly shifting p(σ) under finite-trial sampling.
- *"Plateau"* (σ ∈ [0.005, 0.100], mean SR ≈ +0.15%): not structural; `α σ² · (p_0 - ρ)` simply falls below the finite-trial SE in that range.
- *"Universal direction stochasticity"*: `sign(p_0 - ρ)` is determined per (a, b) setup. We previously called it "stochastic" because we did not derive that sign; now we can.
- *"Plateau + overload" two-regime picture*: a single exponential decay rendered in two qualitative bands by the finite-trial resolution.

**What is retained**:
- The raw measurement data (13 seeds × 12 σ × 200 trials = 31,200 measurements).
- The five-cell regime map predictions (5/5 measured, consistent with closed form).
- The qualitative bridge to noise-as-resource literature.
- The conclusion that SR-based factoring acceleration is precluded.

**What is retracted**:
- The "boundary flip" lexicon as a distinct mechanism.
- The "deterministic flip set within plateau" reading (it is statistical, not deterministic).
- The "stochastic direction" framing as unexplainable — `sign(p_0 - ρ)` explains it.

**Bound on the SR effect (closed-form consequence)**:

`|ΔK_max| = |1/ρ - 1/p_0|` (cap at max_runs M).

For phase noise on Shor-class algorithms, this bound rules out asymptotic SR-based factoring speedup: as N grows, the gap `|p_0 - ρ|` does not produce O(1/log N) speedups. We measure `|Δ| = |p_0 − ρ|` over 12 setups (3 per N, 300–1,000 MC samples each, `shor_n_scaling.py`):

| N | mean p_0 | mean ρ | mean \|Δ\| | min \|Δ\| | max \|Δ\| |
|---:|---:|---:|---:|---:|---:|
| 437 | 0.741 | 0.410 | 0.446 | 0.173 | 0.607 |
| 1147 | 0.443 | 0.402 | 0.211 | 0.118 | 0.260 |
| 2491 | 0.486 | 0.287 | 0.372 | 0.260 | 0.483 |
| 4087 | 0.476 | 0.390 | 0.259 | 0.220 | 0.297 |

Across a 9× range in `N` the gap stays bounded in `[0.12, 0.61]` with no monotone growth toward 1 — the SR swing neither diverges nor vanishes with problem size. Notably the smallest-order setup at each `N` (`r_a ∈ {3, 4}`) has `p_0 < ρ` (a positive-SR cell, `g_∞ > g_0`) while larger-order setups have `p_0 > ρ`, confirming the per-setup `sign(p_0 − ρ)` direction. This is formalized as **Theorem 6 (§3.3.ter)**, which also shows the swing is never tunable to an interior optimum regardless of the sign of `p_0 − ρ`.

**Relation to Yang-Markidis (arXiv:2605.16074, ICS Workshops '26)**: their Eq. (3), §5 states a two-stage noise propagation model `(1-ε)·[p_s∗K_{σ_0}] + ε·Σ_h ν_h·[p_h∗K_{σ_h}]` in which the mixing weight `ε` is a *conceptual* parameter — in their words "the total weight transferred out of the intended family" — left unspecified, neither fitted nor given as a function of σ (full text verified, no appendix; 2026-06-14). Our closed form has the same structure with `ε = 1 − exp(−σ²)`. We do **not** claim this as a new mechanism: the `exp(−σ²)` decay of off-diagonal coherences under i.i.d. Gaussian phase noise is the standard dephasing result (Nielsen–Chuang §8.3). Relative to their work our contribution is narrow and verificational: (i) we *verify numerically* that their qualitative weight coincides with the standard-dephasing factor (R² = 0.95 at N = 437); (ii) we *measure the boundary* of this identification — it holds for phase/depolarizing/bias (R² 0.95–0.996) but breaks for amplitude damping (R² 0.03, structural, §3.3.bis); and (iii) we confirm the same form across five algorithm classes. They, in turn, contribute hardware-level recoverability features from 680 IBM runs. The two are complementary: empirical hardware characterization (theirs) and analytical verification + boundary mapping (ours) — we do not position this as filling an analytical gap they were unable to close.

Reproduce: `python -m experiments.shor_sigma_curve_model` (Shor pure, R²=0.95), `python -m experiments.grover_sigma_curve_model` (Grover, R²=0.88), `python -m experiments.qpe_isolated_sigma` (QPE, R²=0.96), `python -m experiments.simon_sigma_curve` (Simon, R²=0.99), `python -m experiments.hybrid_sigma_curve` (this §3.6 setup, R²=0.91). See `sr_sigma_curve_model.md` for the unified framework.

### 3.7 Joint interpretation

Theorem 1 says: *once* `r_a | L`, success is deterministic regardless of noise. Theorem 2 says: an ideal algorithm reaches this state in `O(log log log N)` bases in expectation. Theorem 3 says: under destructive noise, the actual algorithm reaches it in `E[K_λ^{ideal}] / g_M(η)` bases — i.e., overhead exactly `1/g_M(η)`. **Theorem 3' (§3.3.bis, v0.3.0) extends Theorem 3 to all noise models that admit a coherence-loss decomposition `g_M = (1−ε)·g_0 + ε·g_∞`, including phase noise via `ε = 1 − exp(−σ²)`. Amplitude damping is shown to be outside this class (structural).** Theorem 4 says: under a marginal-distribution assumption, the (C) framework applies coordinate-wise to Regev's multi-base measurements with corresponding reduction in run count. Theorem 5 closes the loop: combining (C)'s `λ(N)` recovery with Regev's `b`-trick yields a complete factoring algorithm that empirically requires ~5× fewer runs than standalone Regev at `N ∈ {437, 1147, 2491, 4087}`. Together they explain the algorithm's empirical robustness: a noise-adjusted logarithmic number of measurements suffice to enter a regime immune to further measurement noise, and the framework composes naturally with Regev's multi-base circuit (modulo marginal assumption).

The §3.6 trial-level "boundary-flip" observation, originally cast as a distinct mechanism, **is now subsumed by Theorem 3' (§3.3.bis) plus §3.6.bis self-correction**: it is the finite-trial K-binning of a smoothly shifting `p(σ) = ρ + (p_0 − ρ)·exp(−σ²)` whose direction is set per-(a, b)-setup by `sign(p_0 − ρ)`. The cross-cell regime-map predictions remain valid (5/5 measured), and they now have a closed-form derivation rather than an empirical taxonomy. Verified across five algorithm classes (Grover R² = 0.88; Shor pure R² = 0.95; QPE isolated R² = 0.96; Simon R² = 0.99; Hybrid (C)+b-trick R² = 0.91) and three noise models (phase R² = 0.95; depolarizing R² = 0.9953; bias zero R² = 0.9963), the closed form is the unifying analytical object behind the empirical observations of §3.3–§3.6.

**Operational consequence.** The bound `|ΔK_max| ≤ |1/g_∞ − 1/g_0|` (§3.3.bis) precludes algorithmic SR-based factoring acceleration under any noise in class (★) — phase, depolarizing, bias zero, modexp. This is a *quantitative cliff edge*: the path to noisy-regime quantum advantage lies in hardware-level noise reduction (QEC) or noise-aware post-processing (Yang–Markidis arXiv:2605.16074), not in algorithm-level SR exploitation. Algorithm and hardware research directions are thereby cleanly separated by an analytical boundary.

## 4. Empirical verification

We verify Theorem 1 across diverse noise models. For each measurement, we log `(r_a, L_before, condition, success)` where `condition := r_a | L_before` and `success := (C output = r_a)`. We count:

- `violations` := count of `condition = True ∧ success = False` (theorem predicts 0)
- `lucky` := `condition = False ∧ success = True` (recovery via convergents alone)
- `covered` := `condition = True`

Noise models tested (each implemented in `noise.py`):

| Category | Model | Parameter |
|---|---|---|
| A (measurement only) | depolarizing | `p` ∈ [0, 1] |
| A | readout bit flip | `p` per bit |
| A | bias-to-zero (adversarial) | `p` |
| B (QFT input) | phase decoherence | `σ` (rad) on amp |
| B | amplitude damping (T1) | `γ` decay rate |
| C (structural) | modular exponentiation error | `q` per `f(x)` value |

**Result.** Across `N ∈ {77, 143, 209}` × 11 noise setups × 500 trials = 16,500 measurements: **violations = 0**. The theorem holds in every case.

The `lucky` count is small (1-4 per setup), confirming most successes flow through the `covered` region. The `missed` count (`condition = False ∧ success = False`) is dominated by early trials (`L = 1`) under high noise.

See `demo.py --verify <N>` to reproduce. We extend the same verification protocol to larger semiprimes `N ∈ {1147, 2491, 4087}` with a representative 4-noise subset and 100 trials per setup — see Appendix C. Per-noise-model breakdown for `N = 77` is in Appendix B.

## 5. Related work

We outline how Theorem 1 relates to existing literature, including which ingredients are folklore and where the gap (if any) lies.

**Knill (1995, LANL tech report) — "On Shor's quantum factor finding algorithm: Increasing the probability of success and tradeoffs involving the Fourier transform modulus."** Proposes running the quantum step *twice for the same base* and taking the lcm of denominators recovered from convergents. This is *single-base lcm*. Our (C) extends this to multi-base lcm and augments the candidate pool with `divisors(L)`.

**McAnally (2001) — "A Refinement of Shor's Algorithm" (arXiv quant-ph/0112055).** Same-base accumulation: enlarges the Fourier transform modulus to `Q ≈ 2wn³` and enumerates *all* convergent denominators in a single measurement, achieving near-certainty per run. Same-base structure, different mechanism (Q enlargement vs L accumulation across bases).

**Bourdon-Williams (2007, arXiv quant-ph/0607148).** Sharp asymptotic lower bound on per-measurement success probability (~94%) via Carmichael-based analysis. Theoretical bound; no algorithmic change.

**Ekerå (2021, arXiv 2007.10044) — "On completely factoring any integer efficiently in a single run of an order-finding algorithm."** Shows the order of a *single* base typically suffices for complete factorization via smoothness extension. Different direction from our multi-base accumulation; complementary.

**A quantum algorithm for computing the Carmichael function (2021, arXiv 2111.02488).** Algorithm 1 of this paper essentially runs the multi-base lcm to recover `λ(N)`. Closest prior work to our framework. The paper *does not augment subsequent measurements' post-processing with `divisors(L)`*; each order is computed independently.

**Ekerå (2024, arXiv 2201.07791) — "On the success probability of quantum order finding."** Comprehensive survey of post-processing techniques. Surveyed approaches include continued fractions, lattice post-processing, offset search, divisor search of recovered `r/d`, Seifert's joint solving. *Multi-base accumulation across different `a` values is not surveyed as a post-processing technique here*.

**Bach-Shallit (1996) — *Algorithmic Number Theory*, Vol 1.** Foundational textbook coverage of `r | λ(N)`, multi-element lcm convergence to `λ(N)`, and factoring from a known exponent (Miller-Rabin reduction).

**Pomerance et al. (2017, arXiv 1707.07193) — *The expected number of elements to generate a finite group*.** Shows `e(G) ≤ d + 2.752...` where `d` is the maximum Sylow generator count. For semiprime `(Z/NZ)*`, this gives constant-bounded expected number of bases to *generate* the group (hence reach `λ(N)` via lcm).

### Three approaches to noise-tolerant Regev factoring

The Regev 2023 algorithm (Regev, JACM 2025) factors n-bit integers via `~√n` runs of a multi-base circuit followed by lattice post-processing. Regev's own analysis requires all runs to be uncorrupted. Three independent approaches to handling corrupted runs have appeared:

**Ragavan-Vaikuntanathan (2023/2025, arXiv 2310.00899) — *Space-Efficient and Noise-Robust Quantum Factoring*.** Extends Regev's lattice post-processing with a *filter-then-LLL* stage: corrupted samples are detected via short-vector search on a constructed lattice and discarded; standard Regev LLL post-processing then runs on the uncorrupted survivors. Requires a "well-spread" noise distribution assumption. Operates entirely within Regev's lattice framework.

**Ekerå-Gärtner (2024) — concurrent with the above.** Shows that under a different (stronger) assumption on the corruption distribution, *standard* Regev post-processing already tolerates a constant fraction of corrupted samples without filtering. Lattice framework preserved; analysis-only contribution.

**This work (Theorem 4, conditional).** Replaces Regev's lattice post-processing with `(C)` *coordinate-wise* — each measurement coordinate `k_i` is processed independently via standard Shor-style continued-fraction recovery plus the divisor search of the accumulated exponent `L`. Output is `λ(N)`, from which factorization follows by standard methods (Miller–Rabin reduction on a known exponent). Requires the *marginal Shor-likeness* of each `k_i` (verified empirically against a joint-constrained model in §3.4). Noise-tolerance inherits Theorem 1's noise invariance per coordinate without needing a filter or distribution assumption.

The three approaches are *orthogonal*: RV adds a filter on top of LLL; EG24 modifies the assumption while keeping the algorithm; our work replaces the post-processing entirely. Direct head-to-head empirical comparison (RV's Algorithm 6.1 vs (C) coordinate-wise on Regev measurements, under matched noise) is a natural next step but is not performed here.

**Position of our contribution.** Every individual ingredient — `r | λ(N)`, multi-base lcm, divisor search given exponent, continued fractions, factor_from_exponent — is folklore or established in the references above. The *clean separation theorem* (noise invariance via the classical verification gate, formalized as Theorem 1) and its corollaries appear to not be stated explicitly; nor does the *application of multi-base lcm post-processing to Regev's coordinate-wise measurements* (Theorem 4) appear in the noise-tolerant-Regev literature surveyed above. We propose Theorems 1–4 as: (i) the folklore packaging that gives a uniform explanation for the empirical robustness of multi-base post-processing (Theorems 1–3), and (ii) an alternative post-processing approach for Regev's framework, distinct from the lattice-based methods of RV and EG24 (Theorem 4).

## 6. Limitations and discussion

**No measurement saving from adaptive base selection.** A natural attempt — pre-filter bases by `a^L ≢ 1 mod N` so that every measured base extends `L` — yields no improvement in measurement count over the random-selection + fast-path scheme (`shor_quantum_adaptive` vs `shor_quantum_multi`; see §9 of README). The fast path already skips classical-only bases; pre-filtering shifts the work but does not eliminate any measurement.

**Marginal value of lattice post-processing.** The Knill-Mosca lattice approach strengthens single-base, multi-measurement convergent recovery. In our multi-base framework, `lcm(r_a₁, ..., r_aₖ)` already accomplishes equivalent recovery across bases. We did not find a regime where additional lattice basis reduction (LLL or simpler) provides measurable benefit over the lcm strategy.

**The theorem is conditional on `r_a | L`.** When `L = 1` (start of algorithm), the condition fails and `(C) ≡ (B)`. Theorem 2 (§3.2) precisely characterizes this transient: `E[K_λ] ≤ 1 + Σ_{ℓ | λ(N)} 1/(ℓ^{s_ℓ} - 1)`, which is `O(log log log N)` for typical semiprimes. Empirically, `K_λ ≤ 14` across all 17,000 trials at `N ≤ 4087` (Appendix D).

**Scope of the framework — group-exponent learning.** The multi-base accumulation underlying Theorems 1–3 is specific to *group-exponent learning* — Shor's order finding is fundamentally an algorithm that learns `λ(N)` (the exponent of `(Z/N)*`) via lcm of orders, with factorization following as a byproduct. The structure does *not* naturally extend to:

- **Discrete logarithm** (single-instance: a specific `x` for a specific `(g, h)`; no natural multi-base lcm).
- **Abelian hidden subgroup problem** (single-instance: a specific `H` for a specific `f`; sublattice accumulation is well-studied but does not gain noise-invariance from a (C)-style divisor verifier).

It does naturally extend, with no algorithmic change, to *multi-prime* factoring (`N = p_1 · ... · p_k` with k ≥ 2): Theorems 1–3 hold verbatim, with `λ(N) = lcm(p_1 - 1, ..., p_k - 1)` and `s_ℓ ∈ {1, ..., k}` indicating how many factors `p_i - 1` achieve the maximum `v_ℓ`.

A useful *composition* is `(C) + Pohlig-Hellman`: use (C) to learn `λ(N)` and hence `ord(g)` for any `g`, then apply Pohlig-Hellman for discrete log on `g`. This is a hybrid algorithm, not an extension of the (C) framework itself.

**Classical verification gate is assumed exact.** Theorem 1 relies on `pow(a, d, N)` returning the correct value. On real hardware, this could in principle be corrupted by classical computation errors, but at standard machine precision this is negligible compared to quantum noise.

**Scale of verification.** Our experiments span `N ∈ {77, 143, 209}` with 500 trials per noise setup (§4) and `N ∈ {1147, 2491, 4087}` with 100 trials per setup (Appendix C), totaling 17,700 measurements. The largest `N = 4087 = 61 · 67` requires a 24-qubit counting register (`Q = 2^24 ≈ 16.8M` amplitudes, 256 MB at complex128) — at the upper end of single-machine state-vector simulation. Scaling to cryptographic `N` requires no algorithmic change; the remaining bottleneck is pure classical simulation memory. Hardware verification on existing NISQ devices is a natural next step.

## 7. Code, reproducibility

All experiments are reproducible via the companion repository:

```
pip install numpy
git clone https://github.com/Hashevolution/shor
cd shor
python demo.py --verify 77 143 209          # reproduces §4 table (Theorem 1)
python verify_large_run.py 1147 2491 4087   # reproduces Appendix C (large N)
python -m experiments.k_lambda_dist         # reproduces Appendix D (Theorem 2)
python -m experiments.g_eta 437             # measures g_M(η) per noise model (Theorem 3)
python -m experiments.k_lambda_alg 437      # reproduces §3.3 table (Theorem 3)
python -m experiments.regev_c               # reproduces §3.4 table (Theorem 4)
python -m experiments.rv_filter_lll factor  # reproduces §3.5 table (Theorem 5 hybrid)
python -m experiments.sr_aop                # reproduces §3.6 single-seed (N, d) AOP grid
python -m experiments.sigma_scan_437        # reproduces §3.6 baseline σ-scan (3 seeds)
python -m experiments.sigma_scan_437_extend # reproduces §3.6 extended seeds (4-13) + backfill
python -m experiments.analyze_histograms    # reproduces §3.6 per-seed flip identification
python -m experiments.sigma_scan_general 1147 2 5 100 compact   # cross-cell verification
python -m experiments.hardware_calibrated   # reproduces Appendix E (hardware proxy)
python demo.py --noise3 77 143              # reproduces 3-noise comparison
python demo.py --adaptive 77 143 209 323    # reproduces §9.1 negative result
```

The codebase is `~1000` lines of numpy. No quantum libraries required.

Raw measurement data for §3.6 is stored in `experiments/sigma_scan_437_d4_results.txt` (initial 3 seeds), `experiments/sigma_scan_437_d4_extended.txt` (seeds 4-13), and `experiments/sigma_scan_437_d4_histograms.txt` (per-seed K-distributions at σ ∈ {0, 0.050} for boundary-flip analysis). All scripts are resumable (immediate per-cell save, skip-existing on re-run).

## 8. Conclusion

We formalize a folklore observation about classical post-processing of multi-base Shor order-finding measurements: augmenting continued-fraction recovery with a divisor search over the accumulated lcm `L` makes the procedure deterministic whenever the order of the current base divides `L` — independent of the measurement distribution, hence robust to all measurement-layer noise models. We prove two complementary quantitative bounds: Theorem 2 on the number of bases `K_λ` needed to reach the noise-invariant regime (refining a corollary of Pomerance et al. 2017), and Theorem 3 characterizing the exact noise overhead `1/g_M(η)` for destructive noise models. We verify Theorem 1 across 17,700 trials at six composite sizes (zero violations), Theorem 2 across 17 semiprimes × 1,000 trials, and Theorem 3 across 9 noise setups on N=437. Theorem 4 extends the framework coordinate-wise to Regev's multi-base measurements (conditional on a verified marginal-distribution assumption); Theorem 5 packages this into a hybrid factoring algorithm that empirically requires ~5× fewer runs than standalone Regev at N ∈ {437, 1147, 2491, 4087}. Two natural attempts to further reduce measurement count (adaptive base selection, lattice post-processing) yield no improvement, indicating Theorem 1 captures the essential structure at this scale.

Additionally (§3.6), at `(N, d) = (437, 4)` we measured `K` with 13 independent base sets × 200 trials × 12 σ values (31,200 trial-measurements total) and identified a universal trial-level mechanism: phase noise σ ∈ [0.005, 0.100] flips a small number of trials at one of the K-bin boundaries (K=1/K=2 in 77% of seeds, K=2/K=3 in 15%) of single-run factoring. The mechanism is robust (13/13 seeds), and the σ-curve follows the classical Benzi–Buchleitner shape (saturation plateau + overload decline). The direction (success ↔ failure) is determined by base-set-specific K-distribution near the active boundary and shows no statistically significant net bias across seeds (mean SR = +0.144%, t = 0.51, p = 0.31). The observation provides a clean conceptual bridge between integer factoring and the broader *noise-as-resource* literature (ENAQT; Plenio–Huelga 2008) but is too small in absolute terms to enable cryptographic advantage.

The work is best viewed as expository and a foundation for explorations of more substantive improvements: noise-resilient implementations on hardware, the `lucky` region under structured measurement distributions, the empirical comparison of (C) coordinate-wise vs RV's filter-then-LLL on matched corrupted Regev data, and a deeper theoretical account of when noise transitions from barrier to weak resource in quantum factoring.

## 9. Companion direction (v0.5.0+): magic across the speedup ladder

The order-finding results above (Theorems 1–6) concern the *classical post-processing* of measurements; the question of *why* the underlying quantum algorithm is fast is orthogonal. From v0.5.0 we add a companion direction that quantifies the **nonstabilizerness ("magic")** an algorithm spends, measured by the stabilizer 2-Rényi entropy $M_2(|\psi\rangle)=-\log_2\big(2^{-n}\sum_{P\in\mathcal P_n}\langle\psi|P|\psi\rangle^4\big)$, computed exactly in $O(n\,4^n)$ via an XOR–Walsh–Hadamard identity (`magic.py`, cross-checked against brute force to $10^{-15}$; 42-assertion regression suite `magic_proofs_check.py`). Full statements/proofs are in `magic-results.md`; the standalone write-up is `magic-paper-draft.md`; prior art in `magic-prior-art.md`.

**Thesis.** The *amount/density* of magic discriminates the *type* of speedup, and magic lives not on the visible circuit surface but in the *nonlinearity of the problem*, which oracles and FFT shortcuts can hide.

| algorithm | speedup | $M_2$ behavior | density $M_2/n$ |
|---|---|---|---|
| Simon / BV | exponential (query) | $0$ (affine oracle = Clifford) | $0$ |
| Grover | quadratic | peak $\to 3$ bits, $0$ at answer | $\to 0$ |
| Shor (comb) | exponential | $\propto t$ | $\sim0.4$–$0.55$ |
| Shor (in-circuit) | exponential | $\to L$ (max) | $\to 1$ |

**Results (new in this work).** (i) **Lemma 1 / Cor. 1–2:** a flat (uniform-amplitude, uniform-phase) state on support $S$ has $M_2=0$ iff $S$ is affine; hence the graph state $2^{-n/2}\sum_x|x\rangle|f(x)\rangle$ is stabilizer iff $f$ is $\mathbb F_2$-affine — Simon keeps $M_2\equiv0$ (magic not needed for *query* advantage) while Shor's multiplicative oracle forces $M_2>0$. (ii) **Prop. 2–3:** single-marked Grover has a closed form whose peak saturates at **exactly 3 bits** ($a^2=\tfrac12$), density $M_2/n\to0$. (iii) **Prop. 2$'$, 5, 5$'$ (coding-theory of marker sets):** the magic of a flat state over $W\subseteq\mathbb F_2^n$ is controlled by the autocorrelation of $\mathbf 1_W$, *not* by minimum Hamming distance ($\{0,1,2,3\}$ and $\{0,1,2,4\}$ share $d_{\min}=1$ but $M_2=0$ vs $1.54$); the exact zero-test is $M_2=0\iff A_W\in\{0,M\}$; generic (Sidon) sets give $M_2=\log_2\frac{M^3}{7M-6}$, and a uniform random $M$-subset has exact $\mathbb E[\xi]\,M^4=(7M^2-6M)+\tfrac{7(M)_4}{N-3}+\tfrac{N(N-1)(N-2)(N-4)(M)_8}{(N)_8}$. (iv) **Prop. 6 (oracle-hiding = $T$-cost):** $M_2(|\psi_f\rangle)>0\iff f$ has a degree-$\ge2$ ANF monomial $\iff U_f$ requires non-Clifford ($T$) gates — the magic a query oracle hides equals the $T$-cost of its gate decomposition (Simon $0$ vs Shor modexp $T_{\rm proxy}=4\to156$, $M_2=1.54\to5.23$).

**Honestly credited to prior work.** Shor's magic↔period quantitative law is Paviglianiti et al. (arXiv:2605.05347). The flat-state SRE **closed form** (Prop. 4, additive energy $M_2=-\log_2(M^{-4}\sum_x E(W\cap(W{\oplus}x)))$) is the uniform-support case of Tarabunga–Castelnovo's Rokhsar–Kivelson/SMF formula (*Quantum* **8**, 1347 (2024), Eq. 8); the $2\log_2 M$ growth rate is the saturation of the standard bound $M_\alpha\le2\log_2 R$. Nearest neighbors differentiated: hypergraph-state magic (phase-encoded / RM(2): arXiv:2308.01886, 2602.23687) vs our support-encoded / RM(1); weight-enumerator SRE tools (Quantum Lego, arXiv:2308.05152); permutation-invariant/Dicke magic (arXiv:2402.08551, symmetric support only). All four code↔magic literatures were full-text compared (see `magic-prior-art.md`).

**Asymptotic tightness of Prop. 5$'$ (D3 finding).** Prop. 5$'$ closes $\mathbb E[\xi]$ exactly, but the Jensen step $\mathbb E[M_2]\ge -\log_2\mathbb E[\xi]$ leaves the gap $J(M,N):=\mathbb E[M_2]+\log_2\mathbb E[\xi]\ge 0$ uncontrolled. An automated rediscovery loop (`experiments/discover_d3_jensen.py`; design in `ai-discovery-engine-design.md`) — which first passes a sanity gate by re-deriving the Sidon constant and the additive-energy closed form from scratch — finds $J\propto 1/N$ in the sparse regime $M^2\ll N$ (measured slope $-0.96$ across $M\in\{6,8,10\}$, $n\in\{10,\dots,13\}$). Consequence: **$-\log_2\mathbb E[\xi]$ is not merely a Jensen lower bound but an asymptotically tight estimate of $\mathbb E[M_2]$, with absolute error $O(1/N)$.** The loop also records the saturation boundary $M^2/N\to 1$ where the $1/N$ law breaks (additive-collision saturation), and reports $\kappa(M):=J\cdot N$ as having no simple closed form in the dictionary $\{M,M^2,M(M-1),M^3,M\log_2 M,1\}$ — kept as an open sub-question rather than over-fit.

**Practical bearing (honest).** The most applied of these is Prop. 6: it converts an oracle function's nonlinearity into a fault-tolerant resource signal, since $T$-count dominates surface-code FTQC cost (one distilled magic state per $T$). Under the synthesis model "degree-$d$ ($\ge2$) ANF monomial $\to(2d{-}3)$ Toffolis $\to 7T$" (`oracle_ftqc_estimate.py`), $T_{\rm est}=0\iff M_2=0\iff f$ affine, and $T_{\rm est}$ tracks $M_2$. This per-output-bit synthesis is an *upper bound* (real modexp is far cheaper via windowed arithmetic — Gidney–Ekerå 2021), so its usable value is the exact affine zero-test, the qualitative magic$\leftrightarrow T$ law, and a quick upper-bound indicator for unstructured oracles. The marker-set results give cheap magic prediction/zero-tests in place of $O(n4^n)$ SRE evaluation; the Grover ladder is foundational/diagnostic (and its 3-bit figure is a *state* quantity, not a full-circuit $T$-count).

## References

1. Shor, P. W. (1994). *Algorithms for quantum computation: discrete logarithms and factoring*. FOCS.
2. Knill, E. (1995). *On Shor's quantum factor finding algorithm: Increasing the probability of success and tradeoffs involving the Fourier transform modulus.* Los Alamos National Laboratory tech report.
3. McAnally, D. (2001). *A Refinement of Shor's Algorithm*. arXiv:quant-ph/0112055.
4. Nielsen, M. A. & Chuang, I. L. (2000). *Quantum Computation and Quantum Information*. Cambridge.
5. Bourdon, P. S. & Williams, H. T. (2007). *Sharp probability estimates for Shor's order-finding algorithm*. arXiv:quant-ph/0607148.
6. Ekerå, M. (2021). *On completely factoring any integer efficiently in a single run of an order-finding algorithm*. arXiv:2007.10044.
7. Ekerå, M. (2024). *On the success probability of quantum order finding*. arXiv:2201.07791. ACM Transactions on Quantum Computing.
8. Bach, E. & Shallit, J. (1996). *Algorithmic Number Theory*, Vol. 1. MIT Press.
9. Erdős, P., Pomerance, C., & Schmutz, E. (1991). *Carmichael's lambda function*. Acta Arithmetica 58(4), 363-385.
10. Pomerance, C. et al. (2017). *The expected number of elements to generate a finite group with d-generated Sylow subgroups*. arXiv:1707.07193.
11. Regev, O. (2023, JACM 2025). *An efficient quantum factoring algorithm*. arXiv:2308.06572.
13. Ragavan, S. & Vaikuntanathan, V. (2023). *Space-Efficient and Noise-Robust Quantum Factoring*. arXiv:2310.00899 / eprint 2023/1559.
14. Ekerå, M. & Gärtner, J. (2024). *Extending Regev's factoring algorithm to compute discrete logarithms.* (See discussion in [13].)
12. *A quantum algorithm for computing the Carmichael function* (2021). arXiv:2111.02488.

**Companion direction (§9) — magic.**
15. Leone, L., Oliviero, S. F. E., & Hamma, A. (2022). *Stabilizer Rényi entropy*. Phys. Rev. Lett. **128**, 050402.
16. Paviglianiti, A., Seclì, M., Tirrito, E., & Savona, V. (2026). *The true cost of factoring: linking magic and number-theoretic complexity in Shor's algorithm*. arXiv:2605.05347.
17. Tarabunga, P. S. & Castelnovo, C. (2024). *Magic in generalized Rokhsar–Kivelson wavefunctions*. Quantum **8**, 1347.
18. Chen, J., Yan, Y., & Zhou, Y. (2024). *Magic of quantum hypergraph states*. Quantum **8**, 1351.
19. Kagamihara, D. & Tsuchiya, S. (2026). *Stabilizer Rényi entropy of 3-uniform hypergraph states*. arXiv:2602.23687.
20. Cao, C., Gullans, M. J., Lackey, B., & Wang, Z. (2024). *Quantum Lego Expansion Pack: Enumerators from Tensor Networks*. arXiv:2308.05152.
21. Passarelli, G., Fazio, R., & Lucignano, P. (2024). *Nonstabilizerness of permutationally invariant systems*. arXiv:2402.08551.
22. *A fast and exact approach for stabilizer Rényi entropy via the XOR–FWHT algorithm* (2026). arXiv:2512.24685.

---

## Appendix A. Reference implementation

```python
# multi_base.py (selected)

def divisors(n):
    """Divisors of n, ascending."""
    small, large = [], []
    i = 1
    while i*i <= n:
        if n % i == 0:
            small.append(i)
            if i != n // i: large.append(n // i)
        i += 1
    return small + large[::-1]

def convergent_denominators(k, Q, max_denom):
    """All convergent denominators of k/Q up to max_denom."""
    a, b = k, Q
    h_prev, h_curr = 0, 1
    k_prev, k_curr = 1, 0
    denoms = []
    while b != 0:
        q = a // b
        a, b = b, a - q*b
        h_prev, h_curr = h_curr, q*h_curr + h_prev
        k_prev, k_curr = k_curr, q*k_curr + k_prev
        if 0 < k_curr <= max_denom:
            denoms.append(k_curr)
        elif k_curr > max_denom:
            break
    return denoms

def minimize_order(a, N, candidate):
    if candidate <= 0 or pow(a, candidate, N) != 1: return 0
    r = candidate
    for p in prime_factors(r):
        while r % p == 0 and pow(a, r // p, N) == 1:
            r //= p
    return r

def C(a, N, k, Q, L):
    """The (C) post-processing procedure."""
    cands = set(convergent_denominators(k, Q, N - 1))
    if L > 1: cands.update(divisors(L))
    valid = [d for d in cands if d > 0 and pow(a, d, N) == 1]
    if not valid: return 0
    return minimize_order(a, N, min(valid))
```

## Appendix B. Per-noise-model violation counts (N = 77)

| Noise | covered/500 | violations | lucky | missed |
|---|---:|---:|---:|---:|
| noise-free | 493 | 0 | 3 | 4 |
| depolarizing p=0.3 | 499 | 0 | 1 | 0 |
| depolarizing p=0.8 | 494 | 0 | 3 | 3 |
| readout flip p=0.3 | 487 | 0 | 3 | 10 |
| bias zero p=0.5 | 499 | 0 | 1 | 0 |
| phase σ=1.0 | 496 | 0 | 3 | 1 |
| phase σ=2.5 | 481 | 0 | 2 | 17 |
| amp damp γ=0.01 | 491 | 0 | 2 | 7 |
| amp damp γ=0.05 | 492 | 0 | 3 | 5 |
| modexp q=0.3 | 493 | 0 | 3 | 4 |
| modexp q=0.8 | 476 | 0 | 3 | 21 |

Identical pattern (violations = 0) for `N = 143` and `N = 209`; see `demo.py --verify` output.

## Appendix C. Extended verification at larger N

We additionally verify Theorem 1 at three larger semiprimes, with a representative subset of noise models (one each from categories A/B/C plus noise-free) and reduced trial count (100 per setup) due to increased per-measurement cost (≈ 0.45 s at `N = 1147`, ≈ 1.9 s at `N = 4087`).

| N            | Noise         | covered/100 | violations | lucky | missed |
|--------------|---------------|------------:|-----------:|------:|-------:|
| 1147 = 31·37 | noise-free    |          91 |          0 |     2 |      7 |
| 1147         | depol p=0.8   |          93 |          0 |     3 |      4 |
| 1147         | phase σ=2.5   |           9 |          0 |     4 |     87 |
| 1147         | modexp q=0.8  |          61 |          0 |     2 |     37 |
| 2491 = 47·53 | noise-free    |          97 |          0 |     2 |      1 |
| 2491         | depol p=0.8   |          96 |          0 |     1 |      3 |
| 2491         | phase σ=2.5   |           1 |          0 |     1 |     98 |
| 2491         | modexp q=0.8  |          59 |          0 |     2 |     39 |
| 4087 = 61·67 | noise-free    |          88 |          0 |     3 |      9 |
| 4087         | depol p=0.8   |          93 |          0 |     2 |      5 |
| 4087         | phase σ=2.5   |          50 |          0 |     3 |     47 |
| 4087         | modexp q=0.8  |           1 |          0 |     1 |     98 |

**Total: 1,200 additional measurements, 0 violations.** Combined with §4 the verification covers 17,700 measurements.

Note that under extreme noise the `covered` fraction can drop sharply — `phase σ = 2.5` at `N = 2491` and `modexp q = 0.8` at `N = 4087` both reach `covered = 1` — because the measurement is so degraded that `L` cannot accumulate. The theorem then applies only vacuously over most trials, but on the few trials where it does apply, success remains 100%. The non-monotonicity across `N` for the same noise (`phase σ = 2.5` covered: 9, 1, 50 at `N = 1147, 2491, 4087`) is a sample-size artifact at 100 trials and reflects which random seeds happened to land on bases with small `r_a` early.

Reproduce: `python verify_large_run.py 1147 2491 4087` (≈ 26 min total on a single CPU).

## Appendix D. Empirical verification of Theorem 2

We measure the empirical distribution of `K_λ` across 17 semiprimes, with 1,000 independent trials each (uniform random base selection, classical order recovery, `L` accumulation until `L = λ(N)`). We compare:

- `mean`, `p99`, `max` — empirical statistics of K_λ
- `thm2(c)` := `1 + Σ_{ℓ | λ(N)} 1/(ℓ^{s_ℓ} - 1)` — the sharp expectation bound
- `thm2(b)` := `log₂ ω(λ(N)) + 2` — the simpler asymptotic bound

| N    | p   | q   | λ(N) | ω(λ) | mean | p99 | max | thm2(c) | thm2(b) |
|------|----:|----:|-----:|-----:|-----:|----:|----:|--------:|--------:|
| 15   |  3  |  5  |    4 |    1 | 1.70 |   6 |   9 |   2.000 |   2.000 |
| 21   |  3  |  7  |    6 |    2 | 1.58 |   4 |   5 |   1.833 |   3.000 |
| 33   |  3  | 11  |   10 |    2 | 1.47 |   4 |   5 |   1.583 |   3.000 |
| 35   |  5  |  7  |   12 |    2 | 2.21 |   7 |  10 |   2.500 |   3.000 |
| 77   |  7  | 11  |   30 |    3 | 1.82 |   5 |   7 |   2.083 |   3.585 |
| 91   |  7  | 13  |   12 |    2 | 1.97 |   7 |   9 |   2.125 |   3.000 |
| 143  | 11  | 13  |   60 |    3 | 2.37 |   7 |  10 |   2.750 |   3.585 |
| 187  | 11  | 17  |   80 |    2 | 2.13 |   7 |  13 |   2.250 |   3.000 |
| 209  | 11  | 19  |   90 |    3 | 1.93 |   5 |   7 |   2.083 |   3.585 |
| 221  | 13  | 17  |   48 |    2 | 2.33 |   8 |  10 |   2.500 |   3.000 |
| 247  | 13  | 19  |   36 |    2 | 2.23 |   7 |   9 |   2.500 |   3.000 |
| 323  | 17  | 19  |  144 |    2 | 2.22 |   7 |   9 |   2.500 |   3.000 |
| 391  | 17  | 23  |  176 |    2 | 2.02 |   8 |  11 |   2.100 |   3.000 |
| 437  | 19  | 23  |  198 |    3 | 1.83 |   5 |   7 |   1.933 |   3.585 |
| 1147 | 31  | 37  |  180 |    3 | 2.42 |   8 |  12 |   2.750 |   3.585 |
| 2491 | 47  | 53  | 1196 |    3 | 2.11 |   7 |  14 |   2.129 |   3.585 |
| 4087 | 61  | 67  |  660 |    4 | 2.26 |   8 |  11 |   2.475 |   4.000 |

In all 17 rows, `mean ≤ thm2(c)`. The bound is consistently tight, within 0.5 in every case. The maximum observed `K_λ` over all 17,000 trials is 14 (occurring at N=2491). Reproduce: `python -m experiments.k_lambda_dist`.

## Appendix E. Hardware-calibrated noise simulation

We additionally simulate (C) at noise levels calibrated to published IBM Quantum Eagle 127-qubit processor specifications (median values, 2024–2025):

- T1 ≈ 150 μs (relaxation), T2 ≈ 100 μs (dephasing)
- 1-qubit gate error ≈ 0.03 %, 2-qubit (ECR) gate error ≈ 1 %
- Readout error ≈ 2 % per bit

For an N=15 Shor circuit (12 qubits, ~50 μs depth, ~30 gates with ~10 two-qubit), these map to combined noise parameters: `readout_flip = 0.02, modexp_error = 0.05, amplitude_damp = 0.002, phase_sigma = 0.3, depolarizing = 0.01` (all applied simultaneously).

**Result.** With all five noise channels active simultaneously:

| Metric | Value |
|---|---|
| trials (Theorem 1) | 500 |
| covered | 499 |
| **violations** | **0** |
| success rate | 100% |
| mean K_λ^alg (Theorem 3) | 3.84 |
| baseline K_λ^alg (noise-free) | 3.04 |
| overhead ratio | 1.26× |

The 1.26× overhead is consistent with Theorem 3's destructive-class prediction `1/g_M(η)`: empirically `g_M ≈ 0.79 · g_0` under the mixed hardware noise. (C) recovers the order `r_a` deterministically once `L = λ(N) = 4` is reached, which happens in ~ 4 base draws.

**Caveat.** This is a numpy-based simulation that maps published hardware noise specifications to our five-channel noise model, not an actual hardware execution. A real hardware demonstration (qiskit + IBM Quantum account) is left to future work; the present simulation is a faithful proxy under the noise-channel mapping. Reproduce: `python -m experiments.hardware_calibrated`.
