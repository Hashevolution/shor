# A Noise-Invariant Determinism Theorem for Multi-Base Post-Processing in Shor's Order Finding

*Draft, 2026-06-11.*
*Companion code: [github.com/Hashevolution/shor](https://github.com/Hashevolution/shor).*

## Abstract

We give three theorems and one conditional compatibility observation about classical post-processing of multiple-base order-finding measurements in Shor's algorithm. (1) **Noise-invariant determinism**: maintaining an accumulated exponent candidate `L` — the least common multiple of orders recovered from previous bases — and augmenting standard continued-fraction post-processing with a divisor search over `L` yields a procedure (which we call (C)) that recovers the order `r_a` of any new base `a` deterministically whenever `r_a | L`, *independent of the measurement distribution*. (2) **Logarithmic coverage time (ideal)**: for a semiprime `N = pq`, the expected number of independent uniform bases `K_λ` required for `L = λ(N)` satisfies `E[K_λ] ≤ 1 + Σ_{ℓ | λ(N)} 1/(ℓ^{s_ℓ} - 1)`, where `s_ℓ ∈ {1, 2}` records the ℓ-Sylow overlap of `(Z/p)*` and `(Z/q)*`. (3) **Noise scaling**: under a class of "destructive" noise models (depolarizing, bias, modexp), the actual algorithm K_λ scales exactly as `E[K_λ^{alg}(η)] = E[K_λ^{ideal}] / g_M(η)` where `g_M(η)` is the per-base extraction probability. Together: a noise-adjusted logarithmic number of measurements suffice to enter the noise-immune regime of (1). We verify (1) across 17,700 measurements (zero violations), (2) across 17,000 trials (mean within bound), and (3) across 9 noise setups on N=437 (mean error ~11% for the destructive class). Theorem 4 demonstrates conditional compatibility with Regev's 2023 multi-base framework (numpy simulation, 200 trials × 4 N), provided each measurement coordinate's marginal is Shor-like. We discuss the theorem's relation to prior work — Knill's lcm trick (1995), McAnally's larger-Q convergent enumeration (2001), the Bach-Shallit textbook treatment of `r | λ(N)`, and the quantum algorithm for computing the Carmichael function (2021) — and conclude that while every individual ingredient is folklore, the *clean statement of the noise-invariance corollary* appears not to be made explicit in surveyed literature. We additionally show that two natural attempts to improve measurement count beyond this scheme — adaptive base selection and lattice-based joint post-processing — yield no further reduction, suggesting the theorem captures the essential structure of the problem at this scale.

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

**Out of scope: structural noise.** For `phase_sigma` and `amplitude_damp` (which are not destructive in the above sense), the per-base recovery probability `g_M(s, η)` depends strongly on the state `s` (specifically, recovery for larger `r_a` is harder under peak smearing). The analogous prediction `E[K_λ^{ideal}] / g_M(L=1, η)` is a *lower* bound on `E[K_λ^{alg}]`, but underestimates the truth by a factor 1.4–2.4 at moderate noise and >5 at extreme noise. Deriving a closed-form `g_M(s, η)` for structural noises is left as future work.

### 3.4 Theorem 4: Compatibility with Regev's multi-base measurement (conditional)

Regev's 2023 algorithm (arXiv:2308.06572) uses `d ≈ √(log N)` bases per quantum circuit and recovers factorization via lattice reduction on `√n + 4` independent measurement vectors. We observe that (C) post-processing applies *coordinate-wise* to Regev's measurements, conditional on a marginal-distribution assumption.

**Assumption (Regev marginal).** Each measurement of Regev's circuit produces a vector `(k_1, ..., k_d) ∈ {0, ..., Q-1}^d`, where the *marginal* distribution of each coordinate `k_i` satisfies `k_i ≈ j_i · Q / r_{a_i}` for some integer `j_i` — i.e., is statistically identical to Shor's single-base measurement distribution for base `a_i`.

(This assumption is consistent with public summaries of Regev's algorithm but ignores potential joint correlations between coordinates that Regev's lattice post-processing may exploit.)

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

**Caveat.** If the Regev marginal assumption fails (e.g., if coordinates are jointly entangled in a way that distorts individual marginals), Theorem 4's empirical numbers do not transfer. A precise analysis of Regev's measurement marginals from the original paper is left to future work.

### 3.5 Joint interpretation

Theorem 1 says: *once* `r_a | L`, success is deterministic regardless of noise. Theorem 2 says: an ideal algorithm reaches this state in `O(log log log N)` bases in expectation. Theorem 3 says: under destructive noise, the actual algorithm reaches it in `E[K_λ^{ideal}] / g_M(η)` bases — i.e., overhead exactly `1/g_M(η)`. Theorem 4 says: under a marginal-distribution assumption, the (C) framework applies coordinate-wise to Regev's multi-base measurements with corresponding reduction in run count. Together they explain the algorithm's empirical robustness: a noise-adjusted logarithmic number of measurements suffice to enter a regime immune to further measurement noise, and the framework composes naturally with Regev's multi-base circuit (modulo marginal assumption).

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

**Position of our contribution.** Every individual ingredient — `r | λ(N)`, multi-base lcm, divisor search given exponent, continued fractions, factor_from_exponent — is folklore or established in the references above. The *clean separation theorem* (noise invariance via the classical verification gate, formalized as Theorem 1) and its corollaries appear to not be stated explicitly. We propose Theorem 1 as the *folklore packaged* result that gives a uniform explanation for the empirical robustness of multi-base post-processing.

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
python demo.py --noise3 77 143              # reproduces 3-noise comparison
python demo.py --adaptive 77 143 209 323    # reproduces §9.1 negative result
```

The codebase is `~700` lines of numpy. No quantum libraries required.

## 8. Conclusion

We formalize a folklore observation about classical post-processing of multi-base Shor order-finding measurements: augmenting continued-fraction recovery with a divisor search over the accumulated lcm `L` makes the procedure deterministic whenever the order of the current base divides `L` — independent of the measurement distribution, hence robust to all measurement-layer noise models. We prove two complementary quantitative bounds: Theorem 2 on the number of bases `K_λ` needed to reach the noise-invariant regime (refining a corollary of Pomerance et al. 2017), and Theorem 3 characterizing the exact noise overhead `1/g_M(η)` for destructive noise models. We verify Theorem 1 across 17,700 trials at six composite sizes (zero violations), Theorem 2 across 17 semiprimes × 1,000 trials, and Theorem 3 across 9 noise setups on N=437, and show that two natural attempts to extend the result (adaptive base selection, lattice post-processing) yield no further measurement savings, indicating Theorem 1 captures the essential structure at this scale.

The work is best viewed as expository and a foundation for explorations of more substantive improvements: noise-resilient implementations on hardware, the `lucky` region under structured measurement distributions, and connections to lattice-based factoring (Regev 2023).

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
11. Regev, O. (2023). *An efficient quantum factoring algorithm*. arXiv:2308.06572.
12. *A quantum algorithm for computing the Carmichael function* (2021). arXiv:2111.02488.

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
