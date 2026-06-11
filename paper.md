# A Noise-Invariant Determinism Theorem for Multi-Base Post-Processing in Shor's Order Finding

*Draft, 2026-06-11.*
*Companion code: [github.com/Hashevolution/shor](https://github.com/Hashevolution/shor).*

## Abstract

We give an elementary observation about classical post-processing of multiple-base order-finding measurements in Shor's algorithm. Maintaining an accumulated exponent candidate `L` — the least common multiple of orders recovered from previous bases — and augmenting standard continued-fraction post-processing with a divisor search over `L` yields a procedure (which we call (C)) that recovers the order `r_a` of any new base `a` deterministically whenever `r_a | L`, *independent of the measurement distribution*. As a corollary, the procedure is robust to any measurement-layer noise model: a measurement returning a value of `k` drawn from any distribution still yields the correct order, provided the accumulated `L` covers `r_a`. We verify this experimentally across 11 noise setups (depolarizing, readout flip, bias, phase decoherence, amplitude damping, modular exponentiation error) at three composite sizes and 500 trials each (16,500 total measurements), observing zero violations. We discuss the theorem's relation to prior work — Knill's lcm trick (1995), McAnally's larger-Q convergent enumeration (2001), the Bach-Shallit textbook treatment of `r | λ(N)`, and the quantum algorithm for computing the Carmichael function (2021) — and conclude that while every individual ingredient is folklore, the *clean statement of the noise-invariance corollary* appears not to be made explicit in surveyed literature. We additionally show that two natural attempts to improve measurement count beyond this scheme — adaptive base selection and lattice-based joint post-processing — yield no further reduction, suggesting the theorem captures the essential structure of the problem at this scale.

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

## 3. Main theorem

**Theorem 1 ((C)-determinism).** If `r_a | L`, then `C(a, N, k, Q, L) = r_a` for every `k ∈ {0, ..., Q-1}`.

**Proof.** Assume `r_a | L`. Then `r_a ∈ divisors(L) ⊆ candidates(k, L)`. By construction, `a^{r_a} ≡ 1 mod N`, so `r_a ∈ valid`. Hence `valid ≠ ∅`. By (F2), every `d ∈ valid` is a positive multiple of `r_a`, so `r_a ≤ d`. Combining with `r_a ∈ valid` gives `min(valid) = r_a`. Finally, `minimize_order(a, N, r_a) = r_a` since `r_a` is already the order. ∎

**Corollary 1 (Noise invariance).** Theorem 1 makes no assumption on the distribution from which `k` is drawn. Hence: under any noise model that affects only the measurement step — including arbitrary distortions of the QFT output distribution — `(C)` returns `r_a` whenever `r_a | L`.

**Corollary 2 (L integrity).** When Algorithm 1 updates `L ← lcm(L, r)`, the value `r` returned by `(C)` is always either `r_a` (success) or `0` (failure). In the success case, `r_a | λ(N)` by (F1), so `L` always remains a divisor of `λ(N)`. Under noise, `L` cannot be corrupted to a value outside the divisors of `λ(N)`.

**Corollary 3 (Failure region).** `(C)` can fail only when `r_a ∤ L`. In that case, success depends on `convergents(k/Q)` containing some positive multiple of `r_a` — the classical "(B) path."

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

See `demo.py --verify <N>` to reproduce.

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

**The theorem is conditional on `r_a | L`.** When `L = 1` (start of algorithm), the condition fails and `(C) ≡ (B)`. The transient before `L` reaches `λ(N)` is governed by classical number theory: for typical semiprime `N = pq`, `E[K_λ] = O(log log N)` bases suffice (essentially proven by Pomerance et al. and the Carmichael-function paper).

**Classical verification gate is assumed exact.** Theorem 1 relies on `pow(a, d, N)` returning the correct value. On real hardware, this could in principle be corrupted by classical computation errors, but at standard machine precision this is negligible compared to quantum noise.

**Scale of verification.** Our experiments are at `N ≤ 437` (small for cryptographic purposes) but with high statistical replication (500 trials per noise setup). Scaling to cryptographic `N` requires no algorithmic change; the bottleneck is classical state-vector simulation. Hardware verification on existing NISQ devices is a natural next step.

## 7. Code, reproducibility

All experiments are reproducible via the companion repository:

```
pip install numpy
git clone https://github.com/Hashevolution/shor
cd shor
python demo.py --verify 77 143 209          # reproduces §4 table
python demo.py --noise3 77 143              # reproduces 3-noise comparison
python demo.py --adaptive 77 143 209 323    # reproduces §9.1 negative result
```

The codebase is `~700` lines of numpy. No quantum libraries required.

## 8. Conclusion

We formalize a folklore observation about classical post-processing of multi-base Shor order-finding measurements: augmenting continued-fraction recovery with a divisor search over the accumulated lcm `L` makes the procedure deterministic whenever the order of the current base divides `L` — independent of the measurement distribution, hence robust to all measurement-layer noise models. We verify this across 16,500 trials in 11 noise setups (zero violations), and show that two natural attempts to extend the result (adaptive base selection, lattice post-processing) yield no further measurement savings, indicating Theorem 1 captures the essential structure at this scale.

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
