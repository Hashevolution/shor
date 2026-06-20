# Magic across the quantum speedup ladder: a Grover saturation, a coding-theory specialization, and oracle-hiding as T-cost

*Working draft (v0.5.0+). Full proofs in `magic-results.md`; prior art in `magic-prior-art.md`;
reproducible code in `experiments/` (`magic.py` + `magic_proofs_check.py`, 42 assertions).*

---

## Contributions vs prior art (honest summary)

After full-text comparison with the four nearest code↔magic literatures, our defensible
contributions and what we **credit to prior work** are:

**New (this work).**
1. **Grover/polynomial-speedup magic** — closed form for single-marked Grover, peak saturating at
   **exactly 3 bits**, density $M_2/n\to0$, returned to $0$ at the solution; the speedup **ladder**
   (Simon $0$ / Grover finite / Shor $\to L$) (Props 2–3, Cor 1–2).
2. **Coding-theoretic *specialization* of marker sets** — the exact zero-test
   $M_2=0\iff A_W\in\{0,M\}$, the **Sidon law** with exact constant $\log_2 7$ and finite-$M$ form
   (Prop 5), the **exact random expectation** $\mathbb E[\xi]$ (Prop 5$'$), the $d_{\min}$ refutation,
   and the **Grover multi-marked application** (Prop 2$'$). The Jensen step in Prop 5$'$ is shown
   asymptotically tight in the sparse regime via an automated rediscovery loop:
   $J(M,N):=\mathbb E[M_2]+\log_2\mathbb E[\xi]\propto 1/N$ for $M^2\ll N$ (§3.1).
3. **Oracle-hiding = $T$-cost** — $M_2(|\psi_f\rangle)>0\iff f$ nonlinear $\iff U_f$ needs $T$ gates
   (Prop 6).

**Credited to prior work (not claimed).**
- Shor's magic↔period quantitative law — Paviglianiti et al., arXiv:2605.05347.
- The flat-state SRE **closed form** (our Prop 4 via additive energy) — uniform-support case of
  Tarabunga–Castelnovo's RK/SMF formula, *Quantum* **8**, 1347 (2024), Eq. (8).
- The $2\log_2 M$ growth rate — saturation of the standard bound $M_\alpha\le2\log_2 R$.
- SRE definition/tool (Leone–Oliviero–Hamma); XOR-FWHT computation (arXiv:2512.24685).

**Nearest neighbors differentiated.** Hypergraph-state magic (phase / RM(2): Chen–Yan–Zhou *Quantum*
**8**, 1351; Kagamihara–Tsuchiya arXiv:2602.23687) vs ours (support / RM(1)); weight-enumerator SRE
tools (Quantum Lego arXiv:2308.05152); permutation-invariant/Dicke magic (arXiv:2402.08551,
symmetric support only).

---

## Abstract

Nonstabilizerness ("magic") is necessary for quantum computational advantage, but how *much*
magic each speedup *type* requires has remained largely qualitative away from Shor's algorithm. We
quantify the **amount and density** of magic, measured by the stabilizer 2-Rényi entropy $M_2$,
across a ladder of algorithms on this repository's state-vector simulator. We show: (i) for the
oracle/intermediate states of Simon-type algorithms the magic is identically zero (affine oracle =
Clifford) despite an exponential *query* speedup; (ii) Grover's quadratic speedup uses only a
**finite** amount of magic — a single-marked Grover state's $M_2$ saturates at exactly $3$ bits and
its density $M_2/n\to0$, and the magic is *returned to zero* at the solution; (iii) Shor's magic
grows with problem size (density $\to1$), reproducing and contrasting with Paviglianiti et al.
[arXiv:2605.05347]. We then develop a **coding-theory of marker sets**: the magic of a flat
("uniform-amplitude, uniform-phase") superposition over a set $W\subseteq\mathbb F_2^n$ is governed
*not* by the minimum Hamming distance but by the **additive energy** of the shifted self-intersections
of $W$ (the uniform-support case of Tarabunga–Castelnovo's Rokhsar–Kivelson SRE formula, which we
credit). Specializing it to a coding-theory of marker sets yields (a) the exact zero-test
$M_2=0\iff A_W\in\{0,M\}$, (b) the law $M_2=\log_2\!\frac{M^3}{7M-6}\approx2\log_2 M-\log_2 7$ for
generic (Sidon) marker sets, and (c) an exact expression for the expected $M_2$ of a random marker set. Finally we unify oracle-hiding with gate cost: the magic a query
oracle hides equals the non-Clifford ($T$) cost of its gate decomposition; both vanish iff the
oracle function is affine. All claims are verified numerically to $\le10^{-13}$ and by a 42-assertion
regression suite.

---

## 1. Introduction

The Gottesman–Knill theorem makes magic a *necessary* resource for quantum advantage: Clifford
circuits on stabilizer states are classically simulable. A sharp recent result (Paviglianiti,
Seclì, Tirrito, Savona, arXiv:2605.05347) links the *amount* of magic in Shor's algorithm to the
number-theoretic hardness of the period. We ask the complementary, cross-algorithm question:

> **Does the *type* of quantum speedup (none / polynomial / exponential) correspond to the
> *amount and density* of magic an algorithm spends?**

Our headline is a **speedup ladder** (Table 1): magic density $M_2/n$ separates Clifford-trivial
($0$), polynomial (finite, density $\to0$), and exponential (density $\to1$) speedups, and magic
lives not on the visible circuit surface but in the *nonlinearity of the problem*, which oracles and
FFT shortcuts can hide.

We use the stabilizer 2-Rényi entropy (Leone–Oliviero–Hamma)
$$M_2(|\psi\rangle)=-\log_2\Big(\tfrac1{2^n}\sum_{P\in\mathcal P_n}\langle\psi|P|\psi\rangle^4\Big),
\qquad 0\le M_2\le n,\quad M_2=0\iff|\psi\rangle\text{ stabilizer},$$
computed exactly in $O(n\,4^n)$ via an XOR–Walsh–Hadamard identity (`magic.py`), cross-checked
against brute-force Pauli enumeration to $10^{-15}$.

---

## 2. Magic across the speedup ladder

**Lemma 1 (flatness ⇒ stabilizer).** A uniform-amplitude, uniform-phase state on support
$S\subseteq\mathbb F_2^n$ has $M_2=0$ iff $S$ is an affine subspace.

**Corollary 1 (oracle-hiding).** The function graph state $|\psi_f\rangle=2^{-n/2}\sum_x|x\rangle|f(x)\rangle$
has $M_2=0$ iff $f$ is $\mathbb F_2$-affine. Hence **Simon** (a linear 2-to-1 oracle suffices) keeps
$M_2\equiv0$ while enjoying an exponential *query* speedup — magic is *not* required for query
advantage — whereas **Shor**'s multiplicative $f(x)=a^x\bmod N$ forces $M_2>0$ growing with size.

**Propositions 2–3 (Grover).** A single-marked Grover state has a closed form for $\sum_P\langle
P\rangle^4$; as $N\to\infty$ the peak magic $\to-\log_2(a^8+(1-a^2)^4)$, maximized at $a^2=\tfrac12$
to **exactly $3$ bits**, so density $M_2/n\to0$, and $M_2$ returns to $0$ at the computational-basis
solution. (Numerics: peak $\to3.000$ at $n=30$.)

**Table 1 — speedup ladder.**

| algorithm | speedup | $M_2$ behavior | density $M_2/n$ | source |
|---|---|---|---|---|
| Simon / BV | exponential (query) | $0$ (affine oracle = Clifford) | $0$ | Cor. 1 |
| **Grover** | quadratic | peak $\to3$ bits, $0$ at answer | $\to0$ | Prop. 2–3 |
| Shor (comb) | exponential | $\propto t$ | $\sim0.4$–$0.55$ | Cor. 2 |
| Shor (in-circuit) | exponential | $\to L$ (max) | $\to1$ | arXiv:2605.05347 |

> **Thesis.** The amount/density of magic discriminates speedup type. Magic lives in the problem's
> nonlinearity (affine vs multiplicative); oracle and FFT black boxes can hide it.

---

## 3. A coding-theory of marker sets

For multiple marked items, the excess magic comes from the *structure of the marker set*
$W\subseteq\mathbb F_2^n$ ($M=|W|$). Writing $|{\rm flat}_W\rangle=M^{-1/2}\sum_{x\in W}|x\rangle$:

**Proposition 4 (additive-energy form; specialization of Tarabunga–Castelnovo).**
$$M_2(|{\rm flat}_W\rangle)=-\log_2\Big(\tfrac1{M^4}\sum_{x\in\mathbb F_2^n}E\big(W\cap(W{\oplus}x)\big)\Big),
\quad E(S)=\#\{(a,b,c,d)\in S^4:a{\oplus}b{\oplus}c{\oplus}d=0\}=\textstyle\sum_v A_S(v)^2,$$
the **additive energy** of the shifted self-intersections $S_x=W\cap(W\oplus x)$ (verified to
$4\times10^{-15}$). **This closed form is not new:** it is the uniform-support specialization of the
4-copy SRE formula of Tarabunga–Castelnovo (*Quantum* **8**, 1347 (2024), Eq. (8)), via the bijection
$(\sigma^{(1)},\sigma^{(2)},\sigma^{(3)},\sigma^{(4)})=(a,b,c,a{\oplus}b{\oplus}c{\oplus}x)$; we credit
it and use it as a tool. Our contribution is the **coding-theoretic specialization**: $M_2=0$ for
$M=1$ and for affine $W$ (Lemma 1), the exact zero-test $M_2=0\iff A_W\in\{0,M\}$ (two-valued
autocorrelation), and that the controlling invariant is the **autocorrelation spectrum** of
$\mathbf 1_W$, *not* the minimum Hamming distance: $\{0,1,2,3\}$ ($M_2=0$) and $\{0,1,2,4\}$
($M_2=1.54$) share $d_{\min}=1$.

**Proposition 5 (generic/Sidon law).** If $W$ is a Sidon ($B_2$) set,
$$M_2=\log_2\frac{M^3}{7M-6}\ \xrightarrow{M\to\infty}\ 2\log_2 M-\log_2 7.$$
A random $W$ is Sidon w.h.p. for $M\ll2^{n/2}$, so its magic concentrates on this value (matched to
$10^{-16}$). The $2\log_2 M$ growth rate is the saturation of the standard bound $M_\alpha\le2\log_2 R$
($R\le M$); the new content is the **exact constant $\log_2 7$ and finite-$M$ form**.

**Proposition 5′ (exact expected magic).** For a uniform random $M$-subset of $\mathbb F_2^n$ ($N=2^n$),
$$\mathbb E[\xi]\,M^4=(7M^2-6M)+\frac{7(M)_4}{N-3}+\frac{N(N-1)(N-2)(N-4)(M)_8}{(N)_8},\quad
\xi=2^{-M_2},$$
where the Sidon term comes from degenerate quadruples and the correction from genuine additive
quadruples; $\Delta\xi\to7(M)_4/(M^4N)$ for $M\ll N$. (Closed form vs Monte Carlo: $\le1.4\times10^{-2}$.)

### 3.1 Asymptotic tightness of Prop. 5$'$ (Jensen-gap law)

Prop. 5$'$ closes $\mathbb E[\xi]$ exactly, but the Jensen step $\mathbb E[M_2]\ge -\log_2\mathbb E[\xi]$
leaves the gap $J(M,N):=\mathbb E[M_2]+\log_2\mathbb E[\xi]\ge 0$ uncontrolled. We run an automated
rediscovery loop (`experiments/discover_d3_jensen.py`; design in `ai-discovery-engine-design.md`)
— which first passes a sanity gate by re-deriving the Sidon constant and the additive-energy closed
form from scratch (`experiments/discover_poc.py`) — and find $J\propto 1/N$ in the sparse regime
$M^2\ll N$ (measured slope $-0.96$ across $M\in\{6,8,10\}$, $n\in\{10,\dots,13\}$). Consequence:
**$-\log_2\mathbb E[\xi]$ is not merely a Jensen lower bound but an asymptotically tight estimate of
$\mathbb E[M_2]$, with absolute error $O(1/N)$.** The loop also records the saturation boundary
$M^2/N\to 1$ where the $1/N$ law breaks (additive-collision saturation) and reports the prefactor
$\kappa(M):=J\cdot N$ as having no simple closed form in the elementary dictionary
$\{M,M^2,M(M{-}1),M^3,M\log_2 M,1\}$ — kept as an open sub-question rather than over-fit.

**Relation to prior art.** The flat-state closed form (Prop. 4) is the uniform-support case of
Tarabunga–Castelnovo's 4-copy SRE formula for Rokhsar–Kivelson/SMF wavefunctions (*Quantum* **8**,
1347 (2024), Eq. (8)), which we credit; the $2\log_2 M$ rate is the standard $M_\alpha\le2\log_2 R$
bound. They apply it to *physical* SMF models (Ising, $J_1$–$J_2$, triangular AFM, spin glass) only.
Our contribution is the **coding-theoretic specialization**: the affine $\iff A_W\in\{0,M\}$
zero-test, the **exact Sidon law** (constant $\log_2 7$, finite-$M$ form, Prop. 5), the **exact random
expectation** (Prop. 5$'$), the $d_{\min}$ refutation, and the **Grover multi-marked application**. The
autocorrelation object also coincides with the geometric term $\Lambda$ of arXiv:2605.05347. This is
distinct from the nearest-neighbor literature on **hypergraph-state magic** (Chen–Yan–Zhou,
*Quantum* **8**, 1351 (2024); Kagamihara–Tsuchiya, arXiv:2602.23687), where the Boolean function is
**phase-encoded** ($\prod CZ_e|+\rangle^n=\sum_x(-1)^{f(x)}|x\rangle$) and the relevant code is
$\mathrm{RM}(2,n)$ (second-order nonlinearity / nonquadraticity); ours is **support-encoded** with
$\mathrm{RM}(1,n)$/affine. Reed–Muller$\leftrightarrow$magic in the *distillation* sense
(Campbell–Anwar–Browne; weight enumerators) is a different meaning, and Dicke/permutation-invariant
flat states cover only the symmetric-support case.

---

## 4. Oracle-hiding equals $T$-cost (Proposition 6)

**Proposition 6.** For the graph state $|\psi_f\rangle$ and oracle $U_f:|x\rangle|y\rangle\mapsto
|x\rangle|y\oplus f(x)\rangle$,
$$M_2(|\psi_f\rangle)>0\iff f\text{ has a degree-}\ge2\text{ ANF monomial}\iff U_f\text{ needs
non-Clifford (Toffoli/}T)\text{ gates},$$
all three vanishing iff $f$ is affine. (Linear ANF terms = CNOT/Clifford; quadratic = Toffoli $=T$.)
Numerically the count of nonlinear ANF monomials $T_{\rm proxy}$ and $M_2$ turn on and grow together
($0/0,1/1.54,2/2.48,3/3.70,4/4.43$); endpoints **Simon** ($T_{\rm proxy}=0,M_2=0$) vs **Shor**
modexp ($T_{\rm proxy}=4\to156$, $M_2=1.54\to5.23$). Thus *the magic a query oracle hides is exactly
the $T$-cost of its gate decomposition* — distinct from the Clifford-conjugation hiding of
Krüger–Mauerer (arXiv:2507.16543).

### 4.1 Practical implications (honest)

The three results sit at different points on the utility scale; we state this plainly.

- **Oracle-hiding $=$ $T$-cost (most practical).** $T$-count is the dominant cost in fault-tolerant
  quantum computing (each $T$ consumes one distilled magic state). Prop. 6 turns the *nonlinearity*
  of an oracle function into a fault-tolerant resource signal: under the explicit synthesis model
  "degree-$d$ ($\ge2$) ANF monomial $\to(2d{-}3)$ Toffolis $\to 7T$ each (or $4T$ with ancilla)"
  (`oracle_ftqc_estimate.py`), $T_{\rm est}=0\iff M_2=0\iff f$ affine, and $T_{\rm est}$ tracks $M_2$
  (Shor modexp ANF-synthesis upper bound: $N{=}15\to4$ Toffolis, $21\to580$, $33\to1375$,
  $35\to4171$). **Honest caveat:** per-output-bit ANF synthesis is an *upper bound* that ignores
  subterm sharing and arithmetic structure; real modexp uses windowed modular arithmetic with far
  lower $T$-count (Gidney–Ekerå 2021: RSA-2048 $\approx2.7\times10^9$ Toffolis). So the usable value is
  (i) the exact affine zero-test, (ii) the qualitative magic$\leftrightarrow T$ law, and (iii) a quick
  upper-bound indicator for *unstructured* oracles — not a replacement for production estimators.
- **Marker-set coding-theory (methodological).** Predicts/zeroes magic from cheap combinatorial
  invariants (autocorrelation, Sidon-ness) instead of an $O(n4^n)$ SRE evaluation; the zero-test
  $A_W\in\{0,M\}$ flags classically-easy (stabilizer) multi-solution structures; brings additive
  combinatorics to bear on magic estimation. The $d_{\min}$ refutation prevents use of the wrong metric.
- **Grover ladder (foundational/diagnostic).** Evidence that magic *density* classifies speedup type;
  useful for understanding the source of advantage and classical-simulability, not for cost reduction.
  Caveat: the 3-bit figure is a *state* quantity, not the gate-level $T$-count of a full Grover circuit.

---

## 5. Positioning and novelty

- **Primary prior results (cited, not claimed):** arXiv:2605.05347 settles Shor's magic↔period law
  (still a preprint); Tarabunga–Castelnovo (*Quantum* **8**, 1347 (2024)) own the flat/SMF SRE closed
  form (our Prop. 4 is its specialization); $2\log_2 M$ is the standard $M_\alpha\le2\log_2 R$ bound.
- **Genuinely new:** the Grover/polynomial-speedup magic trajectory and $3$-bit saturation; the
  **coding-theoretic specialization** of marker sets — affine$\iff A_W\in\{0,M\}$ zero-test, exact
  Sidon law (constant $\log_2 7$), exact random expectation (Prop. 5$'$), $d_{\min}$ refutation,
  Grover application; **asymptotic tightness of Prop. 5$'$** (Jensen-gap law $J\propto 1/N$, §3.1);
  oracle-hiding $=$ $T$-cost (Prop. 6).
- **Nearest neighbor differentiated:** hypergraph-state magic (phase / RM(2)) vs our support / RM(1).
- **Quantum-walk magic** (arXiv:2506.17783 PRR; arXiv:2504.19750 PRB 113,075142) is 1D-lattice
  transport, not search; cited as adjacent.

---

## 6. Reproducibility

- `magic.py` — SRE via XOR-FWHT, brute-force cross-checked ($10^{-15}$).
- `experiments/`: `grover_magic.py`, `shor_comb_magic.py`, `oracle_magic.py`,
  `oracle_tcount_magic.py`, `oracle_ftqc_estimate.py`, `marker_code_magic.py`,
  `marker_code_closed_form.py`, `marker_code_expected.py`.
- `experiments/magic_proofs_check.py` — **42 assertions** (Lemma 1, Cor 1, Props 2–3, 2′) all pass;
  peak convergence $n=10\to2.815$, $20\to2.994$, $30\to3.000$.
- Automated rediscovery loop (§3.1, design in `ai-discovery-engine-design.md`):
  `experiments/discover_poc.py` (D1 sanity gate: re-derives Sidon constant + additive-energy
  closed form), `experiments/discover_d3_jensen.py` (D3: Jensen-gap law $J\propto 1/N$). numpy-only;
  no LLM in the core loop.
- Zenodo record: v0.5.1 — DOI [10.5281/zenodo.20767685](https://doi.org/10.5281/zenodo.20767685).

## 7. Honest scope

Magic is necessary but not sufficient (it works with entanglement). The comb quantity is the
magic *held* by the post-measurement state, distinct from the magic *consumed* in-circuit
(arXiv:2605.05347). Full-text comparison of all four nearest code↔magic literatures is complete:
(A) Quantum Lego (arXiv:2308.05152) is a weight-enumerator/SRE tooling track, (B) hypergraph-state
magic (arXiv:2308.01886, 2602.23687) is phase-encoded / RM(2), (C) permutation-invariant Dicke
(arXiv:2402.08551) covers only the symmetric-support case, (D) Tarabunga–Castelnovo (Quantum **8**,
1347, 2024, Eq. 8) owns the flat-state closed form (credited; our Prop. 4 is its uniform-support
specialization). The surviving defensible contributions are the Grover 3-bit ladder, the
coding-theoretic specialization (Sidon constant, exact $\mathbb E[\xi]$, $d_{\min}$ refutation,
Grover multi-marked application, Jensen-gap tightness — see §3.1), and oracle-hiding $=$ $T$-cost.

## References (verified identifiers)

- Leone, S., Oliviero, S. F. E., Hamma, A. *Stabilizer Rényi entropy*. PRL **128**, 050402 (2022).
- **Tarabunga, P. S. & Castelnovo, C.** *Magic in generalized Rokhsar-Kivelson wavefunctions*.
  *Quantum* **8**, 1347 (2024). Eq. (8) is the 4-copy SRE formula for SMF wavefunctions; our Prop. 4
  is its uniform-support specialization (credited).
- Paviglianiti, A., Seclì, M., Tirrito, E., Savona, V. *The true cost of factoring: linking magic
  and number-theoretic complexity in Shor's algorithm*. arXiv:2605.05347.
- Chen, J., Yan, Y., Zhou, Y. *Magic of hypergraph states*. *Quantum* **8**, 1351 (2024).
- Kagamihara, R., Tsuchiya, A. *Nonstabilizerness of hypergraph states*. arXiv:2602.23687.
- Cao, C., Gullans, M. J., Lackey, B., Wang, Z. *Quantum Lego Expansion Pack: Enumerators from
  Tensor Networks*. arXiv:2308.05152. Weight-enumerator/SRE tooling track.
- Passarelli, G., Fazio, R., Lucignano, P. *Nonstabilizerness of Permutationally Invariant Systems*.
  arXiv:2402.08551. Symmetric-support (Dicke) case.
- Krüger, T., Mauerer, W. *Permutation-agnostic distance, Clifford-conjugation hiding*.
  arXiv:2507.16543. Adjacent but distinct from our oracle-hiding (nonlinearity).
- Mittal, S., Huang, J.-Y. *Quantum-walk magic*. PRR, DOI 10.1103/7rwg-lhpv. 1D transport, adjacent.
- *Quantum-walk magic on PRB*. PRB **113**, 075142 (2026); arXiv:2504.19750.
- *A fast and exact approach for stabilizer Rényi entropy via the XOR-FWHT algorithm*.
  arXiv:2512.24685.
- *State-vector SRE methods*. arXiv:2601.07824.

Full annotated list (including the (A)(B)(C)(D) full-text audit notes) in `magic-prior-art.md`.
