# Release notes — v0.5.1 (2026-06-19)

**DOI**: [10.5281/zenodo.20767685](https://doi.org/10.5281/zenodo.20767685)

## Summary

**v0.5.1** extends the magic (nonstabilizerness) direction of v0.5.0 (DOI
[10.5281/zenodo.20725965](https://doi.org/10.5281/zenodo.20725965)) with a
**coding-theory of marker sets**, an **oracle-hiding = T-cost** identification with a
fault-tolerant resource estimate, and — most importantly — a **full-text prior-art audit of all
four code↔magic literatures** that led us to **honestly credit the flat-state closed form to prior
work**. The order-finding results (Theorems 1–6) are unchanged; the magic direction is now
integrated into the canonical `paper.md` as §9.

## What's new

### Coding-theory of marker sets (Props 4, 5, 5′)
For a flat marker state $|{\rm flat}_W\rangle=|W|^{-1/2}\sum_{x\in W}|x\rangle$:

- **Prop 4 (closed form = additive energy):**
  $M_2=-\log_2\big(M^{-4}\sum_x E(W\cap(W{\oplus}x))\big)$, $E$ = additive energy. Reduces to $0$ for
  $M=1$ and affine $W$; exact zero-test $M_2=0\iff A_W\in\{0,M\}$. **Credited to prior work** (see below).
- **Prop 5 (Sidon law):** generic (Sidon) marker sets give $M_2=\log_2\frac{M^3}{7M-6}\to 2\log_2M-\log_27$.
- **Prop 5′ (exact random expectation):**
  $\mathbb E[\xi]M^4=(7M^2-6M)+\frac{7(M)_4}{N-3}+\frac{N(N-1)(N-2)(N-4)(M)_8}{(N)_8}$ (closed form
  vs Monte Carlo $\le1.4\times10^{-2}$).
- **Metric correction:** minimum Hamming distance does **not** determine magic ($\{0,1,2,3\}$ vs
  $\{0,1,2,4\}$: same $d_{\min}=1$, $M_2=0$ vs $1.54$).

### Oracle-hiding = T-cost (Prop 6) + FTQC estimate
- $M_2(|\psi_f\rangle)>0\iff f$ nonlinear (degree-$\ge2$ ANF) $\iff U_f$ needs non-Clifford ($T$)
  gates; all zero iff $f$ affine.
- **FTQC concretization** (`oracle_ftqc_estimate.py`): degree-$d$ monomial $\to(2d{-}3)$ Toffolis
  $\to7T$ (or $4T$). $T_{\rm est}=0\iff M_2=0$; Shor modexp ANF-synthesis upper bound $N{=}15\to4$
  Toffolis $\sim N{=}35\to4171$.

### Grover ladder (carried from v0.5.0, now in `paper.md` §9)
Simon $0$ / Grover finite (peak exactly 3 bits, density $\to0$) / Shor $\to L$.

## Honest novelty audit (the headline of this release)

We obtained and full-text compared the four nearest code↔magic literatures:

| group | prior work | verdict |
|---|---|---|
| (D) flat/SMF SRE closed form | **Tarabunga–Castelnovo, *Quantum* 8, 1347 (2024), Eq. 8** | **Prop 4 is its uniform-support specialization — credited, not new** |
| (B) hypergraph-state magic | Chen–Yan–Zhou (*Quantum* 8, 1351); Kagamihara–Tsuchiya (2602.23687) | phase-encoded / RM(2) ≠ our support / RM(1) — differentiated |
| (A) weight-enumerator SRE | Quantum Lego (2308.05152) | tensor-network SRE *tool*, not marker-set coding-theory |
| (C) Dicke/permutation-invariant | Passarelli–Fazio–Lucignano (2402.08551) | symmetric-support only |

Also: the "$2\log_2 M$" growth is the saturation of the standard bound $M_\alpha\le2\log_2 R$; Shor's
magic↔period law is Paviglianiti et al. (2605.05347). **Surviving defensible novelty:** the Grover
3-bit ladder, the *coding-theoretic specialization* (zero-test, exact Sidon constant, exact random
expectation, $d_{\min}$ refutation, Grover application), and oracle-hiding = $T$-cost (Prop 6).

## What changed

### New files
- `experiments/marker_code_magic.py` — autocorrelation/$\tau$ indicator, $d_{\min}$ refutation.
- `experiments/marker_code_closed_form.py` — Prop 4 (additive energy), Prop 5 (Sidon law).
- `experiments/marker_code_expected.py` — Prop 5′ (exact $\mathbb E[\xi]$).
- `experiments/oracle_tcount_magic.py` — Prop 6 (state magic ⟺ oracle $T$-cost).
- `experiments/oracle_ftqc_estimate.py` — Prop 6 FTQC $T$-count estimate (honest upper bound).
- `magic-paper-draft.md` — standalone draft (contributions box + practical implications).
- `magic-program-overview.md` — integrated overview + M1–M5 roadmap.
- `HISTORY.md`, `HISTORY-쉬운설명.md`, `M5-쉬운설명.md` — chronology + lay explanations.

### Updated
- `magic-results.md` — Props 4, 5, 5′, 6 (with prior-work credits).
- `magic-prior-art.md` — §5c: full-text audit of (A)(B)(C)(D), novelty table, citations.
- `paper.md` — new §9 "Companion direction: magic across the speedup ladder" + abstract pointer +
  references [15]–[22]; `paper.tex` remains frozen at v0.2.1. **§9 also includes the D3 finding
  paragraph** (asymptotic tightness of Prop 5′, see below).

### Added late in v0.5.1 (post-initial-draft, before Zenodo mint): JAMES-DISCOVER + D3
- `ai-discovery-engine-design.md` — design + feasibility review for an automated discovery loop
  (Generator → Probe → Miner → Adversary → Promoter) layered on the existing magic infrastructure.
  numpy-only; no LLM in the core loop (kept as optional D3+ extension).
- `experiments/discover_poc.py` — D1 sanity gate. Re-discovers the Sidon constant
  $M_2=\log_2(M^3/(7M-6))$ and the additive-energy closed form (Props 4/5) from scratch.
- `experiments/discover_d3_jensen.py` — **D3 result**, closing one of v0.5.1's listed open items.
  Defines the Jensen gap $J(M,N):=\mathbb E[M_2]+\log_2\mathbb E[\xi]\ge 0$ (Prop 5′ closes only
  $\mathbb E[\xi]$, leaving $J$ uncontrolled) and finds $J\propto 1/N$ in the sparse regime
  $M^2\ll N$ (measured slope $-0.96$ across $M\in\{6,8,10\}$, $n\in\{10,\dots,13\}$). Consequence:
  **$-\log_2\mathbb E[\xi]$ is not merely a Jensen lower bound but an asymptotically tight estimate
  of $\mathbb E[M_2]$, with absolute error $O(1/N)$.** Also records the saturation boundary
  $M^2/N\to 1$ where the $1/N$ law breaks (additive-collision saturation), and reports
  $\kappa(M):=J\cdot N$ as having no simple closed form in the elementary dictionary
  $\{M,M^2,M(M-1),M^3,M\log_2 M,1\}$ — kept as an open sub-question rather than over-fit.

### Verification
- Regression suite `magic_proofs_check.py`: **42 assertions pass**.
- SRE tool cross-checked to $10^{-15}$; closed forms (Props 4/5/5′) verified to $\le1.4\times10^{-2}$
  (Monte Carlo) / $10^{-15}$ (deterministic).

## What this is NOT
- Not a new quantum speedup, algorithm, or cryptographic result.
- Not a production FTQC estimator — `oracle_ftqc_estimate.py` is an explicit upper bound (arithmetic
  oracles are far cheaper via windowed synthesis; cf. Gidney–Ekerå 2021).
- Prop 4's closed form is **not** claimed as new (credited to Tarabunga–Castelnovo).

## Honest scope / remaining
- (A) full-text checked for 2308.05152 only; the distillation refs (2510.10852, 1702.06990,
  2501.10163) are *distillation*-side and do not affect our claims.
- ~~Open: the $M^2\gtrsim N$ regime of $\mathbb E[\xi]$ is exact in closed form but its asymptotics are
  characterized only to leading order.~~ **(D3 partial closure:** in the sparse regime $M^2\ll N$,
  Prop 5′ is shown asymptotically tight with absolute error $O(1/N)$. The $M^2/N\to 1$ saturation
  boundary remains open, and the prefactor $\kappa(M)=J\cdot N$ remains open in closed form.**)**

## Citing
Cite v0.5.1 for the marker-set coding-theory + oracle/T-cost, alongside v0.5.0 (magic direction +
ladder) and the order-finding releases (v0.2.1 theorems, v0.4.0 Theorem 6). Prior work credited
in-text: Tarabunga–Castelnovo (*Quantum* 8, 1347), Paviglianiti et al. (2605.05347),
Leone–Oliviero–Hamma (PRL 128, 050402).
