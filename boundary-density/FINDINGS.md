# Boundary density (rho_b) — K-distribution geometry as the structural observable

A working note that formalizes a quantity left implicit in §3.6: the **boundary
density** `rho_b`, the fraction of trials sitting on a flippable K-bin boundary.
It is the structural driver behind the *per-seed SR magnitude* distribution, and it
is computed **without any noise** — directly from the noise-free K-histograms.

Run: `python boundary-density/boundary_density.py` (pure stdlib, reads `../experiments/*_histograms.txt`).

## 1. The headline: `K_baseline` is a misleading proxy for boundary density

Decompose the mean:

```
mean K  =  (dead / ceiling mass, K >= max_runs)
         + (floor mass, K = 1)
         + (live boundary structure, 2 <= K < max_runs = rho_b)
```

Only `rho_b` is the population noise can flip. Across all measured cells
(noise-free baseline, σ = 0):

| structure / cell | meanK (`K_baseline`) | dead% | floor% | **rho_b** | meanK(live) | observed per-seed \|SR\| (paper §3.6) |
|---|---:|---:|---:|---:|---:|---|
| pure_shor (d=1)        | **10.38** | 47.8 | 34.0 | 18.2% | 1.57 | 0–1.10% (narrow) |
| pure_regev (LLL)       | **2.44**  |  6.4 | 76.0 | 17.6% | 1.24 | 1.27–3.95% |
| sigma_scan_437_d4      | 1.84      |  4.1 | 90.2 |  5.7% | 1.06 | 0–1.93% |
| sigma_scan_1147_d2     | 2.92      |  5.2 | 54.6 | **40.2%** | 1.98 | up to +9.44% (high-K rescue) |
| sigma_scan_2491_d2     | 2.30      |  5.0 | 73.8 | 21.2% | 1.37 | −4.89% single → +0.30% multi |

**Result A — the mean lies.** pure_shor `K_baseline = 10.38` vs pure_regev `2.44`
(4× gap), but `rho_b` is ~18% for *both*, and meanK(live) is nearly equal
(1.57 vs 1.24). The entire 4× gap is Shor's ~48% dead-trial mass (b-trick
nontrivial-sqrt failures pinned at max_runs) — unflippable ceiling mass that
inflates the mean while contributing **zero** boundary density. This makes
precise the paper's remark that Pure Shor's small SR is from "dead-trial
dominance."

## 2. The positive signal: rho_b tracks SR magnitude (★)

Ordering cells by `rho_b` largely reproduces the ordering of observed per-seed
|SR|:

- `(437, 4)` rho_b = 5.7% → smallest |SR| (0–1.93%)
- `(1147, 2)` rho_b = 40% → largest |SR| (high-K rescue, up to +9.44%)

This is direct evidence that `rho_b` — a noise-free structural quantity — predicts
where the SR mechanism has room to act. It is exactly the "borderline-population is
the bottleneck" prediction (§3.6 engineered-amplification), now measured as a
geometric property of the distribution rather than inferred from thinning.

## 3. The open nuance: rho_b alone is not sufficient — shape matters too

pure_shor and pure_regev have nearly equal `rho_b` (18.2% vs 17.6%) but Shor's
|SR| is *narrower* (0–1.10% vs 1.27–3.95%). The difference is **where inside the
boundary band** the live mass sits: Regev/`(1147,2)` carry higher meanK(live)
(more spread into K = 2…8, enabling long jumps / high-K rescue), while Shor's
boundary trials concentrate at K = 2–3 just above the floor. So the right
structural observable is `rho_b` **plus the within-band K-spread**, not `rho_b`
alone. This is the concrete next-step question.

## 4. Why this is a clean, conflict-free continuation of §3.6

- It is the direct answer to **open question 2** ("which K-boundary is active is
  not predictable from `K_baseline` alone; finer-grained analysis needed").
- It *refines*, does not contradict, the regime map: ordering by mean K conflates
  dead mass with boundary structure; ordering by `rho_b` does not. (Fits the
  paper's existing "retracted/refined" register.)
- It decouples from the null, expensive σ-sweep: `rho_b` is computed from σ = 0
  histograms alone. Noise is just the kernel that turns `rho_b` into observable SR.

## 5. Plan (pending the in-progress algorithm-generalization experiments)

1. Finish the generalization runs; have each emit `*_histograms.txt` in the
   standard `seed/sigma/K/count` format (this analyzer ingests them as-is).
2. Compute `rho_b` + within-band spread for every structure; build the
   structure → boundary-sensitivity map.
3. **Decision gate:** if `rho_b` (and spread) vary *systematically and
   predictably* with algorithm structure → grow into its own line (separate repo
   `boundary-density`, already scaffolded and bundled). If it scatters like the
   mean K → fold back into shor as a one-paragraph refinement of §3.6.

### Determinants to test
- **base count d / parallelism** → success rate → floor mass vs rho_b trade-off
- **post-processing** (LLL vs continued-fraction vs (C) divisor augmentation; the
  (C) buffer is already shown to absorb borderline trials — mild-thin experiment)
- **dead-trial generator** (b-trick failure rate) → ceiling mass (mean-only, no rho_b)

## Provenance

A standalone scaffold of this study (analyzer + engine + baseline data + README)
is parked as a portable git bundle (`boundary-density.bundle`) pending the
location decision in step 3. This in-repo version reuses shor's existing
`experiments/` data and engine and adds only the analyzer + this note.
