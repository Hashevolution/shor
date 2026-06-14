# boundary-density/

K-distribution geometry analysis for the §3.6 noise-as-resource program.

The net noise (SR) effect is null; the *per-seed magnitude* is algorithm-dependent.
This directory formalizes the structural quantity behind that magnitude — the
**boundary density** `rho_b`, the flippable K-bin population — and shows that the
mean K (`K_baseline`) is a misleading proxy for it (dead-trial mass inflates the
mean without adding boundary density).

- `boundary_density.py` — pure-stdlib analyzer; decomposes every
  `../experiments/*_histograms.txt` into dead / floor / `rho_b`.
- `FINDINGS.md` — results, the rho_b↔|SR| signal, the open shape nuance, and the
  research plan + decision gate for the algorithm boundary-sensitivity study.

```bash
python boundary-density/boundary_density.py
```

No dependencies for the analysis. Regenerate the underlying histograms via the
engine as in `paper.md` §3.6 (`python -m experiments.pure_shor_sr`, etc.).
