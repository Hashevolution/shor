# §3.6 Amplification 단락 초안 (대기 중)

**Status**: 결과 받는 즉시 paper §3.6 에 통합.

## 시나리오 A: 확실한 증폭 (thinned SR > 5%)

```markdown
**Engineered amplification demonstration.** To verify that the small magnitude
of the SR signal at (N, d) = (437, 4) is a property of the algorithm's
*well-functioning* regime rather than a property of the underlying mechanism,
we constructed a deliberately sub-functional variant of the hybrid algorithm
("thinned hybrid"): each measurement's convergent-recovery pool is restricted
to the single smallest convergent-fraction denominator, and the (C)
divisor-search augmentation is disabled. This yields a baseline that succeeds
substantially less often than the original hybrid (K_baseline = [TBD: thinned
K_base] vs 1.92 in the original, a [TBD]× increase). At this sub-functional
baseline, the same phase-noise mechanism produces an amplified SR signal:
mean SR at σ = 0.050 is [TBD: thinned SR%], a factor [TBD: ratio]× larger
than the +0.144% mean in the original (3 seeds × 100 trials). This confirms
the structural prediction that the noise-as-resource paradigm is in its
*function-restoration* regime here, with the magnitude bounded only by how
sub-functional the baseline is (here engineered, not natural). The amplified
demonstration is purely conceptual: the thinned algorithm is intentionally
worse than the original and has no practical advantage. It shows that the
underlying boundary-flip mechanism is the same as in the original hybrid,
just operating on a larger pool of borderline trials.
```

## 시나리오 B: 약한 증폭 (thinned SR ~ 1-3%)

```markdown
**Engineered amplification test.** We tested whether the small magnitude of
the SR signal at (N, d) = (437, 4) reflects the algorithm's *well-functioning*
regime by constructing a deliberately sub-functional variant ("thinned hybrid":
each measurement uses only the smallest convergent-fraction denominator, and
(C) divisor search is disabled). The thinned baseline has K_baseline = [TBD]
(vs 1.92 in the original), confirming that the algorithm is substantially
less efficient without the full convergent pool. The SR magnitude at σ = 0.050
is [TBD: thinned SR%] in the thinned variant vs +0.144% in the original ([TBD]×
amplification). The amplification is modest because the dominant contribution
to K_baseline in the thinned variant is now from trials that fail entirely
(K = max_runs), not from borderline K = 2 or K = 3 trials. The mechanism is
the same; the magnitude is bounded by the *boundary-trial population*, which
does not scale linearly with K_baseline. Larger amplification would require
constructing a baseline that is sub-functional *specifically near the K = 1 /
K = 2 boundary*, an engineering not pursued here.
```

## 시나리오 C: 차이 없음 (thinned SR ≈ normal SR)

```markdown
**Engineered amplification test (null result).** We tested whether the small
magnitude of the SR signal at (N, d) = (437, 4) reflects the algorithm's
*well-functioning* regime by constructing a thinned variant (smallest
convergent only, no (C) divisor search). The thinned baseline has K_baseline
= [TBD] (vs 1.92 in the original), confirming sub-functional operation. The
SR magnitude is, however, comparable to the original ([TBD] vs +0.144%),
showing that *post-processing thinning alone does not amplify the SR signal*.
This suggests the SR magnitude is bounded by the *base-set-specific K-
distribution near the active boundary*, not by the algorithm's overall
success rate. To amplify the effect, the base set's borderline K-population
would need to be enlarged — a base-set-engineering question we leave open.
```

## 시나리오 D: thinned 의 K_base 너무 커서 noise floor

```markdown
**Engineered amplification test (regime change).** Our initial thinned variant
(smallest convergent only, no (C) divisor search) produced K_baseline = [TBD],
which is in the *noise-floor regime* where per-seed variance dominates the
SR effect. We did not detect an amplified mechanism signal at this scale
([TBD] result). A more refined sub-functional baseline that retains K = 1 /
K = 2 boundary population while reducing the overall success rate is needed
to demonstrate amplification.
```

---

## (1147, 2) cross-cell 통합 단락 (별도)

```markdown
**Cross-cell verification at (N, d) = (1147, 2).** To test mechanism
universality across different cells in the active-boundary regime, we ran
the σ-scan with 5 seeds × 100 trials × 5 σ values at (N, d) = (1147, 2),
where K_baseline = [TBD final mean]. Results:
- [TBD: x/5] of 5 seeds show boundary-flip activity (mechanism universality)
- Boundary distribution: [TBD]% K = 1 / K = 2, [TBD]% K = 2 / K = 3,
  [TBD]% K = 3 ↔ K = 1 long-jump
- σ-curve plateau + decline pattern confirmed
- Cross-seed mean SR at σ = 0.050: [TBD ± TBD]% ([TBD]t, [TBD]p)
- σ-curve direction asymmetry [TBD: observed / not observed]

The mechanism observation generalizes to this cell, with the [TBD: same /
shifted] boundary distribution we predicted given the larger K_baseline
relative to (437, 4). [TBD: optional: comparison with (437, 4) result.]
```
