# Contact draft — Yang & Markidis (KTH)

**Recipients**:
- Qingxin Yang (qingxiny@kth.se)
- Stefano Markidis (markidis@kth.se)

**Suggested subject**:
> Analytical complement to arXiv:2605.16074: closed-form σ-curve for noisy Shor's algorithm + possible data collaboration

**Body** (formal academic English):

---

Dear Dr. Yang and Prof. Markidis,

I read your recent paper "When Noisy Quantum Order Finding Remains Recoverable for Shor's Algorithm" (arXiv:2605.16074, ICS Workshops '26) with great interest. Your empirical analysis of 680 IBM precision-register distributions and the four interpretable features (A_peak, H_norm, M_{1,frac}, Δ_{ver,frac}) is, to my knowledge, the first systematic characterization of recoverability in noisy NISQ order finding.

I am writing to share a complementary analytical work that I have just released as v0.3.0 of my Shor-class analysis repository:

- **DOI**: 10.5281/zenodo.20685015
- **GitHub**: https://github.com/Hashevolution/shor
- **Short paper draft**: see `arxiv_draft.md` in the repository (preparing for ICS Workshops '27 or Quantum journal).

**Brief summary of the complementary contribution**:

Under per-amplitude Gaussian phase noise of magnitude σ, the noise-averaged measurement distribution satisfies

  E[|FFT(a·e^{iε})_y|²] = (1 − e^{−σ²}) / Q + e^{−σ²} · P_0(y),

i.e. exactly the structural form of your two-stage noise propagation model
`(1 − ε) · P_s + ε · distractors`, with the analytical identification
**ε = 1 − exp(−σ²)**.

This closed form predicts the success probability of standard continued-fraction post-processing as

  p(σ) = ρ + (p_0 − ρ) · exp(−σ²),

which we have verified across five algorithm classes:

| Algorithm                          | R²    |
|------------------------------------|------:|
| Grover (k iterations)              | 0.88  |
| Shor pure (with b-trick)           | 0.95  |
| QPE isolated (no b-trick)          | 0.96  |
| Simon's algorithm                  | 0.99  |
| Hybrid (C)+b-trick (our setup)     | 0.91  |

Additionally, the same form extends to depolarizing (R² = 0.995) and bias-zero (R² = 0.996) noise (paper's universal form `p = (1−ε)·p_0 + ε·g_∞`), while amplitude damping remains structural (R² = 0.03 — outside the form, requiring model-specific treatment).

Most relevant to your work: the four Yang-Markidis features admit (partial) closed-form predictions under the same noise-averaged distribution:

- **A_peak(σ) = A_peak(0)** (σ-invariant, because q_σ = e^{−σ²}·q_0 cancels in the ratio).
- **M_{1,frac}(σ)** and **Δ_{ver,frac}(σ)** are rational functions in u = exp(−σ²).
- **H_norm(σ)** monotonically increases (no clean closed form, but limits are explicit).

I have verified the A_peak invariance and the H_norm monotonicity in numerical experiments at N = 437.

**Possible collaboration / data sharing**:

I would be very interested in two follow-ups, if you are open to them:

1. **Fitting the closed form to your 680 IBM distributions**. The model predicts ε per run from the broadening of P_σ relative to P_0; this is a single-parameter analytical test of the framework on real hardware noise. If your dataset is shareable (or you can share aggregate statistics), I would be happy to provide the fitting code.

2. **A joint short paper** (or back-to-back) at ICS Workshops '27 or Quantum journal, positioning the empirical (your) and analytical (mine) results as complementary halves of a single framework for noisy NISQ Shor.

I also want to flag — for full transparency — that the closed form led me to **self-correct** a "boundary-flip mechanism" claim I had made earlier in v0.2.1 (DOI 10.5281/zenodo.20681847, §3.6). That earlier framework is now retracted in favor of the closed form; full details are in `release_notes_v0.3.0.md` and §3.6.bis of `paper.md` in the repository.

I would be glad to hear your thoughts, and happy to set up a video call if useful.

Best regards,

[Your name / affiliation]
Hashevolution
https://github.com/Hashevolution/shor

---

**Notes for sender**:

1. Replace "[Your name / affiliation]" with your real name and affiliation. If you prefer to remain anonymous as "Hashevolution," that is academic-norm-acceptable for an arXiv-like work but slightly lowers response rate; consider a real name if comfortable.
2. If you have an institutional email, send from that; cold email from a personal address is also fine but tag the GitHub repo for credibility.
3. Attach a PDF of the v0.3.0 release notes or the arxiv_draft.md (after LaTeX conversion) for quick review.
4. Send to **both** authors — Yang likely handles day-to-day correspondence, Markidis is the senior author who decides on collaborations.
5. Expected response time: 1-3 weeks. No response within 4 weeks → follow-up once is acceptable.

**If they respond positively**:

- Share GitHub repo access (it is already public).
- Ask for IBM dataset (anonymized OK).
- Discuss venue (ICS Workshops '27 deadline is typically March-April; Quantum journal rolling).
- Discuss authorship if collaboration agreed (joint paper with both works or back-to-back works).

**If they decline or no response**:

- Proceed independently. Cite their work in your arXiv submission.
- Approach Quantum journal directly (no endorsement needed).
- Workshop venue: NISQ workshops or quantum information workshops are alternatives.
