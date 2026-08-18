# Exceptional points and cusp singularities in relativistic stellar oscillations with causal bulk relaxation

Reproducibility repository for the manuscript (CQG-115250).

## Summary

Causal bulk relaxation (Israel–Stewart) promotes the stellar oscillation problem
from a quadratic self-adjoint eigenvalue problem to a **cubic non-Hermitian**
spectral problem. Its discriminant geometry is organized by a family of open
cubic normal forms whose exceptional-point (EP2) curves terminate at an
**A₃ (cusp) organizing center located strictly in the physically admissible
transport half-plane** (De\* > 0, Λ\* > 0) — *not* at the Navier–Stokes boundary.
In the closed limit (γ = 0) the cusp sits at De\* = 1/(3√3) ≈ 0.19245,
Λ\* = 8/(3√3) ≈ 1.53960; under an open-system radiative rate γ it migrates
according to the exact law (1+2De·γ)³ = 27De²(1+γ²).

The open cubic normal form, its EP2 curves, its A₃ cusp, and the migration law
are realized **exactly** as the sector decomposition of the quasinormal spectrum
of a radiating wave–relaxation system with outgoing Siegert boundary conditions,
in which every γₘ = (m+½)/k is an output of the boundary-value problem.

The manuscript further derives a criterion for whether dense matter can reach
this structure at all. The ratio of the two transport parameters is fixed by the
equation of state alone, Λ/De = (v_FR/v_β)² − 1, with v_FR and v_β the frozen and
β-equilibrated sound speeds, so an exceptional point requires v_FR ≥ 3 v_β.
Urca, hyperonic and colour-superconducting channels fall short of this by one to
two and a half orders of magnitude. **The astrophysical conclusion of the paper
is therefore negative**; the transferable result is the criterion itself.

## Contents

### `code/`
- **`cqg_core.py`** — single source of truth. The Open Cusp Polynomial
  `P(Ω) = -i·De·Ω³ + (1+2De·γ)·Ω² + i[De(1+γ²)+Λ+2γ]·Ω - (1+γ²)`, with root
  finders, EP2/cusp locators (bordered Newton), and the migration law.
- **`certify_table1.py`** — reproduces Table 1 to 20 significant digits
  (EP2 and A₃ cusp coordinates, γ = 0 and γ = 0.1).
- **`hyperboloidal_pencil.py`** — Chebyshev collocation of the hyperboloidal
  wave–relaxation pencil (Eq. 38); the sector-decomposition tower and the
  Pöschl–Teller QNM check at ζ₀ = 0.
- **`bordered_newton_robustness.py`** — multiprecision bordered-Newton EP2/EP3
  solvers (Eq. 76) and the η-continuation of §7.4 (robustness window, complex
  escape, family of **defective** inter-sector crossings). All Jacobian entries
  are analytic; no finite differences. Requires mpmath — double precision is
  provably insufficient (bordered-Jacobian condition numbers 10¹²–10¹⁵).
- **`symbolic_discriminant.py`** — exact symbolic verification (sympy) of the
  cubic discriminant, the rank-one Hessian at the cusp and its kernel, the
  migration law, the finite-$\gamma$ EP2 parametrisation, and the sign relation
  $\Delta_y=-\Delta_\Omega$. This is the symbolic computation referred to in the
  Data Availability statement.
- **`audit_appendixB.py`** — the numerical audit reported in Appendix B:
  bordered-Jacobian condition numbers, Newton contraction factors, forward
  accuracy against the **closed-form** sector-0 EP2, its dependence on the
  resolution N and on the working precision d, and the exactness check of the
  collocation on the polynomial sector. Self-contained (mpmath only).
  The closed form used as the reference is obtained by eliminating ζ₀ between
  P₀ = 0 and ∂ₛP₀ = 0, which leaves the cubic 2τs³ + (1+τ)s² − V₀ = 0, i.e.
  6s³ + 103s² − 2525 = 0 at V₀ = 101/4, τ = 3/100.
  **Precision note.** Every constant is built by `set_precision(d)`, which raises
  `mp.dps` *and* rebuilds V₀, k and τ as exact rationals. Do not hoist them back
  to module level: `mpf('0.03')` evaluated before `mp.dps` is raised binds τ at
  15 digits, and the bordered solve then saturates at 16 forward digits for every
  N, every d and any number of Newton iterations. See CHANGELOG.md.
- **`make_all_figures.py`** — regenerates figures 1, 2 and 3 from the core.
  The EP2 locus is parametrised by y = |Ω| and De(y) attains its **maximum** at
  y = √3, which is the A₃ cusp, so the locus has two branches over every
  De < De*. The polyline is ordered by y, not by De: ordering by De interleaves
  the two branches and produces a spurious sawtooth between De ≈ 0.105 and
  De ≈ 0.19. Figure 2 shades the parameter plane by the exact sign of the
  discriminant rather than by a threshold on Im(root).
- **`diagnostic_dense_continuation.py`** — diagnostic only, not a figure
  generator. Dense η-continuation at N = 20, 24, 28 used to establish that the
  rapid turn of the branch near η ≈ 0.11 is a genuine feature (the eigenvector
  overlap across it is 0.9997, indistinguishable from neighbouring steps) and
  that continuation in η stalls there. It does **not** estimate η_c: fold
  extrapolation from an η-parametrised run is unstable, because some
  resolutions stop at the turn rather than at the fold. The published value
  η_c = 0.1941(1) comes from `bordered_newton_robustness.py`.
- **`make_fig4_robustness.py`** — figure 4. Draws the archived continuation
  data rather than recomputing, so that what is plotted is exactly what is
  tabulated. **Run with no flag to reproduce the manuscript figure**: it draws
  every object the caption refers to — A, A′, B, R2, the complex escape and the
  cusp trajectory. The `--published` flag is a historical leftover that draws a
  reduced variant (branch A plus the three birthplace markers); that variant is
  *not* the figure in the paper and is kept only for reference.
- **`branches_Ap_B.py`** — the two continuation branches born on the defective
  inter-sector crossings that section 7.4 refers to and that were missing from
  earlier versions of this repository: **A′** from the (2,10) crossing (the
  annihilation partner of the fundamental EP2) and **B** from (0,12). A′
  terminates at η = 0.19414 against η = 0.19403 for branch A: the two
  continuations stall on opposite sides of the same fold, which is how
  η_c = 0.1941(1) is bracketed. Writes `data/branches_Ap_B.json`.
- **`eigenvector_overlap.py`** — the branch-identity test of section 7.4(i).
  Computes the normalised eigenvector overlap between consecutive η on branch
  A. Across the turn (η = 0.105 → 0.13) it gives 0.9997, inside the
  0.9987–0.99999 obtained between neighbouring points elsewhere on the branch
  (the archived extremes are 0.998707 and 0.999989).
  Writes `data/eigenvector_overlap.json`.
- **`make_fig5_postmerger.py`** — figure 5 (sweep through two real-axis EP2
  coalescences at the causal fiducial point: f₀ = 3.0 kHz, k_eff = π/R,
  De = 0.17, γ_GW = 0.03; v_Π² + c_s² < 1 at both degeneracies). Produces
  Λ* = 1.565195 and 1.633207, A₃ cusp Λ* = 1.487038, ζ ≈ 2.55 and
  2.67 × 10³⁰ g cm⁻¹ s⁻¹.
- **`projection_coefficient.wl`** — **superseded**. Computes C_proj with
  k_eff ∼ 1/R and f₀ = 3.5 kHz, the identification that section 9.5 now
  rejects as acausal. Retained only so that the provenance of the earlier
  figure is auditable; set `CAUSAL = True` inside it for the current
  construction. Do not use it to regenerate figure 5.

### `figures/`
Generated PDF figures (output of the scripts above).

### `data/`
Certified data sets. Numerical values are stored **as strings**, since
multiprecision quantities cannot be round-tripped through IEEE doubles without
losing everything past the 17th digit.

| File | Contents | Produced by |
|---|---|---|
| `table1_certified.json` | Table 1 coordinates | `certify_table1.py` |
| `defective_crossings.json` | the family of defective inter-sector crossings, Eq. (51) | closed form |
| `robustness_continuation.json` | branch A, the complex escape, R2 and the cusp trajectory | `bordered_newton_robustness.py`; this is its output `fase2_robustez.json` renamed, with the branch values kept at the working precision of the run (24–26 decimals) |
| `branches_Ap_B.json` | branches A′ and B | `branches_Ap_B.py` |
| `eigenvector_overlap.json` | branch-identity test across the turn | `eigenvector_overlap.py` |

**On precision.** The JSON files store values at the working precision of the run
that produced them, as decimal strings — multiprecision quantities cannot be
round-tripped through IEEE doubles. `robustness_continuation.json` therefore
carries 24–26 decimals per entry, and `branches_Ap_B.json` carries the step and
residual data at the precision at which they were computed.

The digits **quoted in the manuscript** are far fewer, and deliberately so: the
bordered-Jacobian condition number is 10¹³–10¹⁵, so the forward accuracy of a
continuation point is roughly 7 significant figures regardless of how many are
stored. The archive is not truncated to that accuracy on purpose — keeping the
full strings lets a reader check the rounding rather than take it on trust.
Do not read the archived digits as certified: only the leading ones are.

## Requirements

```
numpy
scipy
matplotlib
mpmath      # ESSENTIAL: all certified quantities use 25–30 digit arithmetic
sympy       # for the symbolic verification of the sector-decomposition theorem
```

## Reproducing the key results

```bash
cd code
python certify_table1.py             # Table 1 to 20 digits
python symbolic_discriminant.py      # discriminant, Hessian, migration law (exact)
python audit_appendixB.py            # Appendix B audit (conditioning, digit stability)
python make_all_figures.py           # figures 1, 2, 3
python branches_Ap_B.py              # branches A' and B  (~13 min)
python make_fig4_robustness.py       # figure 4 as published (no flag)
python make_fig5_postmerger.py       # figure 5
python eigenvector_overlap.py        # branch-identity test of 7.4(i)

```

## Certified numbers (spot check)

Twenty significant digits, resolution-independent across N = 24, 32, 40 at
working precision d = 30, with attained residuals ≤ 10⁻²⁶.

| Quantity | Value |
|----------|-------|
| A₃ cusp γ=0: De\* | 0.19245008972987525484 = 1/(3√3) |
| A₃ cusp γ=0: Λ\* | 1.5396007178390020387 = 8/(3√3) |
| A₃ cusp γ=0.1: De\*, Λ\*, y\* | 0.20329171841712946742, 1.3703990989252106724, 1.7063465116142830416 |
| EP2 (De=0.15, γ=0.1): Λ\*, y\* | 1.5095752865775080251, 1.2385130101588585499 |

Defective inter-sector crossings (η = 0), from Eq. (51):

| Sectors (m,M) | s | ζ₀ |
|---|---|---|
| (2,10) | −13/2 | 6601/1300 = 5.0776923076923077 |
| (0,12) | −13/2 | 7.5546153846153846 |
| (0,16) | −17/2 | 7.8005882352941176 |

The (2,10) crossing is the η → 0 limit of the branch that annihilates with the
fundamental EP2 at η_c ≈ 0.1941. This is verified by continuation in both
directions and the data are archived in `data/branches_Ap_B.json`: starting
from the crossing, branch A′ rises to ζ₀* = 5.8316 and stalls at η = 0.19414,
while branch A falls to ζ₀* = 5.8749 and stalls at η = 0.19403. The two
continuations exhaust their step floor on opposite sides of the same fold.

The parentage of branch **R2** is *not* established. Two incompatible
candidates have been proposed — the (0,16) crossing at (−8.5, 7.8006) and the
sector-5 exceptional point at (−8.0872, 2.9682) — but the archived continuation
runs over η ∈ [0.05, 0.25] and ζ₀* along it is **not monotone** — 4.993, 4.791,
4.854, 5.407, 6.128 at η = 0.05, 0.10, 0.15, 0.20, 0.25 — so extrapolating to
η → 0 supports neither candidate. The branch should not be attributed to a
parent until it is continued to η = 0. The manuscript states this explicitly in
section 7.4(iii) and the figure 4 caption declines to attribute it.

## Figure-to-file map

| Manuscript | File | Generated by |
|-----------|------|------|
| Figure 1 | figures/fig1_ep_curve.pdf | `make_all_figures.py` |
| Figure 2 | figures/fig2_cusp_zoom.pdf | `make_all_figures.py` |
| Figure 3 | figures/fig6_certified_scaling.pdf | `make_all_figures.py` |
| Figure 4 | figures/fig5_robustness_branches_full.pdf | `make_fig4_robustness.py` (no flag) |
| (reduced variant, not used in the manuscript) | figures/fig5_robustness_branches.pdf | `make_fig4_robustness.py --published` |
| Figure 5 | figures/fig4_physical_NS_application.pdf | `make_fig5_postmerger.py` |

(The filenames predate a reindexing of the manuscript figures and are kept for
continuity of the data sets.)

## License
MIT.
