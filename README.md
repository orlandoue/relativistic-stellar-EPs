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
- **`audit_appendixB.py`** — the numerical audit reported in Appendix B:
  bordered-Jacobian condition numbers, Newton contraction factors, digit
  stability against resolution N and working precision d, and the exactness
  check of the collocation on the polynomial sector. Self-contained (mpmath only).
- **`make_all_figures.py`** — regenerates figures 1, 2 and 3 from the core.
- **`make_fig4_robustness.py`** — figure 4 (η-continuation of the fundamental
  EP2, with resolution study and error bars).
- **`make_fig5_postmerger.py`** — figure 5 (sweep through two real-axis EP2
  coalescences at the causal fiducial point: f₀ = 3.0 kHz, k_eff = π/R,
  De = 0.17, γ_GW = 0.03; v_Π² + c_s² < 1 at both degeneracies).

### `figures/`
Generated PDF figures (output of the scripts above).

### `data/`
Certified data sets: Table 1 coordinates, the η-continuation branches of §7.4,
and the family of defective inter-sector crossings. Numerical values are stored
**as strings**, since multiprecision quantities cannot be round-tripped through
IEEE doubles without losing everything past the 17th digit.

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
python audit_appendixB.py            # Appendix B audit (conditioning, digit stability)
python make_all_figures.py           # figures 1, 2, 3
python make_fig4_robustness.py       # figure 4
python make_fig5_postmerger.py       # figure 5
python bordered_newton_robustness.py # §7.4 continuation (slow, multiprecision)
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

Defective inter-sector crossings (η = 0), from Eq. (49):

| Sectors (m,M) | s | ζ₀ |
|---|---|---|
| (2,10) | −13/2 | 6601/1300 = 5.0776923076923077 |
| (0,12) | −13/2 | 7.5546153846153846 |
| (0,16) | −17/2 | 7.8005882352941176 |

The (2,10) crossing is the η → 0 limit of the branch that annihilates with the
fundamental EP2 at η_c ≈ 0.1941 (verified by continuation: the partner branch
lands on s = −6.4999999, ζ₀ = 5.0776923).

## Figure-to-file map

| Manuscript | File |
|-----------|------|
| Figure 1 | figures/fig1_ep_curve.pdf |
| Figure 2 | figures/fig2_cusp_zoom.pdf |
| Figure 3 | figures/fig6_certified_scaling.pdf |
| Figure 4 | figures/fig5_robustness_branches.pdf |
| Figure 5 | figures/fig4_physical_NS_application.pdf |

(The filenames predate a reindexing of the manuscript figures and are kept for
continuity of the data sets.)

## License
MIT.
