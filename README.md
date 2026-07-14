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

## Contents

### `code/`
- **`cqg_core.py`** — single source of truth. The Open Cusp Polynomial
  `P(Ω) = -i·De·Ω³ + (1+2De·γ)·Ω² + i[De(1+γ²)+Λ+2γ]·Ω - (1+γ²)`, with root
  finders, EP2/cusp locators (bordered Newton), and the migration law.
- **`certify_table1.py`** — reproduces Table 1 to 20 significant digits in
  40-digit arithmetic (EP2 and A₃ cusp coordinates, γ = 0 and γ = 0.1).
- **`hyperboloidal_pencil.py`** — Chebyshev collocation of the hyperboloidal
  wave–relaxation pencil (Eq. 37); the sector-decomposition tower and the
  Pöschl–Teller QNM check at ζ₀ = 0.
- **`bordered_newton_robustness.py`** — multiprecision bordered-Newton EP2/EP3
  solvers (Eq. 60) and the η-continuation of §7.4 (robustness window, complex
  escape, diabolic family). Requires mpmath (double precision is provably
  insufficient: condition numbers ~10¹⁴).
- **`make_all_figures.py`** — regenerates figures 1, 2, 3, 6 from the core.
- **`make_fig4_postmerger.py`** — figure 4 (post-merger sweep through the two
  exact real-axis EP2 coalescences, De = 0.11, γ = 0.03, Cproj = 5.06×10⁻³²).

### `figures/`
Generated PDF figures (output of the scripts above).

### `data/`
Certified data sets: Table 1 coordinates, the η-continuation branches of §7.4,
and the diabolic-point family, each with residuals at working precision.

## Requirements

```
numpy
scipy
matplotlib
mpmath      # ESSENTIAL: all certified quantities use 25–40 digit arithmetic
sympy       # for the symbolic verification of the sector-decomposition theorem
```

## Reproducing the key results

```bash
cd code
python certify_table1.py            # Table 1 to 20 digits
python make_all_figures.py          # figures 1, 2, 3, 6
python make_fig4_postmerger.py      # figure 4
python bordered_newton_robustness.py  # §7.4 robustness continuation (slow, mp)
```

## Certified numbers (spot check)

| Quantity | Value |
|----------|-------|
| A₃ cusp γ=0: De\* | 0.19245008972987526219 = 1/(3√3) |
| A₃ cusp γ=0: Λ\* | 1.53960071783900209752 = 8/(3√3) |
| A₃ cusp γ=0.1: De\*, Λ\*, y\* | 0.20329171841712948, 1.37039909892521061, 1.70634651161428308 |
| EP2 (De=0.15, γ=0.1): Λ\*, y\* | 1.50957528657750806, 1.23851301015885862 |

## License
MIT (or as required by the journal).
