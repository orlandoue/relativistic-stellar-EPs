# Are the inter-sector crossings defective or diabolic?

Supporting computation for section 7.4(iii) and Eq. (51) of *Exceptional points
and cusp singularities in relativistic stellar oscillations with causal bulk
relaxation*.

## The question

Under the transport-profile deformation ζ(x) = ζ₀ sech²x (1 + η sech²x), the
manuscript states that at η = 0 any two sectors of equal parity m, M share an
eigenvalue at

```
s = -(m+M+1)/2 ,   zeta0 = (1+s*tau)[s^2 + (2m+1)s + m(m+1) + V0] / (-s)
```

and that these crossings are **defective** — algebraically double but
geometrically simple, i.e. Jordan blocks — rather than **diabolic** (semisimple)
in the von Neumann–Wigner / Berry sense. Three members of the family are used as
η = 0 birthplaces for the EP2 branches of figure 4: (2,10), (0,12) and (0,16).

The distinction is not cosmetic. The pencil is *triangular*, not diagonal, in the
monomial basis σᵐ (Theorem, section 7.2). Triangularity fixes the spectrum, since
the characteristic determinant of the restriction to an invariant subspace is
∏ₘ Pₘ(s), but it says nothing by itself about eigenvectors: the off-diagonal
coupling could leave two independent null directions at a crossing, or collapse
them to one. Section 7.4(iii) argues from the recurrence that only one survives.
This directory checks that argument by direct computation in exact arithmetic,
and independently at the level of the full discretised operator.

## Result

All three points are defective. Algebraic multiplicity 2 by construction — two
diagonal factors vanish together — but **geometric multiplicity 1**.

| Point  | s*   | zeta0*                  | alg. mult. | geom. mult. | eigenvector support |
|--------|------|-------------------------|-----------|------------|---------------------|
| (2,10) | −6.5 | 6601/1300 ≈ 5.077692    | 2 | **1** | σ⁰, σ² (sector m = 2) |
| (0,12) | −6.5 | 9821/1300 ≈ 7.554615    | 2 | **1** | σ⁰ (sector m = 0)     |
| (0,16) | −8.5 | 13261/1700 ≈ 7.800588   | 2 | **1** | σ⁰ (sector m = 0)     |

## Mechanism

Script `01` derives symbolically that

```
T(s) sigma^m = P_m(s) sigma^m  -  (1 + s*tau) m(m-1) sigma^(m-2)
```

so in the basis {1, σ², σ⁴, …} the operator is upper bidiagonal: the diagonal
carries the sector cubics Pₘ(s), and the off-diagonal coefficient
−(1+sτ)m(m−1) is non-zero for every even m ≥ 2. That unbroken coupling chain is
what decides the eigenvector structure.

Solving T(s\*)v = 0 by descending the chain from the top: the zero of the
**higher** degree (M = 10, 12 or 16) never anchors a free direction, because the
chain descending from it forces a consistency condition at the rung of the
**lower** zero (m = 2 or 0) that is satisfied only if the whole direction
vanishes. The surviving eigenvector lives entirely in the lower sector. The
higher-sector zero contributes to the algebraic multiplicity but not to the
geometric one.

This is not an artefact of truncating the basis: script `02` repeats the exact
rank computation in rational arithmetic at Mmax = 20, 30, 40, 60 and obtains
nullity 1 in all four cases for all three points. Script `03` repeats the check
on the full Chebyshev-discretised operator, built exactly as
`bordered_newton_robustness.py` builds it: the second-smallest singular value
sits 10⁶–10⁷ above the double-precision floor, confirming independently that
there is one null direction, not two.

## Why this supports the physics of section 7.4

If the crossings were genuinely semisimple, a generic perturbation η > 0 would
open them **linearly** in η, and there would be no reason to expect EP2 branches
with square-root splitting to emerge from them. The manuscript does observe
genuine EP2 branches (A′, B) being born exactly at these points, which is what a
defective starting point predicts. This computation and the continuation results
of section 7.4 are therefore two independent routes to the same conclusion.

None of the quantitative results is affected: η_c = 0.1941(1), the identity of
branch A′ (born at (−6.5, 6601/1300)), the complex escape at η = 0.21 and
η = 0.40, and the certified splitting exponents all stand.

## Running the scripts

```bash
pip install sympy mpmath numpy
python3 01_estructura_triangular_simbolica.py
python3 02_multiplicidad_geometrica_racional.py
python3 03_verificacion_operador_chebyshev.py     # reuses ../bordered_newton_robustness.py
```

- **`01_estructura_triangular_simbolica.py`** — symbolic derivation (sympy, no
  rounding) of the action of T(s) on σᵐ, isolating the off-diagonal coupling
  coefficient. Everything else rests on this.
- **`02_multiplicidad_geometrica_racional.py`** — the central computation: exact
  rank and null space of T(s\*, ζ₀\*) in rational arithmetic at the three declared
  points, with a truncation-stability check and explicit extraction of the
  surviving eigenvector.
- **`03_verificacion_operador_chebyshev.py`** — independent confirmation on the
  full discretised operator rather than on the monomial matrix.
