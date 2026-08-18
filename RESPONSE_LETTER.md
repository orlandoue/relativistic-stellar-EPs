# Response to the Referee — CQG-115250.R2 → R3

Dear Referee,

Thank you for a further careful reading. Your three points concern section 7,
and addressing them has changed that section substantially and, I think, for the
better: what was previously asserted about the role of the boundary conditions is
now proved, and the scope of the construction is stated explicitly rather than
left implicit. Below I answer each point in turn and then list two further
corrections I made on my own initiative, one of which is a numerical error in the
appendix that worked against the paper.

Changed passages are marked in the revised manuscript.

---

## (1) The monomials σᵐ are basis elements, not eigenfunctions; a polynomial ansatz with lower-order terms and the recurrence should be given; and it should be explained why the lower-order equations do not obstruct the construction

You are right, and the previous formulation invited exactly the confusion you
describe. Section 7.2 has been restructured so that triangularity and the
eigenvalue problem are kept apart from the start.

- **Lemma 1** now states only the operator identities
  L₂σᵐ = −(2m+1)σᵐ and L₁σᵐ = m(m−1)σᵐ⁻² − [m(m+1)+V₀]σᵐ, from which T(s) is
  triangular in the monomial basis with the sector cubics on the diagonal. It
  makes no claim about eigenfunctions.
- Immediately after it, the manuscript now says in as many words that
  *"triangularity alone does not identify σᵐ as an eigenfunction of (38): the
  eigenvalue problem is a statement about a specific linear combination of
  monomials, not about the diagonal entries in isolation."*
- A **polynomial ansatz containing all lower-order terms** is then introduced,
  φ(σ) = Σₖ cₖ σ^{k₀+2k} with c₀ ≠ 0, and substitution gives the **two-term
  recurrence**, Eq. (42), which is displayed and used.
- **Lemma 2** answers your question about the lower-order equations directly.
  The point is that when Pₘ(s) = 0 the series closes and *every* order is
  satisfied, not merely the orders below the truncation: for every j > m both
  sides of (42) vanish identically, irrespective of whether Pⱼ(s) itself
  vanishes there. The converse is also stated. Closure at the bottom of each
  parity sector follows because the lowering coefficient m(m−1) vanishes at
  m = 0 and m = 1.

So the solution is constructed by satisfying the equations at all orders, not by
requiring a determinant to vanish. The determinant argument now appears only in
section 7.4(iii), where it is used solely to count algebraic multiplicity on a
finite-dimensional invariant subspace, with an explicit disclaimer to that
effect.

## (2) The truncation is imposed from the outset; the boundary conditions enter there, not through the hyperboloidal coordinate; the restriction is model-specific; and the construction does not exhaust the general solution space

This was the most valuable of your comments, because the previous version did
assume what it should have proved. All four sub-points are now addressed.

**(2a) Truncation is proved, not assumed.** A new Proposition
(*Regularity ⇔ polynomial*) shows that truncation is **necessary**, not merely
sufficient. For s ≠ −1/τ the equation is Fuchsian with regular singular points
{−1, +1, ∞}; the Fuchs relation is verified; the origin is an ordinary point, so
the two series with k₀ = 0, 1 span the full local solution space; the indicial
equation at the finite singular points is r(r+s) = 0. Imposing the outgoing
condition at both endpoints makes the solution entire, and the indicial equation
at infinity, ρ² + (2s+1)ρ + K = 0, has ρ = m as a root **if and only if**
Pₘ(s) = 0. Liouville's theorem then forces a polynomial of degree m. The
polynomial ansatz therefore entails no loss of generality.

**(2b) Where the boundary conditions enter.** Stated explicitly, in your terms:
the role of the boundary conditions is *"encoded through endpoint regularity
acting on the series (41), and **not** through the introduction of the
hyperboloidal coordinate itself, which merely supplies the frame in which that
regularity becomes equivalent to polynomial truncation."*

The revision also makes precise **which** branch regularity selects, and why the
usual justification is inadequate. Writing ψ = e^{−sh}φ with h = ln cosh x and
u = 1 ∓ σ, the two Frobenius branches behave as ψ ∼ e^{−s|x|} (outgoing) and
ψ ∼ e^{+s|x|} (incoming). Excluding u^{−s} is therefore exactly the imposition of
purely outgoing behaviour — and it **cannot** be justified by boundedness, since
for quasinormal s one has Re s < 0 and the incoming branch decays. What excludes
it is that it is incoming, not that it grows. Phrasing the criterion this way
also keeps it meaningful at the resonant values −s ∈ ℤ₊, where a criterion
phrased as analyticity would select the wrong branch.

**(2c) Model specificity.** Conceded explicitly. The reduction to a single
low-degree algebraic condition follows from the **two-term** character of the
recurrence, itself a consequence of L₁ lowering the degree by exactly two — a
property of the Pöschl–Teller potential, equivalently of the underlying
hypergeometric structure. The manuscript now states that we *"do not claim it
survives for potentials or coordinate choices whose associated recurrence is not
two-term."*

**(2d) The general solution space.** Also addressed directly. For s not a root of
any Pₘ the recurrence never terminates and, since the ratio tends to one, the
series has radius of convergence exactly one and realises the u^{−s} branch at
the endpoints. What the Proposition establishes is that this branch is excluded
by regularity: **the polynomial family exhausts the eigenfunctions, not the
solutions.**

**A gap you did not raise, which I found while writing (2a) and have closed.**
The step "free of branch points ⇒ entire" fails at the isolated resonant values
−s = n ∈ ℤ₊, where the outgoing Frobenius solution is the one of smaller exponent
and may carry a logarithm. The revision settles this by an elementary pairing:
the condition Pₘ(−n) = 0 fixes the coupling to Eq. (44), an expression invariant
under m ↦ m′ = 2n−1−m. Since m + m′ = 2n−1 is odd, m and m′ lie in different
parity sectors, so for 0 ≤ m ≤ 2n−1 the two polynomial solutions are independent
and span the whole (second-order) solution space: no logarithmic solution exists
to be excluded. For m ≥ 2n the companion index is negative; there the degree-m
solution factorises as (1−σ²)ⁿ q(σ), vanishing to order exactly n at both
endpoints, which is the **incoming** exponent — so it is not a quasinormal mode,
and the conclusion cannot be asserted. This is stated rather than papered over.

Consistently with that, the **Theorem is now stated for couplings outside the
discrete resonant set**, matching the hypotheses of the Proposition it rests on,
and the manuscript separates the two inclusions: no spectrum outside the union
holds for every ζ₀ (it uses only Lemma 1), while the reverse inclusion is the one
that needs the Proposition. The nearest resonance to the certified EP2 is
quantified: the EP2 sits at (s, ζ₀) = (−6.1926, 7.5479); the nearest resonance
with m ≥ 2n is (n,m) = (6,12) at (−6, 9.191), a distance 1.65 in that plane. No
result in the paper is evaluated at a resonant coupling.

The word **"complete"** has likewise been made a theorem rather than an
expectation, with an explicit statement of what is *not* claimed: the manuscript
now says that we do not assert that the essential spectrum of a closed extension
to L²(−1,1) or to an energy norm is empty, citing your own community's results on
the norm-dependence of hyperboloidal quasinormal-mode spectra.

## (3) References on hyperboloidal coordinates

Corrected. Zenginoğlu, *Class. Quantum Grav.* **25**, 145002 (2008) and
*Phys. Rev. D* **83**, 127502 (2011) are now cited as the foundational works, and
they appear **before** the PRX/PRD papers at both places where the hyperboloidal
reformulation is introduced (sections 7.1 and 7.2).

---

## Corrections made on my own initiative

**A numerical error in Appendix B, which understated the paper's own result.**
While preparing this revision I re-audited the certification code and found that
`audit_appendixB.py` built its constants at module level, before `mp.dps` was
raised. `mpf('0.03')` was therefore parsed at 15 digits, so the script silently
solved the pencil for a relaxation time differing from 3/100 in the 17th
significant digit. The bordered solve consequently saturated at 16 forward digits
regardless of resolution, working precision or iteration count, and R2 reported
that saturation as a conditioning ceiling.

With the constants built as exact rationals, the operator-level solve reproduces
**all twenty** published digits of Table 1 at N = 24, d = 30, and the corrected
picture is stronger than the one it replaces. Since the collocation is exact on
the polynomial sector, the discretisation error is identically zero and the whole
deviation is arithmetic; the forward accuracy therefore *degrades* with
resolution — about 21, 20 and 18 digits at N = 24, 32, 40 — tracking
κ₂ = 5.45×10¹², 3.34×10¹⁴, 7.67×10¹⁵, and the a-priori bound κ₂‖F‖/‖J‖ is tight
to about one digit at every resolution. At fixed N the accuracy tracks the
working precision one-for-one: 16, 21, 31 and 41 digits at d = 25, 30, 40, 50.
Appendix B now reports this, together with the run-to-run fluctuation of the digit
counts, which is one to two units because the terminal residual sits at the
roundoff floor. This is corroborated by `hyperboloidal_pencil.py`, which was never
affected and already returned 21–22 correct digits. The correction is documented
in the repository CHANGELOG; Table 1, the conclusions and the plain-language
summary have been updated accordingly.

**A bridge between sections 7 and 9.** The transition from the exactly solvable
model to the microphysical criterion involved a change of object that was not
declared. Section 9.5 now opens by stating it: the (1+1) Pöschl–Teller system
establishes that the open normal form, its EP2 curves, its cusp and the migration
law are realized by genuine Siegert eigenvalues; what follows inherits that
guarantee but not the specific values of γ, V₀ or k used there.

**Smaller items.** Branch R2 is no longer attributed to a parent in the figure 4
caption, matching the body text; the identification of branch A′ with the (2,10)
crossing at η = 0 is now stated with its actual accuracy (eight digits) and the
reason for it (the bordered EP2 system is singular at a defective inter-sector
crossing, so that endpoint is approached rather than solved for); the persistence
range of the cusp now records that the continuation reaches η ≈ 0.44 while only
η ≤ 0.20 is archived at the certification standard; the symbol γ is no longer
overloaded in section 9.6; the criterion v_FR ≥ 3v_β is stated in the opening of
the abstract rather than only at its end; and the repository README, equation
cross-references, nomenclature and figure-to-file map have been corrected.

---

I hope you find the revised section 7 clearer and the scope of its claims properly
delimited. I am grateful for the three rounds of scrutiny; the section is
considerably more honest than it was.

Yours sincerely,

Orlando Urbina-Gonzalez
