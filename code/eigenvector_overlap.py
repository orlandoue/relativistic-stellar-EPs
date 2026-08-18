"""
eigenvector_overlap.py -- the branch-identity test quoted in section 7.4(i).

Section 7.4(i) states that the same analytic branch is recovered across the
rapid turn near eta ~ 0.11, "verified by eigenvector overlap:
|<phi_a,phi_b>|/(||phi_a|| ||phi_b||) = 0.9997 across the turn
(eta = 0.105 -> 0.13), indistinguishable from the 0.9991-0.9998 obtained
between neighbouring points on either side."

That number appeared only as a remark in a docstring: no script computed it and
no data set contained it.  This script computes it, prints the full table of
neighbouring overlaps, and archives the result.

The eigenvector of the bordered system is the null vector phi of T(s*), i.e.
the eigenfunction at the exceptional point.  Each phi is normalised to unit
maximum entry, as in the continuation; the overlap is the modulus of the
Hermitian inner product of the normalised vectors evaluated at the Chebyshev
nodes.  Because the nodes are shared between runs at equal N, this is a fair
comparison of the same object.

Output: ../data/eigenvector_overlap.json
Usage:  python eigenvector_overlap.py
"""
from mpmath import mp, mpf, mpc, matrix, lu_solve, cos, pi, nstr, sqrt, conj
import json

mp.dps = 25
V0 = mpf('25.25'); TAU = mpf('0.03'); N = 28

# eta grid: the turn is crossed with a coarse step, and neighbouring pairs on
# either side of it provide the control values.
ETAS = ['0.00', '0.05', '0.08', '0.10', '0.105', '0.13', '0.15', '0.175', '0.19']
TURN = ('0.105', '0.13')


def setup(N, eta):
    x = [cos(pi * mpf(j) / N) for j in range(N + 1)]
    c = [(2 if j in (0, N) else 1) * (-1) ** j for j in range(N + 1)]
    D = matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i, j] = mpf(c[i]) / mpf(c[j]) / (x[i] - x[j])
    for i in range(N + 1):
        D[i, i] = -sum(D[i, j] for j in range(N + 1) if j != i)
    D2 = D * D; n = N + 1
    L1 = matrix(n, n); L2 = matrix(n, n)
    for i in range(n):
        w = 1 - x[i] ** 2
        for j in range(n):
            L1[i, j] = w * D2[i, j] - 2 * x[i] * D[i, j]
            L2[i, j] = -2 * x[i] * D[i, j]
        L1[i, i] -= V0; L2[i, i] -= 1
    Cp = [1 + eta * (1 - x[i] ** 2) for i in range(n)]
    return L1, L2, Cp, n


def T_all(s, z0, tau, L1, L2, Cp, n):
    T = matrix(n, n); T1 = matrix(n, n); T2 = matrix(n, n)
    a = 1 + s * tau
    for i in range(n):
        for j in range(n):
            base = -s * L2[i, j] - L1[i, j]
            T[i, j] = a * base
            T1[i, j] = tau * base - a * L2[i, j]
            T2[i, j] = -2 * tau * L2[i, j]
        T[i, i] += a * s * s + z0 * s * Cp[i]
        T1[i, i] += tau * s * s + 2 * a * s + z0 * Cp[i]
        T2[i, i] += 4 * tau * s + 2 * a
    return T, T1, T2


def solve_EP2(s, z0, phi0, L1, L2, Cp, n, itmax=40):
    tol = mpf('10') ** (-(mp.dps - 8))
    phi = phi0.copy(); phi1 = matrix([mpc(0)] * n); c = phi0.copy()
    m = 2 * n + 2
    for it in range(itmax):
        T, T1, T2 = T_all(s, z0, TAU, L1, L2, Cp, n)
        Tphi = T * phi; T1phi = T1 * phi
        Tphi1 = T * phi1; T1phi1 = T1 * phi1; T2phi = T2 * phi
        F = matrix(m, 1)
        for i in range(n):
            F[i] = Tphi[i]; F[n + i] = T1phi[i] + Tphi1[i]
        F[2 * n] = sum(c[i] * phi[i] for i in range(n)) - 1
        F[2 * n + 1] = sum(c[i] * phi1[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        if res < tol:
            return s, z0, phi, res, True
        J = matrix(m, m)
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
                J[n + i, j] = T1[i, j]; J[n + i, n + j] = T[i, j]
            J[i, 2 * n] = T1phi[i]
            J[i, 2 * n + 1] = s * Cp[i] * phi[i]
            J[n + i, 2 * n] = T2phi[i] + T1phi1[i]
            J[n + i, 2 * n + 1] = Cp[i] * phi[i] + s * Cp[i] * phi1[i]
        for j in range(n):
            J[2 * n, j] = c[j]; J[2 * n + 1, n + j] = c[j]
        d = lu_solve(J, -F)
        for i in range(n):
            phi[i] += d[i]; phi1[i] += d[n + i]
        s += d[2 * n]; z0 += d[2 * n + 1]
    return s, z0, phi, res, False


def renorm(phi, n):
    mx = max(abs(phi[i]) for i in range(n))
    return matrix([phi[i] / mx for i in range(n)])


def overlap(a, b, n):
    num = abs(sum(conj(a[i]) * b[i] for i in range(n)))
    na = sqrt(sum(abs(a[i]) ** 2 for i in range(n)))
    nb = sqrt(sum(abs(b[i]) ** 2 for i in range(n)))
    return num / (na * nb)


if __name__ == "__main__":
    s = mpc('-6.192565050794292749611')
    z = mpc('7.547876432887540125414')
    phi = None
    store = {}
    print(f"{'eta':>8} {'s*':>18} {'zeta0*':>18}  ok")
    for ev in ETAS:
        L1, L2, Cp, n = setup(N, mpf(ev))
        seed = phi if phi is not None else matrix([mpc(1, 0)] * n)
        s, z, ph, res, ok = solve_EP2(s, z, seed, L1, L2, Cp, n)
        ph = renorm(ph, n)
        store[ev] = ph
        phi = ph
        print(f"{ev:>8} {nstr(s,12):>18} {nstr(z,12):>18}  {ok}")

    print("\nnormalised eigenvector overlaps between consecutive eta:")
    rows = []
    for a, b in zip(ETAS[:-1], ETAS[1:]):
        o = overlap(store[a], store[b], N + 1)
        flag = "   <-- ACROSS THE TURN" if (a, b) == TURN else ""
        print(f"  {a:>6} -> {b:<6}  {nstr(o, 6)}{flag}")
        rows.append([a, b, nstr(o, 6), (a, b) == TURN])

    off = [float(r[2]) for r in rows if not r[3]]
    across = [float(r[2]) for r in rows if r[3]]
    print(f"\n  across the turn : {across[0]:.4f}")
    print(f"  neighbours off the turn : {min(off):.4f} - {max(off):.4f}")
    verdict = ("indistinguishable" if min(off) <= across[0] <= max(off)
               else "OUTSIDE the neighbouring range -- re-examine")
    print(f"  verdict: {verdict}")

    with open("../data/eigenvector_overlap.json", "w") as f:
        json.dump({
            "_note": ("Branch-identity test of section 7.4(i). Normalised "
                      "eigenvector overlap |<phi_a,phi_b>|/(||phi_a|| ||phi_b||) "
                      "between consecutive eta on branch A, N=28, 25 digits, "
                      "tau=0.03. The pair (0.105, 0.13) straddles the rapid "
                      "turn near eta ~ 0.11."),
            "overlaps": rows,
            "across_the_turn": across[0],
            "neighbour_range": [min(off), max(off)],
            "verdict": verdict}, f, indent=1)
    print("\nWrote ../data/eigenvector_overlap.json")
