"""
branches_Ap_B.py -- the two continuation branches that section 7.4 declares
"tabulated in the repository" but which were absent from the archived data set.

  * Branch A'  : born at eta = 0 on the defective inter-sector crossing of
                 sectors (m,M) = (2,10),  (s, zeta0) = (-13/2, 6601/1300).
                 This is the annihilation partner of the fundamental EP2
                 (branch A) at eta_c.
  * Branch B   : born at eta = 0 on the defective inter-sector crossing of
                 sectors (0,12),  (s, zeta0) = (-13/2, 9821/1300).

Both are continued upward in eta with the same bordered-Newton machinery and
the same step control documented in Appendix B.  Output:
../data/branches_Ap_B.json   (values as strings, truncated to 7 significant
figures, for the same reason as robustness_continuation.json: the bordered
Jacobian condition number is 1e13-1e15, so more digits would be fiction).

Usage:  python branches_Ap_B.py
"""
from mpmath import mp, mpf, mpc, matrix, lu_solve, cos, pi, nstr
import json, time, sys

mp.dps = 25
V0 = mpf('25.25')
TAU = mpf('0.03')          # De = 0.15
N = 28


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
    D2 = D * D
    n = N + 1
    L1 = matrix(n, n); L2 = matrix(n, n)
    for i in range(n):
        w = 1 - x[i] ** 2
        for j in range(n):
            L1[i, j] = w * D2[i, j] - 2 * x[i] * D[i, j]
            L2[i, j] = -2 * x[i] * D[i, j]
        L1[i, i] -= V0
        L2[i, i] -= 1
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


def solve_EP2(s, z0, tau, phi0, L1, L2, Cp, n, itmax=40, tol=None):
    """Bordered system (76), 2n+2 complex unknowns (phi, phi1, s, zeta0)."""
    tol = tol or mpf('10') ** (-(mp.dps - 8))
    phi = phi0.copy(); phi1 = matrix([mpc(0)] * n); c = phi0.copy()
    m = 2 * n + 2
    res = mpf('inf')
    for it in range(itmax):
        T, T1, T2 = T_all(s, z0, tau, L1, L2, Cp, n)
        Tphi = T * phi; T1phi = T1 * phi
        Tphi1 = T * phi1; T1phi1 = T1 * phi1; T2phi = T2 * phi
        F = matrix(m, 1)
        for i in range(n):
            F[i] = Tphi[i]
            F[n + i] = T1phi[i] + Tphi1[i]
        F[2 * n] = sum(c[i] * phi[i] for i in range(n)) - 1
        F[2 * n + 1] = sum(c[i] * phi1[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        if res < tol:
            return s, z0, phi, res, True
        J = matrix(m, m)
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
                J[n + i, j] = T1[i, j]
                J[n + i, n + j] = T[i, j]
            J[i, 2 * n] = T1phi[i]
            J[i, 2 * n + 1] = s * Cp[i] * phi[i]
            J[n + i, 2 * n] = T2phi[i] + T1phi1[i]
            J[n + i, 2 * n + 1] = Cp[i] * phi[i] + s * Cp[i] * phi1[i]
        for j in range(n):
            J[2 * n, j] = c[j]
            J[2 * n + 1, n + j] = c[j]
        d = lu_solve(J, -F)
        for i in range(n):
            phi[i] += d[i]; phi1[i] += d[n + i]
        s += d[2 * n]; z0 += d[2 * n + 1]
    return s, z0, phi, res, False


def renorm(phi, n):
    mx = max(abs(phi[i]) for i in range(n))
    return matrix([phi[i] / mx for i in range(n)])


def crossing(m_, M_):
    """Exact birthplace, Eq. (51)."""
    s = -mpf(m_ + M_ + 1) / 2
    z = (1 + s * TAU) * (s ** 2 + (2 * m_ + 1) * s + m_ * (m_ + 1) + V0) / (-s)
    return s, z


def continue_branch(label, m_, M_, eta_max, h0=mpf('0.02'), hmin=mpf('2e-5'),
                    budget=2400):
    s0, z0 = crossing(m_, M_)
    print(f"\n[{label}] birthplace from Eq.(51), sectors ({m_},{M_}): "
          f"s = {nstr(s0, 10)}, zeta0 = {nstr(z0, 12)}", flush=True)
    s, z = mpc(s0), mpc(z0)
    L1, L2, Cp, n = setup(N, mpf(0))
    phi = matrix([mpc(1, 0)] * n)
    s, z, phi, res, ok = solve_EP2(s, z, TAU, phi, L1, L2, Cp, n)
    print(f"   eta=0.00000 : s*={nstr(s,12)} zeta0*={nstr(z,12)} "
          f"res={nstr(res,2)} ok={ok}", flush=True)
    rows = [(0.0, str(s), str(z), str(res))]
    phi = renorm(phi, n)
    eta = mpf(0); h = h0
    t0 = time.time()
    while eta < eta_max and h > hmin and time.time() - t0 < budget:
        et = eta + h
        L1, L2, Cp, n = setup(N, et)
        try:
            st, zt, pt, res, ok = solve_EP2(s, z, TAU, phi, L1, L2, Cp, n)
        except ZeroDivisionError:
            ok = False
        if ok and abs(st - s) < mpf('0.25') + 3 * h:
            eta, s, z, phi = et, st, zt, renorm(pt, n)
            rows.append((float(eta), str(s), str(z), str(res)))
            print(f"   eta={float(eta):.5f} : s*={nstr(s,12)} "
                  f"zeta0*={nstr(z,12)} res={nstr(res,2)}", flush=True)
            h = min(h * mpf('1.4'), h0)
        else:
            h = h / 2
    print(f"   [{label}] stopped at eta={float(eta):.5f} "
          f"(h={float(h):.2e}, {time.time()-t0:.0f}s)", flush=True)
    return rows


def _c(x):
    return complex(str(x).replace('(', '').replace(')', '').replace(' ', ''))


def trunc(rows):
    out = []
    for r in rows:
        out.append([r[0], nstr(mpc(_c(r[1])), 7), nstr(mpc(_c(r[2])), 7),
                    nstr(mpf(float(r[3])), 2)])
    return out


if __name__ == "__main__":
    out = {}
    # Branch A': the annihilation partner of the fundamental EP2.
    def save():
        out["_note"] = (
            "Bordered-Newton continuation of Eq. (78) with the deformed profile "
            "C = 1 + eta(1-sigma^2), N=28, 25-digit arithmetic, tau=0.03 "
            "(De=0.15), V0=25.25. Birthplaces are the exact defective "
            "inter-sector crossings of Eq. (51). Values truncated to 7 "
            "significant figures. Fields: [eta, s*, zeta0*, residual].")
        with open("../data/branches_Ap_B.json", "w") as f:
            json.dump(out, f, indent=1)

    out["branch_Aprime_from_(2,10)"] = trunc(
        continue_branch("A'", 2, 10, mpf('0.25')))
    save()
    # Branch B: followed as far as the budget allows.
    out["branch_B_from_(0,12)"] = trunc(
        continue_branch("B", 0, 12, mpf('0.60'), budget=2000))
    save()
    print("\nWrote ../data/branches_Ap_B.json")
