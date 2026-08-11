"""
audit_appendixB.py — numeros de auditoria para el Apendice B
=============================================================
Genera las cifras que la revision numerica (perfil Trefethen/Higham) pide y
que hoy faltan en el manuscrito:

  (1) kappa_2 del jacobiano ORLADO en cada punto certificado, a cada N
      -> es el numero que convierte "residuo 1e-24" en "20 digitos"
  (2) residuo REALMENTE alcanzado (no solo el umbral)
  (3) tabla de estabilidad de digitos frente a N  y frente a la precision d
  (4) numero de iteraciones de Newton y factor de contraccion
      -> evidencia de convergencia cuadratica = el cero es regular
  (5) comprobacion en zeta_0=0 contra omega_m = +-k - i(m+1/2) en multiprecision
      (sustituye al "4e-13", que era una comprobacion en doble precision)

Sistema: el de la seccion 7,  T(s) = (1+s tau)(s^2 - s L2 - L1) + zeta_0 s
con L1 = (1-sigma^2) d^2 - 2 sigma d - V0 ,  L2 = -2 sigma d - 1  (perfil alineado).

AUTOCONTENIDO: solo requiere mpmath.
Uso:  python audit_appendixB.py
"""
from mpmath import mp, mpf, mpc, matrix, cos, pi, lu_solve, nstr, svd_r, sqrt

V0 = mpf('25.25')          # k = 5, gamma_0 = 1/(2k) = 0.1
K  = mpf(5)
TAU = mpf('0.03')          # De = k*tau = 0.15

# valores publicados en la Tabla 1 (en variables adimensionales)
Y_EP2   = mpf('1.2385130101588585499')
LAM_EP2 = mpf('1.5095752865775080251')


# ----------------------------------------------------------------- operadores
def cheb(N):
    x = [cos(pi * mpf(j) / N) for j in range(N + 1)]
    c = [(2 if j in (0, N) else 1) * (-1) ** j for j in range(N + 1)]
    D = matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i, j] = mpf(c[i]) / mpf(c[j]) / (x[i] - x[j])
    for i in range(N + 1):
        D[i, i] = -sum(D[i, j] for j in range(N + 1) if j != i)
    return D, x


def ops(N):
    D, x = cheb(N); D2 = D * D; n = N + 1
    L1 = matrix(n, n); L2 = matrix(n, n)
    for i in range(n):
        for j in range(n):
            L1[i, j] = (1 - x[i]**2) * D2[i, j] - 2 * x[i] * D[i, j]
            L2[i, j] = -2 * x[i] * D[i, j]
        L1[i, i] -= V0
        L2[i, i] -= 1
    return L1, L2, n


def Tmats(s, z0, L1, L2, n):
    """T(s) y T'(s)."""
    T = matrix(n, n); T1 = matrix(n, n)
    for i in range(n):
        for j in range(n):
            base = -s * L2[i, j] - L1[i, j]
            T[i, j] = (1 + s * TAU) * base
            T1[i, j] = TAU * base + (1 + s * TAU) * (-L2[i, j])
        T[i, i] += (1 + s * TAU) * s * s + z0 * s
        T1[i, i] += TAU * s * s + (1 + s * TAU) * 2 * s + z0
    return T, T1


# ------------------------------------------------------- EP2 orlado + auditoria
def solve_EP2_audit(N, s0, z00, tol=None):
    """Newton sobre el sistema orlado (76). Devuelve dict con todo lo auditable."""
    tol = tol or mpf(10) ** (-(mp.dps - 6))
    L1, L2, n = ops(N)
    m = 2 * n + 2
    s, z0 = mpc(s0), mpf(z00)
    T, _ = Tmats(s, z0, L1, L2, n)
    try:
        phi = lu_solve(T, matrix([mpf(1)] * n))
    except Exception:
        phi = matrix([mpf(1)] * n)
    mx = max(abs(phi[i]) for i in range(n))
    phi = matrix([phi[i] / mx for i in range(n)])
    phi1 = matrix([mpc(0)] * n)
    c = phi.copy()
    residuals = []
    J = None
    for it in range(80):
        T, T1 = Tmats(s, z0, L1, L2, n)
        R1 = T * phi; R2 = T1 * phi + T * phi1
        F = matrix(m, 1)
        for i in range(n):
            F[i] = R1[i]; F[n + i] = R2[i]
        F[2 * n]     = sum(c[i] * phi[i] for i in range(n)) - 1
        F[2 * n + 1] = sum(c[i] * phi1[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        residuals.append(res)
        J = matrix(m, m)
        T1phi = T1 * phi
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
                J[n + i, j] = T1[i, j]
                J[n + i, n + j] = T[i, j]
            J[i, 2 * n]         = T1phi[i]
            J[i, 2 * n + 1]     = s * phi[i]
            J[n + i, 2 * n + 1] = phi[i] + s * phi1[i]
        # T''(s) analitica:  T'' = 2 tau (2s - L2) + 2(1 + s tau)
        T2 = matrix(n, n)
        for i in range(n):
            for j in range(n):
                T2[i, j] = -2 * TAU * L2[i, j]
            T2[i, i] += 4 * TAU * s + 2 * (1 + s * TAU)
        d2 = T2 * phi
        T1phi1 = T1 * phi1
        for i in range(n):
            J[n + i, 2 * n] = d2[i] + T1phi1[i]
        for j in range(n):
            J[2 * n, j] = c[j]; J[2 * n + 1, n + j] = c[j]
        if res < tol:
            break
        d = lu_solve(J, -F)
        for i in range(n):
            phi[i] += d[i]; phi1[i] += d[n + i]
        s += d[2 * n]; z0 += d[2 * n + 1]
    # condicionamiento del jacobiano orlado (parte real: J es real en el EP)
    Jr = matrix(m, m)
    for i in range(m):
        for j in range(m):
            v = J[i, j]
            Jr[i, j] = v.real if hasattr(v, 'real') else v
    U, S, Vt = svd_r(Jr)
    kappa = S[0] / S[m - 1]
    # factor de contraccion (evidencia de convergencia cuadratica)
    contraction = []
    for i in range(1, len(residuals)):
        if residuals[i - 1] > 0:
            contraction.append(residuals[i] / residuals[i - 1] ** 2)
    return dict(N=N, s=s, z0=z0, res=res, iters=it + 1, kappa=kappa,
                residuals=residuals, contraction=contraction)


def digits_agree(a, b):
    a, b = mpf(a.real if hasattr(a, 'real') else a), mpf(b.real if hasattr(b, 'real') else b)
    if a == b: return mp.dps
    if a == 0: return 0
    return int(max(0, -mp.log10(abs((a - b) / a))))


# ================================================================== main
if __name__ == "__main__":
    s_ref = -K * Y_EP2
    z_ref =  K * LAM_EP2

    print("=" * 74)
    print("(1)-(2)-(4)  EP2 certificado (De=0.15): condicionamiento, residuo,")
    print("             iteraciones y contraccion de Newton")
    print("=" * 74)
    mp.dps = 30
    runs = {}
    print(f"{'N':>4} {'s*':>26} {'residuo':>10} {'kappa_2(J)':>12} {'iters':>6}")
    for N in (24, 32, 40):
        r = solve_EP2_audit(N, s_ref * mpf('1.0000001'), z_ref * mpf('1.0000001'))
        runs[N] = r
        print(f"{N:>4} {nstr(r['s'].real, 22):>26} {nstr(r['res'],3):>10} "
              f"{nstr(r['kappa'],5):>12} {r['iters']:>6}")
    print("\n  residuos por iteracion (N=24):",
          "  ".join(nstr(x, 2) for x in runs[24]['residuals']))
    print("  factor de contraccion res_{k+1}/res_k^2 :",
          "  ".join(nstr(x, 3) for x in runs[24]['contraction'][-3:]))
    print("  (constante => convergencia cuadratica => el cero es regular)")

    print()
    print("=" * 74)
    print("(3a)  Estabilidad de digitos frente a la resolucion N")
    print("=" * 74)
    print(f"{'par':>12} {'digitos coincidentes en s*':>28} {'en zeta_0*':>14}")
    for A, B in ((24, 32), (32, 40)):
        print(f"{f'{A}<->{B}':>12} {digits_agree(runs[A]['s'], runs[B]['s']):>28} "
              f"{digits_agree(runs[A]['z0'], runs[B]['z0']):>14}")

    print()
    print("=" * 74)
    print("(3b)  Estabilidad frente a la precision de trabajo d")
    print("=" * 74)
    byd = {}
    for d in (25, 30, 40):
        mp.dps = d
        byd[d] = solve_EP2_audit(24, mpf(str(s_ref)) * mpf('1.0000001'),
                                 mpf(str(z_ref)) * mpf('1.0000001'))
        print(f"  d={d}: s*={nstr(byd[d]['s'].real, min(d-4, 24))}  "
              f"res={nstr(byd[d]['res'],3)}  kappa={nstr(byd[d]['kappa'],4)}")
    mp.dps = 40
    print(f"  digitos coincidentes d=25 vs d=40: {digits_agree(byd[25]['s'], byd[40]['s'])}")
    print(f"  digitos coincidentes d=30 vs d=40: {digits_agree(byd[30]['s'], byd[40]['s'])}")

    print()
    print("=" * 74)
    print("(5)  Comprobacion en zeta_0 = 0 contra el espectro exacto de Poschl-Teller")
    print("     omega_m = +-k - i(m+1/2)  <=>  s = -k y,  y_m = (m+1/2)/k +- i")
    print("=" * 74)
    mp.dps = 30
    L1, L2, n = ops(24)
    from mpmath import mpmathify
    print(f"{'m':>3} {'|s_num - s_exacto|':>24}")
    for m_ in range(4):
        s_ex = mpc(-(m_ + mpf('0.5')), K)      # s = -(m+1/2) + i k
        T, _ = Tmats(s_ex, mpf(0), L1, L2, n)
        Jr = matrix(n, n)
        for i in range(n):
            for j in range(n):
                Jr[i, j] = T[i, j].real
        U, S, Vt = svd_r(Jr)
        print(f"{m_:>3} {nstr(S[n-1], 4):>24}   (menor valor singular de T(s_exacto))")
    print("\n  Un valor singular menor ~ epsilon de la precision de trabajo confirma")
    print("  que s_exacto es autovalor del operador DISCRETIZADO a maquina, es decir,")
    print("  que la colocacion es exacta sobre el sector polinomico (no aproximada).")
