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
      (sustituye al "4e-13", que era una comprobacion en doble precision).
      NOTA: se calcula sigma_min de la matriz COMPLEJA T(s_exacto) via su
      incrustacion real 2n x 2n. Una version previa tomaba Re T(s_exacto), lo
      que no prueba la singularidad de T; el resultado (~1e-29 frente a
      ||T|| ~ 2e3) no cambia, pero el argumento ahora si lo sostiene.

Sistema: el de la seccion 7,  T(s) = (1+s tau)(s^2 - s L2 - L1) + zeta_0 s
con L1 = (1-sigma^2) d^2 - 2 sigma d - V0 ,  L2 = -2 sigma d - 1  (perfil alineado).

AUTOCONTENIDO: solo requiere mpmath.
Uso:  python audit_appendixB.py
"""
from mpmath import mp, mpf, mpc, matrix, cos, pi, lu_solve, nstr, svd_r, sqrt, log10, findroot

# -----------------------------------------------------------------------------
# PRECISION OF THE CONSTANTS  (this used to be a bug; see CHANGELOG)
# -----------------------------------------------------------------------------
# mpf('...') binds a value at whatever mp.dps happens to be WHEN IT IS EVALUATED.
# A previous version of this file built V0, K and TAU at module level, i.e. before
# __main__ raised mp.dps to 30, so they were parsed at the mpmath default of 15
# digits.  V0 = 25.25 and K = 5 are binary-exact and were unaffected; tau = 0.03
# is NOT, so the script silently solved the problem for a tau differing from
# 3/100 in the 17th significant digit.  The bordered solve then saturated at 16
# forward digits for every N, every working precision d and any number of Newton
# iterations -- an artefact that was mistaken for a conditioning ceiling.
#
# All constants are therefore built as EXACT RATIONALS, at the working precision,
# inside consts().  Never hoist them back to module level.
# -----------------------------------------------------------------------------

def consts():
    """(V0, k, tau) as exact rationals at the CURRENT working precision."""
    return mpf(101) / 4, mpf(5), mpf(3) / 100      # 25.25, 5, 0.03


V0 = K = TAU = None        # populated by set_precision(); never bind these with
                           # mpf('...') at module level -- see the note above.


def set_precision(d):
    """Set the working precision AND rebuild the constants at that precision.

    Always use this instead of assigning mp.dps directly: raising mp.dps without
    rebuilding V0/K/TAU leaves them at whatever precision they were created with,
    which is exactly the failure mode documented above.
    """
    global V0, K, TAU
    mp.dps = d
    V0, K, TAU = consts()
    return V0, K, TAU


set_precision(30)          # module import leaves a consistent, usable state


def exact_ep2(V0, TAU):
    """Closed-form sector-0 EP2.

    Eliminating zeta0 between P_0 = 0 and dP_0/ds = 0 leaves the cubic
        2 tau s^3 + (1 + tau) s^2 - V0 = 0
    (for V0 = 101/4, tau = 3/100 this is 6 s^3 + 103 s^2 - 2525 = 0), whose
    relevant root gives the exact degenerate frequency; zeta0 follows from
    P_0 = 0.  This is the reference against which the operator-level bordered
    solve is measured, and it is exact to arbitrary precision.
    """
    s = findroot(lambda z: 2 * TAU * z ** 3 + (1 + TAU) * z ** 2 - V0, mpf('-6.19'))
    z0 = -(1 + s * TAU) * (s ** 2 + s + V0) / s
    return s, z0


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
                smax=S[0], smin=S[m - 1],
                residuals=residuals, contraction=contraction)


def fwd_digits(a, b):
    """Significant digits of agreement between a and b (b the reference)."""
    a = mpf(a.real if hasattr(a, 'real') else a)
    b = mpf(b.real if hasattr(b, 'real') else b)
    if a == b:
        return float(mp.dps)
    return float(-mp.log10(abs((a - b) / b)))


def digits_agree(a, b):
    a, b = mpf(a.real if hasattr(a, 'real') else a), mpf(b.real if hasattr(b, 'real') else b)
    if a == b: return mp.dps
    if a == 0: return 0
    return int(max(0, -mp.log10(abs((a - b) / a))))


# ================================================================== main
if __name__ == "__main__":
    set_precision(60)
    S_EXACT, Z_EXACT = exact_ep2(V0, TAU)
    print("Referencia exacta (raiz de 2*tau*s^3+(1+tau)*s^2-V0 = 0,")
    print("es decir 6 s^3 + 103 s^2 - 2525 = 0 para V0=101/4, tau=3/100):")
    print("   s*     =", nstr(S_EXACT, 40))
    print("   zeta0* =", nstr(Z_EXACT, 40))
    print("   y*     =", nstr(-S_EXACT / K, 22),
          "   Lambda* =", nstr(Z_EXACT / K, 22))
    print("   (Tabla 1 publica 20 digitos:  1.2385130101588585499 / 1.5095752865775080251)")
    print()
    s_ref, z_ref = S_EXACT, Z_EXACT

    print("=" * 74)
    print("(1)-(2)-(4)  EP2 certificado (De=0.15): condicionamiento, residuo,")
    print("             iteraciones y contraccion de Newton")
    print("=" * 74)
    set_precision(30)
    runs = {}
    print(f"{'N':>4} {'residuo':>10} {'kappa_2(J)':>12} {'iters':>6}"
          f" {'digitos exactos s*':>20} {'zeta0*':>10} {'cota a priori':>14}")
    for N in (24, 32, 40):
        r = solve_EP2_audit(N, s_ref * mpf('1.0000001'), z_ref * mpf('1.0000001'))
        runs[N] = r
        set_precision(60)
        ds = fwd_digits(r['s'], S_EXACT); dz = fwd_digits(r['z0'], Z_EXACT)
        bound = -mp.log10(r['kappa'] * r['res'] / r['smax'] / abs(S_EXACT))
        set_precision(30)
        print(f"{N:>4} {nstr(r['res'],3):>10} {nstr(r['kappa'],5):>12} {r['iters']:>6}"
              f" {ds:>20.1f} {dz:>10.1f} {float(bound):>14.1f}")
    print("\n  La precision directa DECRECE con N siguiendo kappa_2(N): el error de")
    print("  discretizacion es identicamente nulo sobre el sector polinomico, de modo")
    print("  que todo el error es aritmetico. La cota a priori kappa_2*||F||/||J|| es")
    print("  ajustada a ~1 digito en las tres resoluciones.")
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
    set_precision(60)
    for A, B in ((24, 32), (32, 40)):
        print(f"{f'{A}<->{B}':>12} {fwd_digits(runs[A]['s'], runs[B]['s']):>28.1f} "
              f"{fwd_digits(runs[A]['z0'], runs[B]['z0']):>14.1f}")
    print("\n  Coincide con la precision directa de la corrida mas gruesa, como debe")
    print("  ocurrir cuando el error es aritmetico y monotono en kappa_2.")
    set_precision(30)

    print()
    print("=" * 74)
    print("(3b)  Estabilidad frente a la precision de trabajo d")
    print("=" * 74)
    byd = {}
    for d in (25, 30, 40, 50):
        set_precision(d)
        sx, zx = exact_ep2(V0, TAU)
        byd[d] = solve_EP2_audit(24, sx * mpf('1.0000001'), zx * mpf('1.0000001'))
        set_precision(60)
        print(f"  d={d:>2}: res={nstr(byd[d]['res'],3):>10}  kappa={nstr(byd[d]['kappa'],4)}"
              f"   digitos exactos: s*={fwd_digits(byd[d]['s'], S_EXACT):.1f}"
              f"   zeta0*={fwd_digits(byd[d]['z0'], Z_EXACT):.1f}")
    print("\n  La precision directa sigue a la precision de trabajo uno a uno: NO hay")
    print("  techo intrinseco. d>=30 se usa en todo el trabajo porque d=25 ya no")
    print("  alcanza los veinte digitos publicados en la Tabla 1.")

    print()
    print("=" * 74)
    print("(5)  Comprobacion en zeta_0 = 0 contra el espectro exacto de Poschl-Teller")
    print("     omega_m = +-k - i(m+1/2)  <=>  s = -k y,  y_m = (m+1/2)/k +- i")
    print("=" * 74)
    set_precision(30)
    L1, L2, n = ops(24)
    # CORRECCION metodologica: s_exacto es COMPLEJO, de modo que T(s_exacto) es
    # una matriz compleja. Tomar su parte real y calcular sigma_min de esa
    # matriz real NO es una prueba de que T(s_exacto) sea singular. Se usa la
    # incrustacion real 2n x 2n de la matriz compleja,
    #     A = [[Re T, -Im T], [Im T, Re T]],
    # cuyos valores singulares son los de T duplicados, y se normaliza por
    # ||T||_max para que el numero reportado sea relativo y por tanto
    # comparable con el epsilon de la precision de trabajo.
    print(f"{'m':>3} {'sigma_min(T)':>16} {'||T||_max':>12} {'relativo':>12}")
    for m_ in range(4):
        s_ex = mpc(-(m_ + mpf('0.5')), K)      # s = -(m+1/2) + i k
        T, _ = Tmats(s_ex, mpf(0), L1, L2, n)
        A = matrix(2 * n, 2 * n)
        for i in range(n):
            for j in range(n):
                a = T[i, j].real; b = T[i, j].imag
                A[i, j] = a;      A[i, n + j] = -b
                A[n + i, j] = b;  A[n + i, n + j] = a
        U, S, Vt = svd_r(A)
        nrm = max(abs(T[i, j]) for i in range(n) for j in range(n))
        print(f"{m_:>3} {nstr(S[2*n-1], 4):>16} {nstr(nrm, 4):>12} "
              f"{nstr(S[2*n-1]/nrm, 4):>12}")
    print("\n  Un valor singular menor RELATIVO ~ epsilon de la precision de trabajo confirma")
    print("  que s_exacto es autovalor del operador DISCRETIZADO a maquina, es decir,")
    print("  que la colocacion es exacta sobre el sector polinomico (no aproximada).")
