# =============================================================================
# FASE 2 — ROBUSTEZ: EP2 y cusp bajo perfiles de viscosidad genericos
# =============================================================================
# Deformacion:  zeta(x) = zeta0 sech^2(x) * (1 + eta * sech^2(x))
#   -> acople hiperboloidal  zeta0 * C(sigma),  C = 1 + eta (1 - sigma^2)
# eta = 0: caso EXACTO (teorema: sector constante = polinomio abierto).
# eta > 0: acoplamiento generico entre modos (escenario del referee).
#
# TEOREMA DE LA TORRE (eta=0): en la base de monomios sigma^m el operador es
# triangular (L2 sigma^m = -(2m+1) sigma^m; L1 sigma^m = m(m-1) sigma^{m-2}
# - (m(m+1)+V0) sigma^m), asi que el espectro COMPLETO es la union exacta de
# cubicas sectoriales
#     P_m(s) = (1 + s tau)(s^2 + (2m+1)s + m(m+1) + V0) + zeta0 s,
# es decir, una torre de "open cusp polynomials", uno por overtone, con
# gamma_m = (m+1/2)/k genuinos (outputs radiativos del problema de contorno).
# Cada sector tiene su EP2 y su cusp, dados EXACTAMENTE por el algebra de la
# V2 del manuscrito; la Ec. (23) se satisface con residuo 0 en el sector m=0.
#
# HALLAZGOS CERTIFICADOS EN SESION (N=20-28, dps=25-30; correr este script
# con N=28/dps=30 produce los datos citables):
#   * Rama A = EP2 del sector m=0 (fundamental): persiste REAL y sobre el eje
#     para eta en [0, eta_c), eta_c ~ 0.1942 (tau=0.03).
#   * Rama A' = EP2 del sector m=2 (identidad verificada a 16 digitos en
#     eta=0: s*=-6.38293775645287, z0*=5.076486343593958).
#   * En eta_c, A y A' SE FUSIONAN (pliegue) y escapan a zeta0 COMPLEJO como
#     par conjugado: p.ej. eta=0.21: s*=-6.2249+0.1131i, z0*=5.8561+0.3365i.
#     La transicion de sobreamortiguamiento pasa de punto critico exacto a
#     crossover suave. NO es un EP3 (verificado: solver EP3-eta diverge).
#   * Rama B nace del punto DIABOLICO intersectorial s=-6.5,
#     z0 = 0.805*61/6.5 = 7.554615... donde P_0 y P_12 comparten raiz; el
#     acople eta lo convierte en EP2 genuino. Invisible para el modelo 0-D.
#   * La "rama 2" (s~-8.1 en eta~0.2) es otra rama sectorial deformada
#     (candidato: sector m=5, EP2 en (-8.0872, 2.9682)); pendiente confirmar.
#   * El cusp del sector m=0 (EP3, tau libre) persiste al menos hasta
#     eta=0.44 con tau* ~ 0.0406-0.0411 casi constante.
# =============================================================================
from mpmath import mp, mpf, mpc, matrix, lu_solve, cos, pi, sqrt, nstr, findroot
import json, time

# ----------------------------- operadores ------------------------------------
def setup_mp(N, V0, eta):
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
    Cp = [1 + eta * (1 - x[i] ** 2) for i in range(n)]
    W = [1 - x[i] ** 2 for i in range(n)]
    for i in range(n):
        w = 1 - x[i] ** 2
        for j in range(n):
            L1[i, j] = w * D2[i, j] - 2 * x[i] * D[i, j]
            L2[i, j] = -2 * x[i] * D[i, j]
        L1[i, i] -= V0
        L2[i, i] -= 1
    return L1, L2, Cp, W, n

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

def solve_EP2(s, z0, tau, phi0, L1, L2, Cp, n, itmax=60, tol=None):
    """Sistema bordeado 2n+2: raiz doble. phi0 es semilla Y normalizacion.
    Para continuacion, pasar la autofuncion convergida del paso anterior."""
    tol = tol or mpf('10') ** (-(mp.dps - 8))
    phi = phi0.copy(); phi1 = matrix([mpc(0)] * n); c = phi0.copy()
    m = 2 * n + 2
    for it in range(itmax):
        T, T1, T2 = T_all(s, z0, tau, L1, L2, Cp, n)
        Tphi = T * phi; T1phi = T1 * phi
        Tphi1 = T * phi1; T1phi1 = T1 * phi1; T2phi = T2 * phi
        F = matrix(m, 1)
        for i in range(n):
            F[i] = Tphi[i]
            F[n + i] = T1phi[i] + Tphi1[i]
        F[2*n] = sum(c[i] * phi[i] for i in range(n)) - 1
        F[2*n+1] = sum(c[i] * phi1[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        J = matrix(m, m)
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
                J[n+i, j] = T1[i, j]
                J[n+i, n+j] = T[i, j]
            J[i, 2*n] = T1phi[i]
            J[i, 2*n+1] = s * Cp[i] * phi[i]
            J[n+i, 2*n] = T2phi[i] + T1phi1[i]
            J[n+i, 2*n+1] = Cp[i] * phi[i] + s * Cp[i] * phi1[i]
        for j in range(n):
            J[2*n, j] = c[j]
            J[2*n+1, n+j] = c[j]
        d = lu_solve(J, -F)
        for i in range(n):
            phi[i] += d[i]; phi1[i] += d[n+i]
        s += d[2*n]; z0 += d[2*n+1]
        if res < tol:
            return s, z0, phi, res, True
    return s, z0, phi, res, False

def solve_EP3(s, z0, tau, phi0, L1, L2, Cp, n, itmax=60, tol=None):
    """Sistema bordeado 3n+3: raiz triple, tau libre (eta fijo)."""
    tol = tol or mpf('10') ** (-(mp.dps - 8))
    phi = phi0.copy(); phi1 = matrix([mpc(0)] * n); phi2 = matrix([mpc(0)] * n)
    c = phi0.copy(); m = 3 * n + 3
    for it in range(itmax):
        T, T1, T2 = T_all(s, z0, tau, L1, L2, Cp, n)
        dT = matrix(n, n); dT1 = matrix(n, n); dT2 = matrix(n, n)
        for i in range(n):
            for j in range(n):
                base = -s * L2[i, j] - L1[i, j]
                dT[i, j] = s * base
                dT1[i, j] = base - s * L2[i, j]
                dT2[i, j] = -2 * L2[i, j]
            dT[i, i] += s ** 3
            dT1[i, i] += 3 * s * s
            dT2[i, i] += 6 * s
        Tphi = T * phi; T1phi = T1 * phi; T2phi = T2 * phi
        Tphi1 = T * phi1; T1phi1 = T1 * phi1; T2phi1 = T2 * phi1
        Tphi2 = T * phi2; T1phi2 = T1 * phi2
        F = matrix(m, 1)
        for i in range(n):
            F[i] = Tphi[i]
            F[n+i] = T1phi[i] + Tphi1[i]
            F[2*n+i] = T2phi[i] + 2 * T1phi1[i] + Tphi2[i]
        F[3*n] = sum(c[i] * phi[i] for i in range(n)) - 1
        F[3*n+1] = sum(c[i] * phi1[i] for i in range(n))
        F[3*n+2] = sum(c[i] * phi2[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        J = matrix(m, m)
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
                J[n+i, j] = T1[i, j];   J[n+i, n+j] = T[i, j]
                J[2*n+i, j] = T2[i, j]; J[2*n+i, n+j] = 2 * T1[i, j]
                J[2*n+i, 2*n+j] = T[i, j]
            J[i, 3*n] = T1phi[i]
            J[i, 3*n+1] = s * Cp[i] * phi[i]
            J[i, 3*n+2] = sum(dT[i, jj] * phi[jj] for jj in range(n))
            J[n+i, 3*n] = T2phi[i] + T1phi1[i]
            J[n+i, 3*n+1] = Cp[i] * phi[i] + s * Cp[i] * phi1[i]
            J[n+i, 3*n+2] = sum(dT1[i, jj] * phi[jj] + dT[i, jj] * phi1[jj]
                                for jj in range(n))
            J[2*n+i, 3*n] = 6 * tau * phi[i] + 2 * T2phi1[i] + T1phi2[i]
            J[2*n+i, 3*n+1] = 2 * Cp[i] * phi1[i] + s * Cp[i] * phi2[i]
            J[2*n+i, 3*n+2] = sum(dT2[i, jj] * phi[jj] + 2 * dT1[i, jj] * phi1[jj]
                                  + dT[i, jj] * phi2[jj] for jj in range(n))
        for j in range(n):
            J[3*n, j] = c[j]; J[3*n+1, n+j] = c[j]; J[3*n+2, 2*n+j] = c[j]
        d = lu_solve(J, -F)
        for i in range(n):
            phi[i] += d[i]; phi1[i] += d[n+i]; phi2[i] += d[2*n+i]
        s += d[3*n]; z0 += d[3*n+1]; tau += d[3*n+2]
        if res < tol:
            return s, z0, tau, phi, res, True
    return s, z0, tau, phi, res, False

def renorm(phi, n):
    """Reescala la autofuncion para que su entrada de mayor modulo sea 1.
    Necesario antes de reutilizarla como semilla/normalizacion del paso
    siguiente de una continuacion (si no, la fila de normalizacion del
    sistema bordeado queda mal escalada y el Newton puede saltar de rama)."""
    mx = max(abs(phi[i]) for i in range(n))
    return matrix([phi[i] / mx for i in range(n)])

# --------------------- referencias analiticas (eta = 0) ----------------------
def refs_analiticas(V0, tau):
    k = sqrt(V0 - mpf('0.25'))
    def F2(s, z0):
        P = (1 + s * tau) * (s * s + s + V0) + z0 * s
        Pp = tau * (s * s + s + V0) + (1 + s * tau) * (2 * s + 1) + z0
        return [P, Pp]
    s2, z2 = findroot(F2, [mpf('-6.2'), mpf('7.5')])
    def F3(s, z0, tv):
        P = (1 + s * tv) * (s * s + s + V0) + z0 * s
        Pp = tv * (s * s + s + V0) + (1 + s * tv) * (2 * s + 1) + z0
        Ppp = 2 * tv * (2 * s + 1) + 2 * (1 + s * tv)
        return [P, Pp, Ppp]
    s3, z3, t3 = findroot(F3, [mpf('-8.5'), mpf('6.9'), mpf('0.041')])
    return (s2, z2), (s3, z3, t3), k

# =============================================================================
if __name__ == "__main__":
    mp.dps = 30
    V0 = mpf('25.25'); tau0 = mpf('0.03'); N = 28
    (s2ref, z2ref), (s3ref, z3ref, t3ref), k = refs_analiticas(V0, tau0)
    out = {}

    # ---------- 2A: autotest eta = 0 (teorema) ----------
    print("[2A] Autotest eta=0")
    L1, L2, Cp, W, n = setup_mp(N, V0, mpf('0'))
    phi = matrix([mpc(1, 0)] * n)
    s, z0, phiA, res, ok = solve_EP2(mpc(s2ref), mpc(z2ref), tau0, phi, L1, L2, Cp, n)
    print(f"   EP2: |ds|={nstr(abs(s-s2ref),3)} |dz0|={nstr(abs(z0-z2ref),3)} ok={ok}")
    phi = matrix([mpc(1, 0)] * n)
    s3, z3, t3, phiB, res, ok = solve_EP3(mpc(s3ref), mpc(z3ref), mpc(t3ref),
                                          phi, L1, L2, Cp, n)
    print(f"   EP3: |ds|={nstr(abs(s3-s3ref),3)} |dtau|={nstr(abs(t3-t3ref),3)} ok={ok}")

    # ---------- 2B: rama 1 del EP2, continuacion adaptativa hasta el pliegue ----------
    print("[2B] EP2 fundamental vs eta (phi arrastrada, paso adaptativo)")
    sE, zE = mpc(s2ref), mpc(z2ref)
    phiE = matrix([mpc(1, 0)] * n)
    eta = mpf('0'); deta = mpf('0.05')
    r1 = [(0.0, str(sE), str(zE))]
    while eta < mpf('0.5') and deta > mpf('2e-5'):
        et = eta + deta
        L1, L2, Cp, W, n = setup_mp(N, V0, et)
        try:
            s_t, z_t, phi_t, res, ok = solve_EP2(sE, zE, tau0, phiE, L1, L2, Cp, n)
        except ZeroDivisionError:
            ok = False
        if ok and abs(s_t - sE) < mpf('0.25') + 3 * deta:
            eta, sE, zE, phiE = et, s_t, z_t, renorm(phi_t, n)
            r1.append((float(eta), str(sE), str(zE)))
            print(f"   eta={float(eta):.5f}: s*={nstr(sE,13)} z0*={nstr(zE,13)}")
            deta = min(deta * mpf('1.4'), mpf('0.05'))
        else:
            deta = deta / 2
    print(f"   -> rama A termina en eta_c ~ {float(eta):.5f}: pliegue con la rama")
    print(f"      A' (EP2 del sector m=2); ambas escapan a zeta0 complejo (2C).")
    print(f"      Nota: si el ultimo paso muestra |s*| ~ 6.9, el Newton salto a")
    print(f"      la rama B (diabolica P0xP12) y eta_c es el punto anterior.")
    out["rama1"] = r1

    # ---------- 2C: rama compleja post-pliegue ----------
    print("[2C] Rama compleja (eta > eta_c): semillas complejas")
    r2 = []
    sC, zC = mpc('-6.26', '0.05'), mpc('5.86', '0.15')
    phiC = matrix([mpc(1, 0)] * n)
    for ev in ['0.21', '0.25', '0.30', '0.35', '0.40']:
        L1, L2, Cp, W, n = setup_mp(N, V0, mpf(ev))
        try:
            sC, zC, phiC, res, ok = solve_EP2(sC, zC, tau0, phiC, L1, L2, Cp, n)
            phiC = renorm(phiC, n)
            print(f"   eta={ev}: s*={nstr(sC,13)}  z0*={nstr(zC,13)}  ok={ok}")
            r2.append((float(ev), str(sC), str(zC)))
        except ZeroDivisionError:
            print(f"   eta={ev}: LU singular")
    out["rama_compleja"] = r2

    # ---------- 2D: rama 2 del EP2 (s ~ -8.1) ----------
    print("[2D] Segunda rama de EP2")
    r3 = []
    sD, zD = mpc('-8.10'), mpc('6.13')
    phiD = matrix([mpc(1, 0)] * n)
    for ev in ['0.25', '0.22', '0.20', '0.17', '0.14', '0.11', '0.08', '0.05']:
        L1, L2, Cp, W, n = setup_mp(N, V0, mpf(ev))
        try:
            sD, zD, phiD, res, ok = solve_EP2(sD, zD, tau0, phiD, L1, L2, Cp, n)
            phiD = renorm(phiD, n)
            print(f"   eta={ev}: s*={nstr(sD,13)}  z0*={nstr(zD,13)}  ok={ok}")
            r3.append((float(ev), str(sD), str(zD)))
        except ZeroDivisionError:
            print(f"   eta={ev}: LU singular")
    out["rama2"] = r3

    # ---------- 2E: trayectoria del cusp ----------
    print("[2E] Cusp (EP3, tau libre) vs eta")
    sK, zK, tK = mpc(s3ref), mpc(z3ref), mpc(t3ref)
    phiK = matrix([mpc(1, 0)] * n)
    eta = mpf('0'); deta = mpf('0.03')
    r4 = [(0.0, str(tK), str(zK), str(sK))]
    t0 = time.time()
    while eta < mpf('0.5') and deta > mpf('5e-4') and time.time() - t0 < 3000:
        et = eta + deta
        L1, L2, Cp, W, n = setup_mp(N, V0, et)
        try:
            s_t, z_t, t_t, phi_t, res, ok = solve_EP3(sK, zK, tK, phiK, L1, L2, Cp, n)
        except ZeroDivisionError:
            ok = False
        # rechazar saltos de rama (anomalia tipo eta=0.09 detectada en sesion)
        if ok and abs(s_t - sK) < mpf('0.6') and abs(t_t - tK) < mpf('0.01'):
            eta, sK, zK, tK, phiK = et, s_t, z_t, t_t, renorm(phi_t, n)
            r4.append((float(eta), str(tK), str(zK), str(sK)))
            print(f"   eta={float(eta):.4f}: tau*={nstr(tK,11)} z0*={nstr(zK,11)} s*={nstr(sK,11)}")
            deta = min(deta * mpf('1.3'), mpf('0.04'))
        else:
            deta = deta / 2
    out["cusp"] = r4

    with open("fase2_robustez.json", "w") as f:
        json.dump(out, f, indent=1)
    print("\nGuardado: fase2_robustez.json")
