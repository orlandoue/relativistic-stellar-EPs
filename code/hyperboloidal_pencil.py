# =============================================================================
# FASE 1 — Puntos excepcionales del sistema abierto con relajacion causal
# =============================================================================
# Modelo (problema de contorno con ondas salientes genuinas, tipo Siegert):
#
#     d_t^2 phi - d_x^2 phi + V0 sech^2(x) phi + Pi = 0
#     tau d_t Pi + Pi = zeta0 sech^2(x) d_t phi
#
# En frecuencia el termino viscoso es  -i w zeta(x)/(1 - i w tau) phi:
# la estructura cubica EXACTA del manuscrito, pero como problema QNM.
# gamma radiativo genuino:  QNM exactos (zeta0=0):  w_n = ±k - i(n+1/2),
# k = sqrt(V0 - 1/4),  gamma = (1/2)/k  (controlable con V0).
#
# Formulacion hiperboloidal (t = tau_h + ln cosh x, sigma = tanh x):
# regularidad en sigma = ±1 reemplaza las condiciones de Siegert.
# Pencil cubico reducido (n x n), s = -i w:
#     T(s) = (1 + s tau)(s^2 I - s L2 - L1) + zeta0 s I
#     L1 = (1-s^2)D2 - 2 s D - V0 I ,   L2 = -2 s D - I   (s = sigma aqui)
#
# VERIFICACIONES YA REALIZADAS (2026-07-05, sesion con Claude):
#   * ansatz constante reproduce s^2+s+V0=0 (residuo 1e-11)
#   * fundamental exacto a 3e-13 en float64
#   * degradacion con el tono = no-normalidad fisica (usar mp cerca de EPs)
#   * EP2 certificado en V0=25.25, tau=0.03 (De=0.15, gamma=0.1):
#       zeta0* = 7.5478764328875401254...   (Im ~ 1e-322: real)
#       w*     = -6.1925650507942927 i      (eje imaginario)
#     coincidente en 19 cifras entre N=24,32,40; residuo ~1e-29.
# =============================================================================
import numpy as np
import scipy.linalg as sla
import matplotlib.pyplot as plt
from mpmath import mp, mpf, mpc, matrix, lu_solve, cos, pi, sqrt, nstr
import time, json

# ------------------------------ float64 -------------------------------------
def cheb_np(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x

def build_A_np(N, V0, z0, tau):
    D, x = cheb_np(N); D2 = D @ D; n = N + 1
    I = np.eye(n); S = np.diag(x); W = np.diag(1 - x ** 2)
    L1 = W @ D2 - 2 * S @ D - V0 * I
    L2 = -2 * S @ D - I
    Z = np.zeros((n, n))
    A = np.block([[Z, I, Z],
                  [L1, L2, -I],
                  [Z, (z0 / tau) * I, -(1.0 / tau) * I]])
    return A, x

def spectrum_np(N, V0, z0, tau):
    A, _ = build_A_np(N, V0, z0, tau)
    return 1j * sla.eigvals(A)          # omega = i s

# ------------------------------ mpmath --------------------------------------
def setup_mp(N, V0):
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
        L1[i, i] -= V0
        L2[i, i] -= 1
    return L1, L2, n

def T_all(s, z0, tau, L1, L2, n):
    """T(s), T'(s), T''(s) del pencil cubico."""
    T = matrix(n, n); T1 = matrix(n, n); T2 = matrix(n, n)
    a = 1 + s * tau
    for i in range(n):
        for j in range(n):
            base = -s * L2[i, j] - L1[i, j]
            T[i, j]  = a * base
            T1[i, j] = tau * base - a * L2[i, j]
            T2[i, j] = -2 * tau * L2[i, j]
        T[i, i]  += a * s * s + z0 * s
        T1[i, i] += tau * s * s + 2 * a * s + z0
        T2[i, i] += 4 * tau * s + 2 * a
    return T, T1, T2

def newton_eig_mp(s, phi, z0, tau, L1, L2, n, itmax=30, tol=None):
    """Refina un autopar simple (phi, s) de T(s)phi = 0 a zeta0 fijo."""
    tol = tol or mpf('10') ** (-(mp.dps - 6))
    c = phi.copy()
    for it in range(itmax):
        T, T1, _ = T_all(s, z0, tau, L1, L2, n)
        F = T * phi
        T1phi = T1 * phi
        J = matrix(n + 1, n + 1)
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
            J[i, n] = T1phi[i]
            J[n, i] = c[i]
        rhs = matrix(n + 1, 1)
        for i in range(n):
            rhs[i] = -F[i]
        rhs[n] = 1 - sum(c[i] * phi[i] for i in range(n))
        d = lu_solve(J, rhs)
        for i in range(n):
            phi[i] += d[i]
        s += d[n]
        if abs(d[n]) < tol:
            return s, phi, True
    return s, phi, False

def solve_EP2_mp(s, z0, phi0, tau, L1, L2, n, itmax=30, tol=None):
    """Sistema bordeado del EP2:  T phi = 0,  T' phi + T phi1 = 0,
    c.phi = 1, c.phi1 = 0.  Incognitas COMPLEJAS (phi, phi1, s, z0).
    Si converge con Im(z0) ~ 0 => EP en el dominio fisico real."""
    tol = tol or mpf('10') ** (-(mp.dps - 8))
    phi = phi0.copy(); phi1 = matrix([mpc(0)] * n); c = phi0.copy()
    m = 2 * n + 2
    res = mpf('inf')
    for it in range(itmax):
        T, T1, T2 = T_all(s, z0, tau, L1, L2, n)
        Tphi = T * phi; T1phi = T1 * phi
        Tphi1 = T * phi1; T1phi1 = T1 * phi1; T2phi = T2 * phi
        F = matrix(m, 1)
        for i in range(n):
            F[i] = Tphi[i]
            F[n + i] = T1phi[i] + Tphi1[i]
        F[2 * n] = sum(c[i] * phi[i] for i in range(n)) - 1
        F[2 * n + 1] = sum(c[i] * phi1[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        J = matrix(m, m)
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
                J[n + i, j] = T1[i, j]
                J[n + i, n + j] = T[i, j]
            J[i, 2 * n] = T1phi[i]
            J[i, 2 * n + 1] = s * phi[i]                  # dT/dz0 = s I
            J[n + i, 2 * n] = T2phi[i] + T1phi1[i]
            J[n + i, 2 * n + 1] = phi[i] + s * phi1[i]    # dT'/dz0 = I
        for j in range(n):
            J[2 * n, j] = c[j]
            J[2 * n + 1, n + j] = c[j]
        d = lu_solve(J, -F)
        for i in range(n):
            phi[i] += d[i]; phi1[i] += d[n + i]
        s += d[2 * n]; z0 += d[2 * n + 1]
        if res < tol:
            return s, z0, phi, phi1, res, True
    return s, z0, phi, phi1, res, False

def solve_EP3_mp(s, z0, tau, phi0, V0, L1, L2, n, itmax=40, tol=None):
    """Sistema bordeado del CUSP (EP3): agrega T''phi + 2T'phi1 + Tphi2 = 0
    y libera tau como incognita ademas de zeta0. Incognitas:
    (phi, phi1, phi2, s, z0, tau) -> 3n+3 ecuaciones complejas.
    Si converge con Im(z0), Im(tau) ~ 0 => cusp A3 en el dominio fisico."""
    tol = tol or mpf('10') ** (-(mp.dps - 8))
    phi = phi0.copy(); phi1 = matrix([mpc(0)] * n); phi2 = matrix([mpc(0)] * n)
    c = phi0.copy(); m = 3 * n + 3
    for it in range(itmax):
        T, T1, T2 = T_all(s, z0, tau, L1, L2, n)
        # derivadas parciales respecto de tau (T depende de tau explicitamente):
        # T   = (1+s tau)(s^2 I - s L2 - L1) + z0 s I  -> dT/dtau  = s(s^2 I - sL2 - L1)
        # T1  -> dT1/dtau = (s^2 I - sL2 - L1) + s(2sI - L2)
        # T2  -> dT2/dtau = 2(2s I - L2) + 2 s I = (6s I - 2 L2)... ver abajo
        dTdt = matrix(n, n); dT1dt = matrix(n, n); dT2dt = matrix(n, n)
        for i in range(n):
            for j in range(n):
                base = -s * L2[i, j] - L1[i, j]
                dTdt[i, j] = s * base
                dT1dt[i, j] = base - s * L2[i, j]
                dT2dt[i, j] = -2 * L2[i, j]
            dTdt[i, i] += s * s * s
            dT1dt[i, i] += 3 * s * s          # d/dtau[tau s^2 + 2(1+s tau)s] = 3s^2
            dT2dt[i, i] += 6 * s              # d/dtau[4 tau s + 2(1+s tau)] = 6s
        # T''' (tercera derivada en s):  T''' = 6 tau I
        F = matrix(m, 1)
        Tphi = T * phi; T1phi = T1 * phi; T2phi = T2 * phi
        Tphi1 = T * phi1; T1phi1 = T1 * phi1; T2phi1 = T2 * phi1
        Tphi2 = T * phi2; T1phi2 = T1 * phi2
        for i in range(n):
            F[i] = Tphi[i]
            F[n + i] = T1phi[i] + Tphi1[i]
            F[2 * n + i] = T2phi[i] + 2 * T1phi1[i] + Tphi2[i]
        F[3 * n] = sum(c[i] * phi[i] for i in range(n)) - 1
        F[3 * n + 1] = sum(c[i] * phi1[i] for i in range(n))
        F[3 * n + 2] = sum(c[i] * phi2[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        J = matrix(m, m)
        for i in range(n):
            for j in range(n):
                J[i, j] = T[i, j]
                J[n + i, j] = T1[i, j];      J[n + i, n + j] = T[i, j]
                J[2*n + i, j] = T2[i, j];    J[2*n + i, n + j] = 2 * T1[i, j]
                J[2*n + i, 2*n + j] = T[i, j]
            # columnas s, z0, tau
            J[i, 3*n] = T1phi[i]
            J[i, 3*n + 1] = s * phi[i]
            J[i, 3*n + 2] = (dTdt * phi)[i] if False else sum(dTdt[i, jj] * phi[jj] for jj in range(n))
            J[n + i, 3*n] = T2phi[i] + T1phi1[i]
            J[n + i, 3*n + 1] = phi[i] + s * phi1[i]
            J[n + i, 3*n + 2] = sum(dT1dt[i, jj] * phi[jj] + dTdt[i, jj] * phi1[jj] for jj in range(n))
            J[2*n + i, 3*n] = 6 * tau * phi[i] + 2 * (T2 * phi1)[i] if False else (6 * tau * phi[i] + 2 * T2phi1[i] + T1phi2[i])
            J[2*n + i, 3*n + 1] = 2 * phi1[i] + s * phi2[i]
            J[2*n + i, 3*n + 2] = sum(dT2dt[i, jj] * phi[jj] + 2 * dT1dt[i, jj] * phi1[jj] + dTdt[i, jj] * phi2[jj] for jj in range(n))
        for j in range(n):
            J[3*n, j] = c[j]; J[3*n + 1, n + j] = c[j]; J[3*n + 2, 2*n + j] = c[j]
        d = lu_solve(J, -F)
        for i in range(n):
            phi[i] += d[i]; phi1[i] += d[n + i]; phi2[i] += d[2 * n + i]
        s += d[3 * n]; z0 += d[3 * n + 1]; tau += d[3 * n + 2]
        if res < tol:
            return s, z0, tau, res, True
    return s, z0, tau, res, False

# =============================================================================
# DRIVER
# =============================================================================
if __name__ == "__main__":
    mp.dps = 30

    # ---------------- FASE A: validacion exacta ----------------
    V0f = 25.25; tauf = 0.03
    k = np.sqrt(V0f - 0.25); gam = 0.5 / k; De = k * tauf
    print(f"[A] Validacion. k={k}, gamma={gam}, De={De}")
    w = spectrum_np(110, V0f, 0.0, tauf)
    for nn in range(3):
        wex = k - 1j * (nn + 0.5)
        print(f"    n={nn}: err = {np.min(np.abs(w - wex)):.2e}")

    # ---------------- FASE B: paisaje rapido float64 ----------------
    # (solo orientativo; cerca de EPs el float64 es ruido -> ver FASE C)
    print("[B] Paisaje float64 (orientativo)...")
    zs = np.linspace(0, 14, 57)
    fund = []
    prev = k - 0.5j
    for z in zs:
        w = spectrum_np(90, V0f, z, tauf)
        prev = w[np.argmin(np.abs(w - prev))]
        fund.append(prev)
    fund = np.array(fund)
    plt.figure(figsize=(6, 4))
    plt.plot(fund.real, fund.imag, ".-")
    plt.xlabel("Re w"); plt.ylabel("Im w"); plt.grid(alpha=.3)
    plt.title("Trayectoria del fundamental (float64, orientativa)")
    plt.savefig("faseB_paisaje.png", dpi=140)

    # ---------------- FASE C: certificacion del EP2 ----------------
    print("[C] Certificacion multiprecision del EP2...")
    resultados = {}
    for N in (24, 32, 40):
        L1, L2, n = setup_mp(N, mpf(str(V0f)))
        s = mpc('-5.417', '-2.544')          # guess desde el paisaje
        z0 = mpc('6.95', '0')
        phi = matrix([mpc(1, 0)] * n)
        t0 = time.time()
        s_ep, z_ep, phi, phi1, res, ok = solve_EP2_mp(s, z0, phi,
                                                      mpf(str(tauf)), L1, L2, n)
        om = mpc(0, 1) * s_ep
        print(f"    N={N}: res={nstr(res,3)} [{time.time()-t0:.0f}s]")
        print(f"       zeta0* = {nstr(z_ep, 22)}")
        print(f"       w*     = {nstr(om, 22)}")
        resultados[N] = (str(z_ep), str(om))
        if N == 32:
            z_star_mp, s_star_mp = z_ep, s_ep   # valores certificados en memoria
    with open("faseC_EP2_certificado.json", "w") as f:
        json.dump(resultados, f, indent=1)

    # ---------------- FASE D: exponente de escalamiento ----------------
    # alrededor del EP certificado: dw ~ eps^(1/2) para EP2 genuino.
    print("[D] Test de escalamiento |eps|^(1/2)...")
    N = 32
    L1, L2, n = setup_mp(N, mpf(str(V0f)))
    z_star = z_star_mp.real          # certificado en Fase C (Im ~ 1e-300)
    s_star = s_star_mp
    # los dos autovalores cercanos para zeta0 = z* (1 - eps):
    datos = []
    for ee in range(2, 12, 2):
        eps = mpf(10) ** (-ee)
        z0 = z_star * (1 - eps)
        # guess: separacion sqrt: s ~ s* ± C sqrt(eps)
        for sgn in (+1, -1):
            sg = s_star + sgn * mpc(0, 1) * sqrt(eps) * 3
            phi = matrix([mpc(1, 0)] * n)
            sg, phi, ok = newton_eig_mp(sg, phi, z0, mpf(str(tauf)), L1, L2, n)
            if sgn == +1:
                s_plus = sg
            else:
                datos.append((float(eps), float(abs(s_plus - sg))))
    print("    eps        |ds|       pendiente local")
    for i, (e, d) in enumerate(datos):
        pend = ""
        if i > 0:
            pend = f"{(np.log(d) - np.log(datos[i-1][1]))/(np.log(e) - np.log(datos[i-1][0])):.4f}"
        print(f"    {e:.1e}  {d:.6e}   {pend}")
    # pendiente ~ 0.5 => EP2 certificado dinamicamente.

    # ---------------- FASE E: mapeo EP*(De, gamma) ----------------
    # el test decisivo de las Ecs. (22)-(23) del manuscrito: variar
    # tau (De) y V0 (gamma) y seguir el EP2 certificado por continuacion.
    print("[E] Mapeo EP*(De, gamma) -- editar la grilla segun convenga")
    N = 28
    L1, L2, n = setup_mp(N, mpf(str(V0f)))
    filas = []
    # dos ramas monotonas desde el EP certificado en tau=0.03:
    rama_arriba = [0.045, 0.06, 0.09]
    rama_abajo  = [0.02, 0.012]
    for rama in (rama_arriba, rama_abajo):
        s_g, z_g = s_star_mp, z_star_mp   # reiniciar en el punto certificado
        for tv in rama:
            phi = matrix([mpc(1, 0)] * n)
            s_ep, z_ep, phi, phi1, res, ok = solve_EP2_mp(s_g, z_g, phi,
                                                          mpf(str(tv)), L1, L2, n)
            Dev = k * tv
            filas.append((tv, Dev, complex(z_ep), complex(mpc(0, 1) * s_ep),
                          float(res), bool(ok)))
            print(f"    tau={tv}: De={Dev:.4f}  zeta0*={complex(z_ep):.10f}  "
                  f"w*={complex(mpc(0,1)*s_ep):.10f}  ok={ok}")
            s_g, z_g = s_ep, z_ep          # continuacion dentro de la rama
    np.save("faseE_mapa_EP.npy", np.array(filas, dtype=object))

    print("\nListo. Siguientes fases (con Claude): comparacion cuantitativa con")
    print("las Ecs. (22)-(23) del manuscrito, busqueda del segundo EP (colision")
    print("con la rama de relajacion) y del cusp EP3 con solve_EP3_mp.")
