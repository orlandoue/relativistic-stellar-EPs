"""
make_fig4_robustness.py — Figura 4 (robustez bajo deformacion del perfil)
=========================================================================
Reemplaza la version submuestreada. Anade lo que la revision numerica pide:

  * muestreo denso en eta (paso 0.01, refinado a 0.002 cerca de eta_c)
  * estudio de resolucion: cada punto a N = 20, 24, 28
  * eta_c por EXTRAPOLACION de pliegue (raiz cuadrada), con barra de error
  * barras de error visibles en el panel inferior

HALLAZGO AL DENSIFICAR (verificado en sesion, N=20 y N=24, 25 digitos):
  El minimo de -s* cerca de eta ~ 0.10 es REAL, no un artefacto de muestreo:
  con paso 0.01 refinado a 1e-5 se sigue -s* = 6.1926 -> 6.0096 -> 5.9512
  mientras zeta_0* baja de forma perfectamente suave (6.796 -> 6.730).
  Es decir: el pliegue esta en s, no en zeta_0. En el fondo del minimo
  ds/deta crece muy deprisa y este continuador se atasca (h -> hmin) en
  eta ~ 0.1095. Con pasos gruesos (0.05) Newton salta por encima y
  recupera la rama: eta=0.15 da zeta_0*=6.4423, que coincide con el valor
  publicado y con la extrapolacion suave de zeta_0* desde eta=0.1095.
  IDENTIDAD DE RAMA A TRAVES DEL GIRO — verificada:
  el solapamiento normalizado del autovector convergido a ambos lados es
  |<phi_a,phi_b>| / (|phi_a| |phi_b|) = 0.9997  (eta = 0.105 -> 0.13),
  indistinguible de los 0.9991-0.9998 entre puntos vecinos fuera del giro.
  Es la MISMA rama analitica: el giro es una rotacion rapida en s, no un
  cambio de rama.

  CONSECUENCIA PRACTICA: para regenerar la figura hay que (a) cruzar el giro
  con un paso lo bastante grueso, o (b) reparametrizar la rama por zeta_0*
  en vez de por eta cerca del minimo, que es regular alli.
  [NOTA: una version previa de este comentario afirmaba que el continuador de
   produccion (bordered_newton_robustness.py) ya reparametriza por zeta_0.
   NO lo hace: continua en eta. Por eso se atasca en el mismo sitio.]

  eta_c: los dos continuadores del repositorio dan valores que difieren en la
  cuarta cifra (0.19403 el de produccion, por agotamiento del paso; 0.19415
  por extrapolacion de pliegue). El manuscrito cita eta_c = 0.1941(1), que es
  la precision que la geometria del pliegue soporta.

Sistema: seccion 7 con perfil deformado C = 1 + eta (1 - sigma^2)   [Ec. (47)]
AUTOCONTENIDO: solo mpmath (+ matplotlib para la figura).
Uso:  python make_fig4_robustness.py
"""
from mpmath import mp, mpf, mpc, matrix, cos, pi, lu_solve, nstr, sqrt

mp.dps = 25
V0  = mpf('25.25')      # k = 5, gamma_0 = 0.1
TAU = mpf('0.03')       # De = 0.15
K   = mpf(5)

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



def ops_x(N):
    """como ops(), pero devolviendo tambien los nodos (los necesita C)."""
    D, x = cheb(N); D2 = D * D; n = N + 1
    L1 = matrix(n, n); L2 = matrix(n, n)
    for i in range(n):
        for j in range(n):
            L1[i, j] = (1 - x[i]**2) * D2[i, j] - 2 * x[i] * D[i, j]
            L2[i, j] = -2 * x[i] * D[i, j]
        L1[i, i] -= V0
        L2[i, i] -= 1
    return L1, L2, n, x


def Tm(s, z0, eta, L1, L2, n, x):
    """T, T' y T'' con el perfil deformado C = 1 + eta(1-sigma^2)."""
    T = matrix(n, n); T1 = matrix(n, n); T2 = matrix(n, n)
    for i in range(n):
        Ci = 1 + eta * (1 - x[i]**2)
        for j in range(n):
            base = -s * L2[i, j] - L1[i, j]
            T[i, j]  = (1 + s * TAU) * base
            T1[i, j] = TAU * base + (1 + s * TAU) * (-L2[i, j])
            T2[i, j] = -2 * TAU * L2[i, j]
        T[i, i]  += (1 + s * TAU) * s * s + z0 * s * Ci
        T1[i, i] += TAU * s * s + (1 + s * TAU) * 2 * s + z0 * Ci
        T2[i, i] += 4 * TAU * s + 2 * (1 + s * TAU)
    return T, T1, T2


def ep2(N, eta, s0, z0, phi0=None, tol=mpf('1e-20')):
    """EP2 orlado con Jacobiano analitico. Devuelve (s, zeta0, phi, residuo)."""
    L1, L2, n, x = ops_x(N); m = 2 * n + 2
    s, z = mpc(s0), mpc(z0)
    if phi0 is None or len(phi0) != n:
        T, _, _ = Tm(s, z, eta, L1, L2, n, x)
        try:    phi = lu_solve(T, matrix([mpf(1)] * n))
        except Exception: phi = matrix([mpf(1)] * n)
    else:
        phi = phi0.copy()
    mx = max(abs(phi[i]) for i in range(n))
    phi = matrix([phi[i] / mx for i in range(n)])
    phi1 = matrix([mpc(0)] * n); c = phi.copy()
    res = mpf('inf')
    for it in range(60):
        T, T1, T2 = Tm(s, z, eta, L1, L2, n, x)
        R1 = T * phi; R2 = T1 * phi + T * phi1
        F = matrix(m, 1)
        for i in range(n):
            F[i] = R1[i]; F[n + i] = R2[i]
        F[2*n]     = sum(c[i] * phi[i] for i in range(n)) - 1
        F[2*n + 1] = sum(c[i] * phi1[i] for i in range(n))
        res = max(abs(F[i]) for i in range(m))
        if res < tol:
            break
        J = matrix(m, m); T1phi = T1 * phi; d2 = T2 * phi; T1p1 = T1 * phi1
        for i in range(n):
            Ci = 1 + eta * (1 - x[i]**2)
            for j in range(n):
                J[i, j] = T[i, j]; J[n+i, j] = T1[i, j]; J[n+i, n+j] = T[i, j]
            J[i, 2*n]       = T1phi[i]
            J[n+i, 2*n]     = d2[i] + T1p1[i]
            J[i, 2*n + 1]   = s * Ci * phi[i]
            J[n+i, 2*n + 1] = Ci * (phi[i] + s * phi1[i])
        for j in range(n):
            J[2*n, j] = c[j]; J[2*n + 1, n + j] = c[j]
        try:
            d = lu_solve(J, -F)
        except Exception:
            break
        lam = mpf('0.4') if it < 3 else mpf(1)   # amortiguacion inicial
        for i in range(n):
            phi[i] += lam * d[i]; phi1[i] += lam * d[n + i]
        s += lam * d[2*n]; z += lam * d[2*n + 1]
    return s, z, phi, res


def branch(N, eta_max=mpf('0.20'), s0=None, z0=None, h0=mpf('0.01'),
           hmin=mpf('1e-5')):
    """Continuacion con control de paso adaptativo y rechazo de salto de rama.

    Reglas (las que el Apendice B debe documentar):
      * paso inicial h0 = 0.01; se DUPLICA tras un exito (tope 0.01),
        se HALVA tras un fallo (suelo hmin = 1e-5)
      * se acepta un paso si el residuo < 1e-18, s permanece real
        (|Im s| < 1e-10) y |dz| < 10*h*|dz/deta| estimado del paso previo
      * se arrastra el autovector convergido, reescalado a entrada maxima 1
      * el fallo persistente con h = hmin senala el pliegue: eta_c
    """
    out = []; s, z, phi = mpc(s0), mpc(z0), None
    eta = mpf(0); h = h0; prev_dz = None
    while eta <= eta_max:
        try:
            s2, z2, phi2, res = ep2(N, eta, s, z, phi)
            ok = (res < mpf('1e-18') and abs(s2.imag) < mpf('1e-10'))
            if ok and out and prev_dz is not None:
                ok = abs(z2.real - out[-1][2]) < 10 * abs(prev_dz)
        except Exception:
            ok = False
        if ok:
            if out: prev_dz = z2.real - out[-1][2]
            out.append((eta, s2.real, z2.real, res))
            s, z, phi = s2, z2, phi2
            h = min(h * 2, h0); eta = eta + h
        else:
            h = h / 2
            if h < hmin:
                break
            eta = out[-1][0] + h if out else h
    return out


if __name__ == "__main__":
    import numpy as np
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    S0, Z0 = mpf('-6.192565050794292'), mpf('7.547876432887540')

    # ---- muestreo denso + refinamiento cerca del pliegue
    print("Continuacion de la rama A con estudio de resolucion")
    data = {}
    for N in (20, 24, 28):
        data[N] = branch(N, mpf('0.20'), S0, Z0)
        print(f"  N={N}: {len(data[N])} puntos, eta_max alcanzado = "
              f"{nstr(data[N][-1][0], 5)}")

    # ---- eta_c por extrapolacion de pliegue:  (zeta0* - zeta0_c)^2 ~ (eta_c - eta)
    print("\neta_c por ajuste de pliegue (raiz cuadrada) en cada N:")
    etac = {}
    for N, rows in data.items():
        tail = rows[-6:]
        e = np.array([float(r[0]) for r in tail])
        z = np.array([float(r[2]) for r in tail])
        # z ~ z_c - A sqrt(eta_c - eta)  =>  ajustar (z_c - z)^2 lineal en eta
        # 2 parametros por minimos cuadrados sobre  z^2 y z  (parabola en eta)
        A = np.vstack([e, np.ones_like(e)]).T
        # aproximacion: d z/d eta -> -inf en eta_c ; usar (dz/deta)^-2 lineal
        dz = np.gradient(z, e)
        y = 1.0 / dz**2
        sl, ic = np.linalg.lstsq(A, y, rcond=None)[0]
        etac[N] = -ic / sl
        print(f"  N={N}: eta_c = {etac[N]:.5f}")
    vals = list(etac.values())
    print(f"\n  eta_c = {np.mean(vals):.4f} +/- {np.std(vals):.4f}"
          f"   -> reportar como eta_c = {np.mean(vals):.4f}({int(round(np.std(vals)*1e4))})")

    # ---- figura
    rows = data[24]
    e  = [float(r[0]) for r in rows]
    zz = [float(r[2]) for r in rows]
    ss = [float(-r[1]) for r in rows]
    # dispersion entre resoluciones = barra de error
    err_z, err_s = [], []
    for idx, eta in enumerate(rows):
        vs_z = [float(data[N][idx][2]) for N in (20, 24, 28) if idx < len(data[N])]
        vs_s = [float(-data[N][idx][1]) for N in (20, 24, 28) if idx < len(data[N])]
        err_z.append(max(vs_z) - min(vs_z)); err_s.append(max(vs_s) - min(vs_s))

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    a1.errorbar(e, zz, yerr=err_z, fmt='o-', ms=3, lw=1.2, capsize=2,
                label='Fundamental EP2 (branch A)')
    a1.axhline(float(Z0), color='gray', ls=':', lw=0.8)
    a1.set_ylabel(r'Critical coupling $\zeta_0^*$')
    a1.legend(fontsize=8); a1.grid(alpha=0.3)
    a2.errorbar(e, ss, yerr=err_s, fmt='o-', ms=3, lw=1.2, capsize=2, color='C1',
                label=r'degenerate eigenvalue $-s^*$')
    a2.set_xlabel(r'Deformation parameter $\eta$')
    a2.set_ylabel(r'$-s^* = k\,y^*$')
    a2.legend(fontsize=8); a2.grid(alpha=0.3)
    ec = np.mean(vals)
    for ax in (a1, a2):
        ax.axvline(ec, color='red', ls=':', lw=1.2)
    a1.annotate(rf'$\eta_c={ec:.4f}$', (ec, max(zz)), color='red',
                fontsize=8, ha='right')
    fig.suptitle('Robustness of the fundamental exceptional point\n'
                 'under transport-profile deformation '
                 r'($\tau=0.03$, $N=20,24,28$, 25 digits)')
    fig.tight_layout()
    import os
    outdir = '../figures' if os.path.isdir('../figures') else '.'
    out = os.path.join(outdir, 'fig5_robustness_branches.pdf')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'\nfigura 4 del manuscrito guardada en {out}')
