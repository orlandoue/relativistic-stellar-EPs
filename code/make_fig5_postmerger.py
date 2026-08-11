"""
Figure 5: post-merger remnant QNM evolution — CAUSAL fiducial point.
=====================================================================
VERSION AUTOCONTENIDA: no requiere cqg_core.py (solo numpy + matplotlib).
En el manuscrito final esta es la FIGURA 5 (el archivo conserva el nombre
fig4_physical_NS_application.pdf por continuidad con los datos archivados).
Las funciones roots_np y find_ep2 estan definidas abajo.

Solo necesita numpy + matplotlib (ya vienen en Colab).

PARAMETROS ACTUALIZADOS tras la correccion de causalidad de la sec. 9.5
  f0     : 3.5 kHz -> 3.0 kHz     (impuesto por c_s <= 0.341 c)
  k_eff  : 1/R     -> pi/R        (modo fundamental con nodo en superficie)
  Cproj  : 5.32e-32 -> 6.12629e-31
  De     : 0.11    -> 0.17        (cerca del cusp De* = 0.19576)
  zeta   : ~3e31   -> ~2.6e30

CAUSALIDAD (Israel-Stewart: v_Pi^2 + c_s^2 <= c^2):
  c_s = 2 f0 R = 0.300 c ; ambos EP2 con suma ~0.84-0.87 < 1  -> CAUSAL
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===================================================================
# NUCLEO ALGEBRAICO (antes en cqg_core.py)
# ===================================================================
# Open Cusp Polynomial, forma expandida:
#   (1 - i De Om)(Om^2 + 2i g Om - (1+g^2)) + i Lam Om = 0
#   = -i De Om^3 + (1+2 De g) Om^2 + i[2g + De(1+g^2) + Lam] Om - (1+g^2)

def poly_coeffs(De, g, Lam):
    """Coeficientes del Open Cusp Polynomial en Omega (orden descendente)."""
    return [-1j*De,
            1.0 + 2.0*De*g,
            1j*(2.0*g + De*(1.0 + g**2) + Lam),
            -(1.0 + g**2)]

def roots_np(De, g, Lam):
    """Las tres raices Omega del polinomio."""
    return np.roots(poly_coeffs(De, g, Lam))

def _De_of_y(y, g):
    """Parametrizacion EP2 a gamma finito: De(y;gamma)."""
    return (y**2 - (1.0 + g**2)) / (2.0 * y**2 * (y - g))

def _Lam_of_y(y, g):
    """Parametrizacion EP2 a gamma finito: Lambda(y;gamma)."""
    f = 3.0*y**2 - 4.0*g*y + 1.0 + g**2
    return 2.0*(y - g) - _De_of_y(y, g)*f

def find_ep2(De, g, y_seed):
    """Localiza un EP2: devuelve (y*, Lambda*) resolviendo De(y;g) = De
       por Newton escalar desde y_seed. Omega* = -i y*."""
    y = float(y_seed)
    for _ in range(200):
        f = _De_of_y(y, g) - De
        h = 1e-8*max(1.0, abs(y))
        df = (_De_of_y(y + h, g) - _De_of_y(y - h, g)) / (2*h)
        if df == 0:
            break
        step = f/df
        y -= step
        if abs(step) < 1e-14:
            break
    return y, _Lam_of_y(y, g)

# ===================================================================
# PARAMETROS FIDUCIALES (causales)
# ===================================================================
De, g, f0 = 0.17, 0.03, 3.0e3
Cproj = 6.12629e-31          # = 1/(rho_bar * L^2 * omega_0), L = R/pi
C_LIGHT = 2.99792458e10      # cm/s
R_CM = 15.0e5                # cm
cs = 2.0 * f0 * R_CM         # velocidad de fase del modo = omega_0/k_eff

# EP2 exactos de este (De, gamma). Semillas: y ~ 1.33 y y ~ 2.50
y1, L1 = find_ep2(De, g, 1.33)
y2, L2 = find_ep2(De, g, 2.50)
z1, z2 = L1/Cproj, L2/Cproj

# cusp de referencia (De* = 0.19576, Lambda* = 1.48704)
LAM_CUSP = 1.4870382
ZETA_CUSP = LAM_CUSP/Cproj

print(f"c_s = {cs/C_LIGHT:.4f} c   (debe ser < 0.341 c)")
print(f"EP2 #1: y* = {y1:.6f}  Lambda* = {L1:.6f}  zeta = {z1:.4e} g/cm/s")
print(f"EP2 #2: y* = {y2:.6f}  Lambda* = {L2:.6f}  zeta = {z2:.4e} g/cm/s")
print(f"cusp  : Lambda* = {LAM_CUSP:.6f}  zeta = {ZETA_CUSP:.4e} g/cm/s")
for lab, L in (("#1", L1), ("#2", L2)):
    vPi2 = (L/De) * cs**2
    tot = vPi2/C_LIGHT**2 + (cs/C_LIGHT)**2
    print(f"  EP2 {lab}: Lambda/De = {L/De:.4f}, v_Pi = {np.sqrt(vPi2)/C_LIGHT:.4f} c,"
          f" v_Pi^2+c_s^2 = {tot:.4f} {'CAUSAL' if tot <= 1 else '*** VIOLA ***'}")

# ===================================================================
# BARRIDO
# ===================================================================
base = np.linspace(2.0e30, 3.4e30, 600)
dense = np.concatenate([np.linspace(z1*0.985, z1*1.015, 500),
                        np.linspace(z2*0.985, z2*1.015, 500)])
zetas = np.sort(np.concatenate([base, dense]))
Lams = Cproj * zetas

F = np.zeros((len(Lams), 3), dtype=complex)
for i, L in enumerate(Lams):
    r = roots_np(De, g, L)
    F[i] = r[np.argsort(r.imag)]
Ref = F.real * f0
gap = np.minimum(np.abs(F[:, 0] - F[:, 1]), np.abs(F[:, 1] - F[:, 2])) * f0

# ===================================================================
# FIGURA
# ===================================================================
SC = 1e30
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

for j, col in enumerate(['C0', 'C1', 'C2']):
    ax1.plot(zetas/SC, Ref[:, j], '.', color=col, ms=2, label=f'Mode {j+1}')

for ax in (ax1, ax2):
    ax.axvline(z1/SC, color='red', ls=':', lw=1.3, alpha=0.7)
    ax.axvline(z2/SC, color='purple', ls=':', lw=1.3, alpha=0.7)
    ax.axvline(ZETA_CUSP/SC, color='k', ls='--', lw=1.0, alpha=0.5)
    ax.set_xlabel(r'Bulk viscosity ($10^{30}$ g cm$^{-1}$ s$^{-1}$)')
    ax.grid(alpha=0.3)

ax1.set_ylabel('Re($f$) (Hz)')
ax1.set_title('Real frequency across the two coalescences')
ax1.legend(fontsize=8)
ax1.twiny().set_xlim(Cproj*base[0], Cproj*base[-1])

ax2.semilogy(zetas/SC, gap + 1, '-', color='darkgreen', lw=1.5,
             label='min eigenvalue gap')
ax2.annotate(rf'EP2 #1  $\Lambda$={L1:.3f}', (z1/SC, 4e3),
             fontsize=8, color='red', ha='right')
ax2.annotate(rf'EP2 #2  $\Lambda$={L2:.3f}', (z2/SC, 4e3),
             fontsize=8, color='purple', ha='left')
ax2.annotate(rf'$A_3$ cusp  $\Lambda_*$={LAM_CUSP:.3f}', (ZETA_CUSP/SC, 1.2e3),
             fontsize=8, color='k', ha='center')
ax2.set_ylabel('Eigenvalue gap (Hz)')
ax2.set_title('Gap collapses at both exact EP2 coalescences')
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3, which='both')
ax2.twiny().set_xlim(Cproj*base[0], Cproj*base[-1])

fig.suptitle(r'Post-merger remnant QNM evolution: causal fiducial point '
             rf'($f_0$=3.0 kHz, $\tau$={De/(2*np.pi*f0)*1e3:.4f} ms, '
             rf'$\gamma_{{GW}}$={g}, $v_\Pi^2+c_s^2<1$)')
fig.tight_layout()

# En Colab: guarda en el directorio actual. Cambia la ruta si hace falta.
import os
outdir = '../figures' if os.path.isdir('../figures') else '.'
outpath = os.path.join(outdir, 'fig4_physical_NS_application.pdf')
fig.savefig(outpath, dpi=150, bbox_inches='tight')
print(f"fig4 saved -> {outpath}")
