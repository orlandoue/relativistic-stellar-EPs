import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# PARAMETRIZACIÓN FÍSICA ASOCIADA A QCD Y REMANENTE
# ==========================================
f0 = 3500.0                       # Frecuencia central compresiva g-mode (Hz)
omega0 = 2 * np.pi * f0           # Frecuencia angular fundamental (rad/s)
tau = 0.005e-3                    # Tiempo de relajación causal causal (s)
De = omega0 * tau                 # Deborah number del remanente (~0.11)
gamma_GW = 0.03                   # Acoplamiento abierto de pérdidas por radiación

# Dominio físico de barrido para la viscosidad bulk colosal de QCD
zeta_vals = np.linspace(1.5, 4.7, 1500)  # Unidades de 10^31 g/cm/s
# Mapeo lineal exacto al parámetro de acoplamiento adimensional Lambda
lambda_vals = (8.0 / 15.0) * zeta_vals - (4.0 / 15.0)

# Arrays de almacenamiento espectral
re_roots = np.zeros((len(zeta_vals), 3))
im_roots = np.zeros((len(zeta_vals), 3))

prev_frequencies = None

# ==========================================
# SOLVER CÚBICO CON SEGUIMIENTO DE HOJAS DE RIEMANN
# ==========================================
for i, Lambda in enumerate(lambda_vals):
    # Coeficientes exactos del polinomio cúbico de sistema abierto (Ecuación 7)
    coeffs = [
        -1j * De,
        1.0 + 2.0 * De * gamma_GW,
        1j * (De * (1.0 + gamma_GW**2) - Lambda + 2.0 * gamma_GW),
        -(1.0 + gamma_GW**2)
    ]

    # Raíces adimensionales de autovalor complejas Omega
    roots = np.roots(coeffs)
    # Conversión directa a frecuencias físicas de oscilación (Hz)
    frequencies = roots * f0

    # Algoritmo de clasificación por proximidad euclidiana (Evita saltos numéricos)
    if prev_frequencies is not None:
        idx_ordered = []
        available_indices = list(range(3))
        for p_f in prev_frequencies:
            distances = [np.abs(f - p_f) for f in frequencies]
            best_match = min(available_indices, key=lambda idx: distances[idx])
            idx_ordered.append(best_match)
            available_indices.remove(best_match)
        frequencies = frequencies[idx_ordered]

    prev_frequencies = frequencies
    re_roots[i, :] = np.real(frequencies)
    im_roots[i, :] = -np.imag(frequencies)  # Tasa de amortiguamiento -Im(f)

# ==========================================
# CONSTRUCCIÓN GRÁFICA BI-PANEL
# ==========================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.0), dpi=300)

# Umbrales críticos de colisión e inversión de estabilidad
zeta_overdamping = 3.30
zeta_collision = 4.50

# ----------- PANEL IZQUIERDO: FRECUENCIAS REALES -----------
ax1.plot(zeta_vals, re_roots[:, 0], color='darkblue', linewidth=1.8, label='Propagating Mode (+)')
ax1.plot(zeta_vals, re_roots[:, 1], color='cyan', linestyle='--', linewidth=1.5, label='Propagating Mode (-)')
ax1.plot(zeta_vals, re_roots[:, 2], color='darkgreen', linestyle=':', linewidth=1.5, label='Relaxation Mode')

ax1.axvline(zeta_overdamping, color='pink', linestyle='--', alpha=0.8, label='Overdamping Threshold (Closed Ref.)')
ax1.axvline(zeta_collision, color='grey', linestyle=':', alpha=0.8, label='Relaxation Collision (Closed Ref.)')

ax1.set_xlim(1.5, 4.7)
ax1.set_ylim(-3500, 3500)
ax1.grid(True, linestyle=':', alpha=0.5)
ax1.set_xlabel(r'Bulk Viscosity $\zeta$ ($10^{31}$ g $\mathrm{cm}^{-1}$ $\mathrm{s}^{-1}$)')
ax1.set_ylabel(r'Mode Frequency $\mathrm{Re}(f)$ (Hz)')
ax1.set_title('Avoided Crossing in Mode Frequencies')
ax1.legend(loc='upper right', fontsize=8, framealpha=0.9)

# Eje gemelo superior para Lambda
ax1_top = ax1.twiny()
ax1_top.set_xlim(ax1.get_xlim())
ax1_top.set_xticks([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
ax1_top.set_xticklabels(['0.53', '0.80', '1.06', '1.33', '1.60', '1.86', '2.13'])
ax1_top.set_xlabel(r'Dimensionless Viscous Coupling $\Lambda$')

# ----------- PANEL DERECHO: TASAS DE AMORTIGUAMIENTO -----------
ax2.plot(zeta_vals, im_roots[:, 0], color='orange', linewidth=1.8, label='Fluid Branch 1 (+)')
ax2.plot(zeta_vals, im_roots[:, 1], color='darkgreen', linestyle='--', linewidth=1.5, label='Fluid Branch 2 (-)')
ax2.plot(zeta_vals, im_roots[:, 2], color='purple', linewidth=1.8, label='Relaxation Branch')

ax2.axvline(zeta_overdamping, color='pink', linestyle='--', alpha=0.8)
ax2.axvline(zeta_collision, color='grey', linestyle=':', alpha=0.8)

ax2.set_xlim(1.5, 4.7)
ax2.set_ylim(-5000, 42000)
ax2.grid(True, linestyle=':', alpha=0.5)
ax2.set_xlabel(r'Bulk Viscosity $\zeta$ ($10^{31}$ g $\mathrm{cm}^{-1}$ $\mathrm{s}^{-1}$)')
ax2.set_ylabel(r'Decay Rate $-\mathrm{Im}(f)$ (Hz)')
ax2.set_title('Avoided Crossing in Decay Rates')
ax2.legend(loc='center left', fontsize=8, framealpha=0.9)

# Eje gemelo superior para Lambda
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
ax2_top.set_xticks([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
ax2_top.set_xticklabels(['0.53', '0.80', '1.06', '1.33', '1.60', '1.86', '2.13'])
ax2_top.set_xlabel(r'Dimensionless Viscous Coupling $\Lambda$')

plt.suptitle(r'Post-Merger Remnant QNM Evolution: Open-System Avoided Crossings ($\tau=0.005$ ms, $\gamma_{GW}=0.03$)', y=0.98)
plt.tight_layout()

# Guardar en alta calidad editorial
plt.savefig('fig4_avoided_crossings.pdf', format='pdf', bbox_inches='tight', dpi=300)
print("SUCCESS: Figure 4 exported to 'fig4_avoided_crossings.pdf'.")
plt.show()
