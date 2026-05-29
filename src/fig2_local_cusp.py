import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# CONFIGURACIÓN EDITORIAL ESTÁNDAR (CQG)
# ==========================================
plt.rcParams.update({
    "text.usetex": False,  # True si tu sistema local posee distribución LaTeX
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9
})

# ==========================================
# 1. RIGOR MATEMÁTICO: PARAMETRIZACIÓN DEL CUSP
# ==========================================
# Dividimos el dominio de 'y' en dos sectores respecto a y_crit = sqrt(3)
# para mapear de forma independiente ambas hojas tangenciales del despliegue
y_lower = np.linspace(1.0, np.sqrt(3), 2000)
y_upper = np.linspace(np.sqrt(3), 4.5, 2000)

# Ecuaciones paramétricas exactas del discriminante cúbico (Eq. 30)
De_lower = (y_lower**2 - 1) / (2 * y_lower**3)
Lambda_lower = -((y_lower**2 + 1)**2) / (2 * y_lower**3)

De_upper = (y_upper**2 - 1) / (2 * y_upper**3)
Lambda_upper = -((y_upper**2 + 1)**2) / (2 * y_upper**3)

# Coordenadas analíticas exactas de la singularidad A3
De_crit = 1.0 / (3.0 * np.sqrt(3))
Lambda_crit = -8.0 / (3.0 * np.sqrt(3))

# ==========================================
# 2. ARQUITECTURA GRÁFICA Y LIENZO
# ==========================================
fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=300)

# Pintamos el fondo para denotar que estamos dentro del cuadrante virtual/no físico
ax.set_facecolor('#f0f4f8')

# Graficación de las ramas del esqueleto topológico
line_trace, = ax.plot(De_lower, Lambda_lower, color='black', linestyle='-', linewidth=1.8)
line_cont, = ax.plot(De_upper, Lambda_upper, color='black', linestyle='-', linewidth=1.8)

# Punto crítico: Centro organizador A3 (Cúspide triple defectiva)
cusp_point, = ax.plot(De_crit, Lambda_crit, color='#e63946', marker='o', markersize=7,
                      markeredgecolor='black', markeredgewidth=0.5, linestyle='None')

# ==========================================
# 3. ROTULACIÓN CIENTÍFICA EN INGLÉS CORREGIDA
# ==========================================
ax.set_xlim(0.12, 0.24)
ax.set_ylim(-1.8, -1.3)
ax.grid(True, linestyle=':', alpha=0.5, color='#cccccc')

ax.set_xlabel(r'Deborah number $De$')
ax.set_ylabel(r'Viscous coupling parameter $\Lambda$')
ax.set_title(r'Local Geometry of the $A_3$ Cusp Organizing Center', pad=12)

# Corrección de rigurosidad: Identificación exacta de las hojas discriminantes
custom_handles = [line_cont, line_trace, cusp_point]
custom_labels = [
    r'Cusp branch (upper sheet, $\Delta = 0$)',
    r'Cusp branch (lower sheet, $\Delta = 0$)',
    r'True $A_3$ Cusp Center'
]
ax.legend(handles=custom_handles, labels=custom_labels, loc='upper left',
          frameon=True, facecolor='white', framealpha=0.95, edgecolor='#dddddd')

plt.tight_layout()

# Exportación vectorial directa sin pérdidas
plt.savefig('fig2_local_cusp.pdf', format='pdf', bbox_inches='tight', dpi=300)
print("SUCCESS: Figure 2 exported to 'fig2_local_cusp.pdf'.")
plt.show()
