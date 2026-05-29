import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ==========================================
# 1. CONFIGURACIÓN EDITORIAL Y ENTORNO (CQG)
# ==========================================
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.titlesize": 13,
    "legend.fontsize": 9
})

# ==========================================
# 2. RIGOR MATEMÁTICO: DOMINIO Y ECUACIONES
# ==========================================
y = np.linspace(1.0, 6.0, 5000)

# Ecuaciones paramétricas del Discriminante de la familia cúbica (Ecuación 30)
De = (y**2 - 1) / (2 * y**3)
Lambda = -((y**2 + 1)**2) / (2 * y**3)

# Coordenadas exactas del Centro Organizador de Cúspide A3
De_crit = 1.0 / (3.0 * np.sqrt(3))
Lambda_crit = -8.0 / (3.0 * np.sqrt(3))

# ==========================================
# 3. ARQUITECTURA GRÁFICA Y LIENZO
# ==========================================
fig, ax = plt.subplots(figsize=(7.5, 5.0), dpi=300)

# ¡CORRECCIÓN CRÍTICA!: Colores pastel más definidos y alphas calibrados
# para que los cuadrantes de fondo sean perfectamente visibles
ax.axhspan(0, 3.5, color='#ffcccc', alpha=0.35)  # Hemisferio Físico (Rosa)
ax.axhspan(-2.5, 0, color='#cce5ff', alpha=0.35) # Hemisferio Virtual (Azul)

# Captura de líneas y puntos esenciales
ep_line, = ax.plot(De, Lambda, color='#1a1a1a', linestyle='--', linewidth=1.8)

cusp_point, = ax.plot(De_crit, Lambda_crit, color='#e63946', marker='o', markersize=7,
                      markeredgecolor='black', markeredgewidth=0.5, linestyle='None')

ns_point, = ax.plot(0, -2, color='#1d3557', marker='s', markersize=7,
                    markeredgecolor='black', markeredgewidth=0.5, linestyle='None')

# Líneas de guía cartesianas (el "ecuador" en Lambda = 0)
ax.axhline(0, color='#555555', linestyle='-', linewidth=0.8, alpha=0.6)
ax.axvline(0, color='#555555', linestyle='-', linewidth=0.8, alpha=0.6)

# ==========================================
# 4. CONSTRUCCIÓN DE LA LEYENDA (MÁXIMA SIMETRÍA)
# ==========================================
# Las cajitas de la leyenda usan exactamente el mismo color base sólido para mantener la armonía
patch_physical = mpatches.Patch(facecolor='#ffcccc', edgecolor='none')
patch_virtual = mpatches.Patch(facecolor='#cce5ff', edgecolor='none')

custom_handles = [patch_physical, patch_virtual, ep_line, cusp_point, ns_point]
custom_labels = [
    r'Physically Admissible Quadrant ($\Lambda > 0$)',
    r'Non-Physical / Virtual Quadrant ($\Lambda < 0$)',
    r'Imaginary axis EP2 branch ($\Delta = 0, \Lambda \leq 0$)',
    r'True $A_3$ Cusp Center ($De_*, \Lambda_*$) = (0.19, -1.54)',
    r'Navier-Stokes regular limit ($De = 0, \Lambda = -2$)'
]

ax.legend(handles=custom_handles, labels=custom_labels, loc='upper left',
          frameon=True, facecolor='white', framealpha=0.95, edgecolor='#dddddd')

# ==========================================
# 5. ROTULACIÓN Y LÍMITES ASINTÓTICOS
# ==========================================
ax.set_xlim(-0.02, 0.30)
ax.set_ylim(-2.5, 3.0)
ax.grid(True, linestyle=':', alpha=0.5, color='#cccccc')

ax.set_xlabel(r'Deborah number $De$')
ax.set_ylabel(r'Viscous coupling parameter $\Lambda$')
ax.set_title('Global Exceptional-Point Skeleton and the $A_3$ Organizing Center', pad=12)

plt.tight_layout()

# ==========================================
# 6. EXPORTACIÓN VECTORIAL DIRECTA
# ==========================================
plt.savefig('fig1_ep_curve.pdf', format='pdf', bbox_inches='tight', dpi=300)
print("SUCCESS: Plot generated with clearly visible shaded backgrounds.")
plt.show()
