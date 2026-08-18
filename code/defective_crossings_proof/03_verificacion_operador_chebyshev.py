"""
03_verificacion_operador_chebyshev.py
========================================
PASO 3 de la demostracion: CONFIRMACION INDEPENDIENTE, a nivel del operador
COMPLETO (no ya la matriz truncada a monomios del script 02, sino la
discretizacion por colocacion de Chebyshev en sigma in [-1,1] usada para
certificar los EP2/EP3 del articulo).

Reutiliza setup_mp / T_all de ../bordered_newton_robustness.py: mismo codigo,
mismos operadores L1, L2, C y misma discretizacion, sin ningun atajo propio.

Metodo: para cada punto defectivo declarado, se arma la matriz operador
T(s*, zeta0*) discretizada (N=20-22 puntos de Chebyshev) EN PRECISION
MULTIPLE (mpmath, dps=25) y se calcula su descomposicion en valores
singulares (SVD, numpy en doble precision sobre la matriz ya evaluada).
Si la degeneracion fuese diabolica (semisimple, geometricamente doble), DOS
valores singulares caerian a cero de maquina; si es defectiva, solo UNO.
"""
import sys, os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
try:
    from bordered_newton_robustness import setup_mp, T_all
except ImportError:
    raise SystemExit(
        "Could not import bordered_newton_robustness. Run this script from\n"
        "code/defective_crossings_proof/ inside the repository tree."
    )

from mpmath import mp, mpf
import numpy as np

mp.dps = 25
V0 = mpf('25.25')   # k=5,  fiducial del paper (Tabla 1 / Sec. 7.3)
tau = mpf('0.03')   # De=0.15 fiducial

PUNTOS = {
    '(2,10)': (mpf('-6.5'),  mpf(6601) / 1300, 20),
    '(0,12)': (mpf('-6.5'),  mpf(9821) / 1300, 20),
    '(0,16)': (mpf('-8.5'), mpf(13261) / 1700, 22),
}

print("="*78)
print("Multiplicidad geometrica a nivel del OPERADOR COMPLETO (Chebyshev)")
print("="*78)
for etiqueta, (s, z0, N) in PUNTOS.items():
    n = N + 1
    L1, L2, Cp, W, nn = setup_mp(N, V0, mpf('0'))   # eta=0, perfil alineado
    T, T1, T2 = T_all(s, z0, tau, L1, L2, Cp, n)
    Tnp = np.array([[complex(T[i, j]) for j in range(n)] for i in range(n)])
    sv = np.linalg.svd(Tnp, compute_uv=False)

    print(f"\n{etiqueta}:  s*={float(s):.6f}  zeta0*={float(z0):.6f}  "
          f"(discretizacion N={N}, {n} nodos de Chebyshev)")
    print("  5 menores valores singulares de T(s*):")
    for x in sv[-5:]:
        print(f"    {x:.4e}   (razon a sv_max = {x/sv[0]:.3e})")

    # el segundo menor valor singular debe estar CLARAMENTE resuelto (no en
    # cero de maquina) si la degeneracion es defectiva (mult. geom. = 1)
    ratio1, ratio2 = sv[-1]/sv[0], sv[-2]/sv[0]
    # sv[-1] cae al piso de precision doble (~1e-16); sv[-2] debe estar muy
    # por encima de ese piso para contar como "claramente resuelto, no nulo"
    if ratio1 < 1e-13 and ratio2 > 1e-11:
        veredicto = (f"DEFECTIVO: solo 1 valor singular en cero de maquina "
                      f"(el segundo esta {ratio2/ratio1:.1e} veces por encima del piso)")
    elif ratio1 < 1e-13 and ratio2 < 1e-13:
        veredicto = "AMBIGUO a esta N/precision: subir N o dps"
    else:
        veredicto = "revisar (ningun valor singular claramente nulo)"
    print(f"  ==> {veredicto}")

print("""
Nota metodologica: esta es una confirmacion NUMERICA (doble precision sobre
la matriz ya evaluada en multiprecision) del resultado EXACTO y auto-
consistente-con-cualquier-truncacion del script 02. La evidencia primaria y
rigurosa es la de aritmetica racional exacta (02); este script demuestra que
el mismo patron persiste en el operador tal como se lo discretiza en este
repositorio para certificar los EP2/EP3 del articulo.
""")
