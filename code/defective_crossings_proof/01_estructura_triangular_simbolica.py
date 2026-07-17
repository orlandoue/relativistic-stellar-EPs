"""
01_estructura_triangular_simbolica.py
======================================
PASO 1 de la demostracion: establecer, de forma simbolica y exacta, la accion
del pencil hiperboloidal T(s) sobre los monomios sigma^m.

Contexto (paper, Sec. 7.2, Teorema "Sector tower"):
    T(s) phi = [(1+s*tau)(s^2 - s*L2 - L1) + zeta0*s*C] phi = 0
    L1 = (1-sigma^2) d^2/dsigma^2 - 2 sigma d/dsigma - V0
    L2 = -2 sigma d/dsigma - 1

El teorema del paper prueba que, en la base {sigma^m}, T(s) es TRIANGULAR y que
su determinante (para todo N) es el producto de las cubicas sectoriales
P_m(s) = (1+s tau)(s^2+(2m+1)s+m(m+1)+V0) + zeta0*s.

Ese resultado (ya verificado independientemente en la auditoria previa) NO basta
para decidir si un cruce accidental P_m(s*)=P_M(s*)=0 (m != M) es un genuino
cruce "diabolico" (dos autovectores independientes) o una degeneracion
"defectiva" (un solo autovector, bloque de Jordan). Eso depende de los
terminos FUERA de la diagonal, que aqui se derivan explicitamente.

Este script deriva, de forma puramente simbolica (sympy, sin ningun redondeo):
  (a) L1 sigma^m y L2 sigma^m exactos,
  (b) T sigma^m expresado en la base {sigma^m, sigma^(m-2)},
  (c) la identificacion de los coeficientes diagonal y fuera-de-diagonal,
que son el insumo del script 02.
"""
import sympy as sp

sigma, s, tau, zeta0, V0, m = sp.symbols('sigma s tau zeta0 V0 m', )

# ---- operadores L1, L2 actuando sobre sigma^m (derivacion directa) ----
f = sigma**m
L2_f = sp.expand(-2*sigma*sp.diff(f, sigma) - f)
L1_f = sp.expand((1-sigma**2)*sp.diff(f, sigma, 2) - 2*sigma*sp.diff(f, sigma) - V0*f)

print("="*78)
print("(a) Accion de L1, L2 sobre sigma^m")
print("="*78)
print("L2 sigma^m =", sp.simplify(L2_f), "   [esperado: -(2m+1) sigma^m]")
print("L1 sigma^m =", sp.simplify(L1_f),
      "   [esperado: m(m-1) sigma^(m-2) - (m(m+1)+V0) sigma^m]")

check_L2 = sp.simplify(L2_f - (-(2*m+1)*sigma**m))
check_L1 = sp.simplify(L1_f - (m*(m-1)*sigma**(m-2) - (m*(m+1)+V0)*sigma**m))
print("residuo L2:", check_L2, "   residuo L1:", check_L1)

print()
print("="*78)
print("(b) T(s) sigma^m = (1+s tau)[s^2 - s L2 - L1] sigma^m + zeta0 s sigma^m")
print("    (perfil alineado C=1, eta=0; el caso eta>0 se trata en la Sec. 7.4")
print("    del paper y no es necesario para localizar los puntos 'diabolicos')")
print("="*78)
T_f = (1+s*tau)*(s**2*sigma**m - s*(-(2*m+1)*sigma**m) - (m*(m-1)*sigma**(m-2) - (m*(m+1)+V0)*sigma**m)) \
      + zeta0*s*sigma**m
T_f = sp.expand(T_f)
coef_diag = sp.simplify(T_f.coeff(sigma**m))
coef_off  = sp.simplify(T_f.coeff(sigma**(m-2)))
print("Coeficiente de sigma^m       (diagonal)       =", sp.factor(coef_diag))
print("  [debe ser exactamente P_m(s) = (1+s tau)(s^2+(2m+1)s+m(m+1)+V0)+zeta0 s]")
Pm_formula = (1+s*tau)*(s**2+(2*m+1)*s+m*(m+1)+V0) + zeta0*s
print("  residuo vs P_m(s):", sp.simplify(coef_diag - Pm_formula))
print()
print("Coeficiente de sigma^(m-2)   (fuera de diagonal, ACOPLAMIENTO) =", coef_off)
print("  [debe ser exactamente  -(1+s tau) m(m-1) ]")
print("  residuo:", sp.simplify(coef_off - (-(1+s*tau)*m*(m-1))))

print()
print("="*78)
print("CONCLUSION ESTRUCTURAL (insumo del script 02)")
print("="*78)
print("""
T sigma^m = P_m(s) sigma^m  -  (1+s tau) m(m-1) sigma^(m-2)

=> En la base ordenada por grado creciente {1, sigma^2, sigma^4, ...}, la matriz
   de T(s) (restringida a potencias pares) es TRIANGULAR SUPERIOR CON ANCHO DE
   BANDA 1 (bidiagonal): la unica entrada fuera de la diagonal conecta el grado
   m con el grado m-2, con coeficiente -(1+s tau) m(m-1), que es DISTINTO DE
   CERO para todo m par >= 2 (salvo en la resonancia espuria s=-1/tau, ajena a
   los puntos de cruce del paper).

   Esta cadena de acoplamiento ininterrumpida es la que decide, en el script
   02, si un cruce accidental P_m(s*)=P_M(s*)=0 dejara UNO o DOS autovectores
   independientes.
""")
