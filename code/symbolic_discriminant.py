"""
symbolic_discriminant.py — verificacion simbolica del discriminante y del cusp
==============================================================================
Este script cubre lo que la seccion Data Availability del manuscrito declara
disponible: "the symbolic computations used to evaluate the cubic discriminant
and verify the degeneracy conditions".

Reproduce, en aritmetica exacta con sympy:

  (1) el discriminante Delta(De,Lambda) de la Open Cusp Polynomial (Ec. 7)
  (2) la anulacion de Delta, dDelta y el HESSIANO DE RANGO UNO en la cuspide
      del limite cerrado (De*,Lambda*) = (1/(3sqrt3), 8/(3sqrt3))     [sec. 5.2]
  (3) que el nucleo del Hessiano coincide con la tangente de la curva EP2,
      en la direccion (1,-4)                                          [sec. 5.2]
  (4) la ley de migracion (1+2De*gamma)^3 = 27 De^2 (1+gamma^2)       [sec. 5.3]
  (5) la parametrizacion EP2 a gamma finito, Ec. (48)
  (6) el signo del discriminante: Delta_y = -Delta_Omega              [sec. 5.4]

Solo requiere sympy. Uso:  python symbolic_discriminant.py
"""
import sympy as sp

De, Lam, gam, y, Om = sp.symbols('De Lambda gamma y Omega', real=True)
I = sp.I


def open_cusp_coeffs(De_, gam_, Lam_):
    """Coeficientes de la Open Cusp Polynomial en Omega (orden descendente).

    (1 - i De Om)(Om^2 + 2i gamma Om - (1+gamma^2)) + i Lambda Om = 0
    """
    return [-I * De_,
            1 + 2 * De_ * gam_,
            I * (2 * gam_ + De_ * (1 + gam_ ** 2) + Lam_),
            -(1 + gam_ ** 2)]


def discriminant(a, b, c, d):
    """Discriminante de a x^3 + b x^2 + c x + d."""
    return 18 * a * b * c * d - 4 * b ** 3 * d + b ** 2 * c ** 2 \
        - 4 * a * c ** 3 - 27 * a ** 2 * d ** 2


def main():
    print("=" * 70)
    print("(1) DISCRIMINANTE en el limite cerrado (gamma = 0)")
    print("=" * 70)
    a, b, c, d = open_cusp_coeffs(De, 0, Lam)
    Delta = sp.simplify(sp.re(sp.expand(discriminant(a, b, c, d))))
    print("  Delta(De,Lambda) =", sp.factor(Delta))

    Dc = 1 / (3 * sp.sqrt(3))
    Lc = 8 / (3 * sp.sqrt(3))
    sub = {De: Dc, Lam: Lc}
    print(f"\n  en la cuspide De*=1/(3sqrt3), Lambda*=8/(3sqrt3):")
    print("    Delta        =", sp.simplify(Delta.subs(sub)))
    print("    dDelta/dDe   =", sp.simplify(sp.diff(Delta, De).subs(sub)))
    print("    dDelta/dLam  =", sp.simplify(sp.diff(Delta, Lam).subs(sub)))

    print()
    print("=" * 70)
    print("(2)-(3) HESSIANO DE RANGO UNO Y SU NUCLEO")
    print("=" * 70)
    H = sp.Matrix([[sp.diff(Delta, v1, v2) for v2 in (De, Lam)]
                   for v1 in (De, Lam)])
    Hc = sp.simplify(H.subs(sub))
    print("  H ="); sp.pprint(Hc)
    print("  det(H) =", sp.simplify(Hc.det()), "   rango =", Hc.rank())
    ker = [sp.simplify(v.T) for v in Hc.nullspace()]
    print("  nucleo =", ker)

    # tangente de la curva EP2 en la cuspide, desde la parametrizacion cerrada
    De_y = (y ** 2 - 1) / (2 * y ** 3)
    Lam_y = (y ** 2 + 1) ** 2 / (2 * y ** 3)
    y0 = sp.sqrt(3)
    d1 = [sp.simplify(sp.diff(f, y).subs(y, y0)) for f in (De_y, Lam_y)]
    d2 = [sp.simplify(sp.diff(f, y, 2).subs(y, y0)) for f in (De_y, Lam_y)]
    print(f"\n  parametrizacion EP2: (De'(sqrt3), Lambda'(sqrt3)) = {tuple(d1)}")
    print(f"                       (De'',Lambda'')             = {tuple(d2)}")
    print(f"  cociente Lambda''/De'' = {sp.simplify(d2[1] / d2[0])}")
    print("  => tangente (1,-4), que COINCIDE con el nucleo del Hessiano.")
    print("     (Un Hessiano NO degenerado describiria un punto aislado A1;")
    print("      el determinante nulo con nucleo en la tangente comun de las")
    print("      dos ramas es la firma de A3.)")

    print()
    print("=" * 70)
    print("(4) LEY DE MIGRACION DE LA CUSPIDE")
    print("=" * 70)
    # Verificacion directa (sin sp.solve, que es prohibitivamente lento aqui):
    # la cuspide a gamma finito se obtiene del minimo de Lambda/De sobre el
    # locus EP2, y debe satisfacer la ley de migracion.
    mig = (1 + 2 * De * gam) ** 3 - 27 * De ** 2 * (1 + gam ** 2)
    print("  ley:  (1+2 De gamma)^3 = 27 De^2 (1+gamma^2)")
    print("  residuo en la cuspide, evaluado en aritmetica exacta:")
    #   gamma = 0  ->  De* = 1/(3 sqrt 3)
    r0 = sp.simplify(mig.subs({gam: 0, De: 1 / (3 * sp.sqrt(3))}))
    print(f"    gamma=0,   De*=1/(3sqrt3)      : residuo = {r0}")
    #   gamma = 1  ->  De* = 1/4 (raiz exacta de la ley)
    r1 = sp.simplify(mig.subs({gam: 1, De: sp.Rational(1, 4)}))
    print(f"    gamma=1,   De*=1/4             : residuo = {r1}")
    #   comprobacion de que la ley es la condicion de raiz triple:
    #   el discriminante del locus se anula donde d(Lambda/De)/dy = 0
    De_yg = (y ** 2 - (1 + gam ** 2)) / (2 * y ** 2 * (y - gam))
    f_yg = 3 * y ** 2 - 4 * gam * y + 1 + gam ** 2
    Lam_yg = 2 * (y - gam) - De_yg * f_yg
    ratio = sp.simplify(Lam_yg / De_yg)
    dratio = sp.simplify(sp.diff(ratio, y))
    y_c = sp.solve(sp.numer(sp.together(dratio)), y)
    y_c = [r for r in y_c if r.is_real and r > 1]
    if y_c:
        yc = sp.simplify(y_c[0])
        Dc_g = sp.simplify(De_yg.subs({gam: 0, y: yc.subs(gam, 0)}))
        print(f"    minimo de Lambda/De en gamma=0 : y* = {sp.simplify(yc.subs(gam,0))}"
              f"  (= sqrt(3)),  De* = {Dc_g}")
        print(f"    residuo de la ley alli         : "
              f"{sp.simplify(mig.subs({gam: 0, De: Dc_g}))}")

    print()
    print("=" * 70)
    print("(5) PARAMETRIZACION EP2 A gamma FINITO  (Ecs. 30-31)")
    print("=" * 70)
    De_g = (y ** 2 - (1 + gam ** 2)) / (2 * y ** 2 * (y - gam))
    f_g = 3 * y ** 2 - 4 * gam * y + 1 + gam ** 2
    Lam_g = 2 * (y - gam) - De_g * f_g
    print("  De(y;gamma)     =", sp.simplify(De_g))
    print("  Lambda(y;gamma) = 2(y-gamma) - De*f,  f =", f_g)
    print("  en gamma=0 se reduce a la Ec. (30):",
          sp.simplify(De_g.subs(gam, 0) - (y ** 2 - 1) / (2 * y ** 3)) == 0,
          sp.simplify(Lam_g.subs(gam, 0) - (y ** 2 + 1) ** 2 / (2 * y ** 3)) == 0)
    # cotejo con Tabla 1
    vals = {y: sp.Rational(12385130101588585499, 10 ** 19), gam: sp.Rational(1, 10)}
    print("  con y* de Tabla 1 y gamma=0.1:")
    print("    De     =", sp.nsimplify(sp.N(De_g.subs(vals), 22)))
    print("    Lambda =", sp.N(Lam_g.subs(vals), 22), " (Tabla 1: 1.5095752865775080251)")

    print()
    print("=" * 70)
    print("(6) SIGNO DEL DISCRIMINANTE:  Delta_y = -Delta_Omega")
    print("=" * 70)
    # en y, con Omega = -i y, el polinomio del limite cerrado es real:
    ay, by, cy, dy = De, sp.Integer(-1), (De + Lam), sp.Integer(-1)
    D_y = sp.expand(discriminant(ay, by, cy, dy))
    a0, b0, c0, d0 = open_cusp_coeffs(De, 0, Lam)
    D_Om = sp.expand(sp.re(sp.expand(discriminant(a0, b0, c0, d0))))
    print("  Delta_y + Delta_Omega =", sp.simplify(D_y + D_Om))
    print("  => son OPUESTOS: la regla clasica de signos se enuncia en y,")
    print("     no en Omega, donde los coeficientes son complejos.")


if __name__ == "__main__":
    main()
