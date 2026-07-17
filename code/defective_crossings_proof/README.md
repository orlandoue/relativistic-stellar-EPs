# ¿Son "diabólicos" los cruces intersectoriales del paper? — No: son defectivos

**Manuscrito auditado:** CQG-115250, *"Exceptional points and cusp singularities
in relativistic stellar oscillations with causal bulk relaxation"*, Sec. 7.4(iii)
y Ec. (44).

## 1. Qué afirma el paper y por qué importa

En la Sec. 7.4(iii), al describir la deformación del perfil de transporte
ζ(x) = ζ₀ sech²x (1+η sech²x), el paper señala que en **η = 0** (sistema
exactamente resoluble, perfil alineado) dos sectores cualesquiera de igual
paridad m, M comparten un autovalor en

```
s = -(m+M+1)/2 ,   ζ₀ = (1+sτ)[s²+(2m+1)s+m(m+1)+V0] / (-s)
```

y los llama **"a semisimple (diabolic) degeneracy"**. Tres miembros de esta
familia se usan como puntos de nacimiento de las ramas EP2 de la Figura 4:
`(2,10)`, `(0,12)`, `(0,16)`.

"Diabólico" es un término técnico preciso (von Neumann–Wigner, Berry): un
cruce **semisimple**, donde la multiplicidad algebraica y la geométrica del
autovalor doble coinciden (=2) — es decir, existen **dos** autovectores
genuinamente independientes que se cruzan sin interactuar. Es la contraparte
"benigna" de un punto excepcional (EP2), que en cambio es **defectivo**:
multiplicidad algebraica 2 pero geométrica 1 (un solo autovector, bloque de
Jordan, el operador no diagonaliza ahí).

La pregunta que motivó esta demostración: ¿el pencil hiperboloidal del paper,
que **no** es diagonal sino solo *triangular* en la base de monomios σᵐ
(Teorema de la Sec. 7.2), realmente deja dos autovectores independientes en
esos cruces, o el acoplamiento fuera de la diagonal los colapsa a uno solo?

## 2. Resultado

**Los tres puntos son defectivos (bloque de Jordan), no diabólicos.**
Multiplicidad algebraica = 2 (por construcción: dos factores del determinante
se anulan a la vez), pero **multiplicidad geométrica = 1** en los tres casos,
verificado de tres formas independientes (aritmética racional exacta,
insensible a la truncación; y SVD del operador discretizado completo, tal
como lo construye el propio repositorio del autor).

| Punto | s* | ζ₀* | mult. algebraica | mult. geométrica | Soporte del autovector |
|---|---|---|---|---|---|
| (2,10) | −6.5 | 6601/1300 ≈ 5.077692 | 2 | **1** | σ⁰, σ² (sector m=2) |
| (0,12) | −6.5 | 9821/1300 ≈ 7.554615 | 2 | **1** | σ⁰ (sector m=0) |
| (0,16) | −8.5 | 13261/1700 ≈ 7.800588 | 2 | **1** | σ⁰ (sector m=0) |

## 3. Por qué ocurre (el mecanismo, en una idea)

El script `01` deriva de forma simbólica que

```
T(s) σᵐ = P_m(s) σᵐ  −  (1+sτ) m(m−1) σᵐ⁻²
```

es decir, T(s) en la base {1, σ², σ⁴, …} es una matriz **triangular superior
de ancho de banda 1** (bidiagonal): la diagonal son las cúbicas sectoriales
P_m(s) (esto es lo que ya prueba el Teorema del paper, y es correcto — el
**espectro**, es decir el conjunto de raíces de det T(s) = ∏ₘ P_m(s), sí es
la unión exacta de las raíces sectoriales). Pero el coeficiente fuera de la
diagonal, −(1+sτ)m(m−1), es **genéricamente distinto de cero** para todo m
par ≥ 2. Esa cadena de acoplamiento ininterrumpida es la que decide la
estructura de **autovectores**, no solo de autovalores.

Al resolver T(s*)v = 0 recorriendo la recursión bidiagonal de arriba (grado
alto) hacia abajo: el cero de **mayor** grado (M=10, 12 ó 16) nunca ancla una
dirección libre genuina, porque la cadena que desciende desde él termina
forzando una condición de consistencia en el peldaño del cero de **menor**
grado (m=2 ó 0) que sólo se satisface si esa dirección entera se anula. El
único autovector que sobrevive vive exclusivamente en el sector de **menor**
grado. El cero del sector de mayor grado es, en este sentido preciso, un
"cero fantasma" del determinante: contribuye a la multiplicidad algebraica
pero no a la geométrica.

Esto no es un artefacto de truncar la base en algún Mmax arbitrario: el
script `02` repite el cálculo de rango exacto (aritmética racional) con
Mmax = 20, 30, 40, 60 y la nulidad da **1** en los cuatro casos, para los
tres puntos. El script `03` repite la verificación a nivel del operador
completo (discretización de Chebyshev, exactamente como la construye
`bordered_newton_robustness.py` del repositorio del autor para certificar
sus propios EP2/EP3): el segundo valor singular más pequeño queda resuelto
~10⁶–10⁷ veces por encima del piso de precisión doble, confirmando de forma
independiente que sólo hay **una** dirección nula, no dos.

## 4. Por qué esto no debilita al paper — al contrario

Este hallazgo **no** contradice ningún resultado físico o cuantitativo del
paper: η_c=0.19415, la identidad de la rama A′ (nace exactamente en
`(−6.5, 6601/1300)`, confirmado por continuación independiente), el escape
complejo a η=0.21 y η=0.40, y los exponentes de *splitting* certificados,
fueron todos reproducidos de forma independiente en la auditoría previa y
siguen siendo correctos.

De hecho, el resultado aquí es **consistente con, y refuerza**, la física del
paper: si los cruces fuesen genuinamente diabólicos (semisimples), no habría
ninguna garantía de que una perturbación genérica η>0 genere inmediatamente
una rama EP2 (√ε) a partir de ellos — un cruce semisimple típicamente se
abre linealmente en ε bajo perturbación genérica. Que el paper *sí* observe
el nacimiento de ramas EP2 genuinas (B, R2) exactamente en estos puntos es
precisamente lo que se espera si el punto de partida ya era defectivo (como
aquí se demuestra), no diabólico.

**Corrección sugerida (cosmética, no de fondo):** en el abstract, Ec. (44),
Sec. 7.4(iii) y el pie de la Figura 4, sustituir "a semisimple (diabolic)
degeneracy" por algo del tipo *"an algebraically double, geometrically
simple inter-sector crossing, at which the two sector eigenvalue branches
cross with vanishing splitting to leading order but only one physical
eigenmode"* — o, más brevemente, "a defective inter-sector crossing". Ningún
número del paper cambia.

## 5. Cómo correr los scripts

```bash
# Script 1 y 2: solo necesitan sympy, no requieren el repo del autor
pip install sympy --break-system-packages
python3 01_estructura_triangular_simbolica.py
python3 02_multiplicidad_geometrica_racional.py

# Script 3: requiere el repositorio del autor y mpmath/numpy
git clone https://github.com/orlandoue/relativistic-stellar-EPs
pip install mpmath numpy --break-system-packages
REPO_PATH=./relativistic-stellar-EPs python3 03_verificacion_operador_chebyshev.py
```

- **`01_estructura_triangular_simbolica.py`** — deriva simbólicamente (sympy,
  sin ningún redondeo) la acción de T(s) sobre σᵐ y aísla el coeficiente de
  acoplamiento fuera de la diagonal. Es el fundamento de todo lo demás.
- **`02_multiplicidad_geometrica_racional.py`** — el cálculo central: rango y
  espacio nulo *exactos* (aritmética racional) de T(s*,ζ₀*) para los tres
  puntos declarados, con verificación de estabilidad ante la truncación y
  extracción explícita del autovector superviviente.
- **`03_verificacion_operador_chebyshev.py`** — confirmación independiente a
  nivel del operador discretizado completo (no ya la matriz de monomios
  truncada), reutilizando literalmente el código de discretización del
  repositorio del autor.
