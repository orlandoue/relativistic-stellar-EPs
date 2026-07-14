"""
=============================================================================
CQG-115250 -- nucleo de calculo certificado (CONVENCION CORRECTA DEL PAPER).
Open Cusp Polynomial Form (Ec. open_polynomial del manuscrito):
  P(Om;De,Lam,g) = -i De Om^3 + (1+2De g) Om^2
                   + i[De(1+g^2)+Lam+2g] Om - (1+g^2) = 0
Limite g=0: -i De Om^3 + Om^2 + i(De+Lam)Om - 1 = 0  (forma cerrada, Ec.cubic).
Verificado contra Tabla 1 a 5-6 cifras. Fuente unica de verdad de las figuras.
=============================================================================
"""
import numpy as np
from mpmath import mp, mpf, mpc, polyroots, findroot
mp.dps = 30

def coeffs(De, g, Lam):
    """[Om^3, Om^2, Om^1, Om^0] del Open Cusp Polynomial (numpy complex)."""
    return [-1j*De, (1+2*De*g), 1j*(De*(1+g**2)+Lam+2*g), -(1+g**2)]

def roots_np(De, g, Lam):
    return np.roots(coeffs(De, g, Lam))

def roots_mp(De, g, Lam):
    De,g,Lam = mpf(str(De)),mpf(str(g)),mpf(str(Lam))
    c=[-1j*De,(1+2*De*g),1j*(De*(1+g**2)+Lam+2*g),-(1+g**2)]
    return polyroots(c, maxsteps=300, extraprec=200)

def min_gap(De, g, Lam):
    r=roots_np(De,g,Lam)
    return min(abs(r[i]-r[j]) for i in range(3) for j in range(i+1,3))

def find_ep2(De, g, seed_Om=-1.2j, seed_L=1.5):
    De_,g_ = mpf(str(De)),mpf(str(g))
    def eqs(Om,Lam):
        c=[-1j*De_,(1+2*De_*g_),1j*(De_*(1+g_**2)+Lam+2*g_),-(1+g_**2)]
        P =c[0]*Om**3+c[1]*Om**2+c[2]*Om+c[3]
        Pp=3*c[0]*Om**2+2*c[1]*Om+c[2]
        return [P,Pp]
    Om,L=findroot(lambda a,b:eqs(a,b),(mpc(seed_Om),mpc(seed_L)))
    return complex(Om), complex(L)

def find_cusp(g):
    """cusp exacto (raiz triple P=P'=P''=0) a gamma dado."""
    g_=mpf(str(g))
    def eqs(De,Lam,Om):
        c=[-1j*De,(1+2*De*g_),1j*(De*(1+g_**2)+Lam+2*g_),-(1+g_**2)]
        P  =c[0]*Om**3+c[1]*Om**2+c[2]*Om+c[3]
        Pp =3*c[0]*Om**2+2*c[1]*Om+c[2]
        Ppp=6*c[0]*Om+2*c[1]
        return [P,Pp,Ppp]
    De,Lam,Om=findroot(lambda a,b,c:eqs(a,b,c),(mpf('0.19'),mpf('1.5'),mpc(0,-1.7)))
    return complex(Om), float(De.real), float(Lam.real)

def migration_De(g):
    """De* del cusp via ley de migracion (1+2De g)^3 = 27 De^2 (1+g^2)."""
    g_=mpf(str(g))
    return float(findroot(lambda D:(1+2*D*g_)**3-27*D**2*(1+g_**2), mpf('0.19')))

def ep2_curve(g, y_range=(1.02, 4.0), n=500):
    """(De,Lam) sobre la curva EP2, parametrizada por y=|Om| (Om=-i y, eje imag)."""
    De_l,Lam_l,y_l=[],[],[]
    g_=mpf(str(g))
    for y in np.linspace(*y_range,n):
        Om=mpc(0,-y)
        try:
            def eqs(De,Lam):
                c=[-1j*De,(1+2*De*g_),1j*(De*(1+g_**2)+Lam+2*g_),-(1+g_**2)]
                P =c[0]*Om**3+c[1]*Om**2+c[2]*Om+c[3]
                Pp=3*c[0]*Om**2+2*c[1]*Om+c[2]
                return [P,Pp]
            De,Lam=findroot(lambda a,b:eqs(a,b),(mpf('0.15'),mpf('1.5')))
            if 0<float(De.real)<1.5 and float(Lam.real)>0 and abs(De.imag)<1e-6:
                De_l.append(float(De.real)); Lam_l.append(float(Lam.real)); y_l.append(y)
        except Exception:
            pass
    return np.array(De_l),np.array(Lam_l),np.array(y_l)

if __name__=="__main__":
    print("=== verificacion contra Tabla 1 (convencion correcta) ===")
    Om,De,Lam=find_cusp(0.0)
    print(f"cusp g=0:   De*={De:.5f}(0.19245) Lam*={Lam:.5f}(1.53960) y*={abs(Om):.5f}(1.73205)")
    Om,De,Lam=find_cusp(0.1)
    print(f"cusp g=0.1: De*={De:.5f}(0.20329) Lam*={Lam:.5f}(1.37040) y*={abs(Om):.5f}(1.70634)")
    Om,L=find_ep2(0.15,0.1,-1.238j,1.51)
    print(f"EP2 0.15,0.1: Lam*={L.real:.5f}(1.50958) y*={abs(Om):.5f}(1.23851)")
    # migracion
    print(f"migracion g=0.1: De*={migration_De(0.1):.5f}(0.20329)")
