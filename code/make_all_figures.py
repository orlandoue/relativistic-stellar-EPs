"""
Regenerate all six figures of the manuscript from the certified core (cqg_core.py).
Every figure derives from the single validated Open Cusp Polynomial convention.
Usage: python make_all_figures.py   (outputs to ../figures/)
"""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpmath import mp, mpf, mpc, findroot
from cqg_core import find_cusp, find_ep2, roots_np, roots_mp
mp.dps = 40
OUT = "../figures/"

def ep2_closed(y):
    """EP2 of the closed (gamma=0) pencil, root on imaginary axis at Om=-i y."""
    y_ = mpf(str(y)); Om = mpc(0, -y_)
    def eqs(De, Lam):
        P  = -1j*De*Om**3 + Om**2 + 1j*(De+Lam)*Om - 1
        Pp = -3j*De*Om**2 + 2*Om + 1j*(De+Lam)
        return [P, Pp]
    try:
        De, Lam = findroot(lambda a,b: eqs(a,b), (mpf('0.17'), mpf('1.5')))
        if abs(De.imag) < 1e-8: return float(De.real), float(Lam.real)
    except Exception: pass
    return None

# ---------- FIG 1: global EP2 curve ----------
def fig1():
    ys = np.linspace(1.02, 4.5, 300)
    pts = [ep2_closed(y) for y in ys]; pts = [p for p in pts if p and 0<=p[0]<1.2 and p[1]>0]
    pts = np.array(pts); idx = np.argsort(pts[:,0])
    Om,Dc,Lc = find_cusp(0.0); Om1,Dc1,Lc1 = find_cusp(0.1)
    fig,ax = plt.subplots(figsize=(6.5,5))
    ax.plot(pts[idx,0],pts[idx,1],'-',color='C0',lw=2,label=r'EP2 curve ($\Delta=0$)')
    ax.plot(0,2,'s',color='C2',ms=10,label='Navier-Stokes EP2 $(0,2)$')
    ax.plot(Dc,Lc,'*',color='red',ms=18,label=f'$A_3$ cusp $\\gamma=0$')
    ax.plot(Dc1,Lc1,'D',color='darkred',ms=8,label=f'migrated cusp $\\gamma=0.1$')
    ax.set_xlabel('Deborah number $De$'); ax.set_ylabel(r'Viscous coupling $\Lambda$')
    ax.set_title('Global exceptional-point (EP2) curve'); ax.legend(fontsize=9); ax.grid(alpha=0.3)
    ax.set_xlim(-0.02,0.65); ax.set_ylim(0,3.2)
    fig.tight_layout(); fig.savefig(OUT+'fig1_ep_curve.pdf',dpi=150,bbox_inches='tight'); plt.close()

# ---------- FIG 2: cusp zoom with discriminant regions ----------
def fig2():
    Om,Dc,Lc = find_cusp(0.0); ystar = abs(Om)
    ys = np.linspace(float(ystar)-0.6, float(ystar)+0.6, 400)
    pts = [ep2_closed(y) for y in ys]; pts=[p for p in pts if p and p[1]>0]; pts=np.array(pts)
    DD=np.linspace(0.12,0.28,200); LL=np.linspace(1.2,1.9,200); DEg,LAg=np.meshgrid(DD,LL)
    nreal=np.zeros_like(DEg)
    for i in range(DEg.shape[0]):
        for j in range(DEg.shape[1]):
            r=roots_np(DEg[i,j],0.0,LAg[i,j]); nreal[i,j]=np.sum(np.abs(r.imag)<0.05)
    fig,ax=plt.subplots(figsize=(6.5,5.5))
    ax.contourf(DEg,LAg,nreal,levels=[0.5,1.5,2.5,3.5],colors=['#e8f0ff','#ffe8e8','#e8ffe8'],alpha=0.5)
    idx=np.argsort(pts[:,0]); ax.plot(pts[idx,0],pts[idx,1],'-',color='C0',lw=2,label=r'EP2 locus ($\Delta=0$)')
    ax.plot(Dc,Lc,'*',color='red',ms=20,label=f'$A_3$ cusp')
    ax.set_xlabel('Deborah number $De$'); ax.set_ylabel(r'Viscous coupling $\Lambda$')
    ax.set_title(r'Local structure near the $A_3$ cusp'); ax.legend(fontsize=10)
    ax.set_xlim(0.12,0.28); ax.set_ylim(1.2,1.9); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(OUT+'fig2_cusp_zoom.pdf',dpi=150,bbox_inches='tight'); plt.close()

# ---------- FIG 3 & 6: certified splitting over the certified range ----------
def fig_splitting(fname, title):
    De_ep=0.15; g=0.1
    Om_ep,L_ep=find_ep2(De_ep,g,-1.238j,1.51); Lstar=L_ep.real
    Om_c,Dcusp,Lcusp=find_cusp(0.1)
    def se(eps):
        r=[complex(x) for x in roots_mp(De_ep,g,Lstar*(1+eps))]
        return sorted(abs(r[i]-r[j]) for i in range(3) for j in range(i+1,3))[0]
    def sc(eps):
        r=[complex(x) for x in roots_mp(Dcusp,g,Lcusp*(1+eps))]
        return sorted(abs(r[i]-r[j]) for i in range(3) for j in range(i+1,3))[-1]
    # CERTIFIED RANGE: eps in [1e-10, 1e-2] (eight decades, matches Table 1)
    epss=np.logspace(-10,-2,25)
    yE=np.array([float(se(e)) for e in epss]); yC=np.array([float(sc(e)) for e in epss])
    # local exponent in the asymptotic regime only
    mask=epss<=1e-4
    expE=np.polyfit(np.log(epss[mask]),np.log(yE[mask]),1)[0]
    expC=np.polyfit(np.log(epss[mask]),np.log(yC[mask]),1)[0]
    fig,ax=plt.subplots(figsize=(6.5,5))
    ax.loglog(epss,yE,'o',color='C0',ms=6,label=f'EP2 (slope {expE:.4f}, $\\to 1/2$)')
    ax.loglog(epss,yC,'s',mfc='none',mec='red',ms=7,label=f'$A_3$ cusp (slope {expC:.4f}, $\\to 1/3$)')
    ax.loglog(epss,yE[0]*(epss/epss[0])**0.5,'--',color='C0',alpha=0.5,lw=1,label=r'$\epsilon^{1/2}$')
    ax.loglog(epss,yC[0]*(epss/epss[0])**(1/3),'--',color='red',alpha=0.5,lw=1,label=r'$\epsilon^{1/3}$')
    ax.set_xlabel(r'relative detuning $|\epsilon|$'); ax.set_ylabel(r'eigenvalue splitting $|\delta\Omega|$')
    ax.set_title(title); ax.legend(fontsize=9,loc='lower right'); ax.grid(alpha=0.3,which='both')
    fig.tight_layout(); fig.savefig(OUT+fname,dpi=150,bbox_inches='tight'); plt.close()
    return expE, expC

if __name__ == "__main__":
    print("Generating figures from certified core...")
    fig1(); print("  fig1_ep_curve.pdf")
    fig2(); print("  fig2_cusp_zoom.pdf")
    eE,eC = fig_splitting('fig3_splitting.pdf','Non-analytic scaling of the eigenvalue splitting')
    print(f"  fig3_splitting.pdf  (asymptotic exponents: EP2={eE:.5f}, cusp={eC:.5f})")
    fig_splitting('fig6_certified_scaling.pdf','Certified non-analytic splitting at the degeneracies')
    print("  fig6_certified_scaling.pdf")
    print("Note: fig4 (post-merger) and fig5 (robustness) require the")
    print("post-merger sweep and the bordered-Newton continuation; see")
    print("make_fig4_postmerger.py and bordered_newton_robustness.py.")
