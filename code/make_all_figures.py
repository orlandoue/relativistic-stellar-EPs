"""
Regenerate figures 1, 2 and 3 of the manuscript from the certified core (cqg_core.py).
Every figure derives from the single validated Open Cusp Polynomial convention.
Usage: python make_all_figures.py   (outputs to ../figures/)

FIX (v2), relative to the previously archived version
-----------------------------------------------------
The EP2 locus of the closed pencil is parametrised by y = |Omega| through
    De(y) = (y^2-1)/(2y^3),      Lambda(y) = (y^2+1)^2/(2y^3),
and De(y) has a MAXIMUM at y = sqrt(3), which is the A3 cusp.  The locus is
therefore a genuine cusp with TWO branches over every De < De*, one with
y < sqrt(3) and one with y > sqrt(3).

The previous version sorted the sampled points by De (`np.argsort(pts[:,0])`)
before plotting.  That interleaves the two branches and makes matplotlib
zig-zag between them, producing the spurious sawtooth "fan" visible in the
published figures 1 and 2 between De ~ 0.105 and De ~ 0.19.  The fix is to
order by the natural parameter y and draw the two branches as two curves.
Nothing about the underlying numbers changes; only the polyline ordering.
"""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpmath import mp, mpf, mpc, findroot
from cqg_core import find_cusp, find_ep2, roots_np, roots_mp
mp.dps = 40
OUT = "../figures/"

SQRT3 = float(np.sqrt(3.0))


def ep2_closed(y):
    """EP2 of the closed (gamma=0) pencil at Omega = -i y, from the exact
    parametrisation (30).  Closed form: no root-finding needed, hence no
    possibility of the Newton iteration jumping branch."""
    y = float(y)
    return (y**2 - 1.0) / (2.0 * y**3), (y**2 + 1.0)**2 / (2.0 * y**3)


def ep2_branches(y_lo=1.02, y_hi=4.5, n=600):
    """Return (lower, upper) branches of the closed EP2 locus, each ordered by
    the natural parameter y so that the polyline never jumps between them.
    Lower branch: 1 < y < sqrt(3).  Upper branch: y > sqrt(3).  They meet, with
    a common tangent, at the A3 cusp y = sqrt(3)."""
    y_low = np.linspace(y_lo, SQRT3, n // 2)
    y_up = np.linspace(SQRT3, y_hi, n // 2)
    low = np.array([ep2_closed(y) for y in y_low])
    up = np.array([ep2_closed(y) for y in y_up])
    return low, up


# ---------- FIG 1: global EP2 curve ----------
def fig1():
    low, up = ep2_branches(1.02, 4.5, 600)
    Om, Dc, Lc = find_cusp(0.0)
    Om1, Dc1, Lc1 = find_cusp(0.1)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot(low[:, 0], low[:, 1], '-', color='C0', lw=2,
            label=r'EP2 curve ($\Delta=0$), $1<y<\sqrt{3}$')
    ax.plot(up[:, 0], up[:, 1], '--', color='C0', lw=2,
            label=r'EP2 curve ($\Delta=0$), $y>\sqrt{3}$')
    ax.plot(0, 2, 's', color='C2', ms=10, label='Navier-Stokes EP2 $(0,2)$')
    ax.plot(Dc, Lc, '*', color='red', ms=18, label=r'$A_3$ cusp $\gamma=0$')
    ax.plot(Dc1, Lc1, 'D', color='darkred', ms=8, label=r'migrated cusp $\gamma=0.1$')
    ax.set_xlabel('Deborah number $De$')
    ax.set_ylabel(r'Viscous coupling $\Lambda$')
    ax.set_title('Global exceptional-point (EP2) curve')
    ax.legend(fontsize=8.5); ax.grid(alpha=0.3)
    ax.set_xlim(-0.02, 0.65); ax.set_ylim(0, 3.2)
    fig.tight_layout()
    fig.savefig(OUT + 'fig1_ep_curve.pdf', dpi=150, bbox_inches='tight')
    plt.close()


# ---------- FIG 2: cusp zoom with exact discriminant regions ----------
def disc_closed(De, Lam):
    """Discriminant of the closed cubic in Omega, real-valued.  Expanded form,
    verified symbolically against 18abcd-4b^3d+b^2c^2-4ac^3-27a^2d^2 with
    a=-i De, b=1, c=i(De+Lambda), d=-1."""
    return (4*De**4 + 12*De**3*Lam + 12*De**2*Lam**2 + 8*De**2
            + 4*De*Lam**3 - 20*De*Lam - Lam**2 + 4)


def fig2():
    Om, Dc, Lc = find_cusp(0.0)
    low, up = ep2_branches(1.30, 2.60, 800)
    DD = np.linspace(0.12, 0.28, 400)
    LL = np.linspace(1.2, 1.9, 400)
    DEg, LAg = np.meshgrid(DD, LL)
    # Exact sign of the discriminant instead of a threshold on Im(root).
    # Delta < 0  -> three distinct real y (three overdamped modes)
    # Delta > 0  -> one real y + a propagating conjugate pair
    S = np.sign(disc_closed(DEg, LAg))
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.contourf(DEg, LAg, S, levels=[-1.5, 0, 1.5],
                colors=['#e8f0ff', '#ffe8e8'], alpha=0.6)
    ax.plot(low[:, 0], low[:, 1], '-', color='C0', lw=2,
            label=r'EP2 locus ($\Delta=0$), $y<\sqrt{3}$')
    ax.plot(up[:, 0], up[:, 1], '--', color='C0', lw=2,
            label=r'EP2 locus ($\Delta=0$), $y>\sqrt{3}$')
    ax.plot(Dc, Lc, '*', color='red', ms=20, label=r'$A_3$ cusp')
    # The Delta<0 region (three distinct real y, i.e. three overdamped modes)
    # is the NARROW wedge between the two branches; everywhere else Delta>0
    # (one overdamped mode plus a damped oscillatory doublet).
    ax.annotate(r'$\Delta<0$: three overdamped modes',
                xy=(0.160, 1.691), xytext=(0.128, 1.78), fontsize=8, color='0.2',
                arrowprops=dict(arrowstyle='->', color='0.4', lw=0.8))
    ax.text(0.124, 1.27, r'$\Delta>0$: one overdamped mode $+$ oscillatory doublet',
            fontsize=8, color='0.25')
    ax.set_xlabel('Deborah number $De$')
    ax.set_ylabel(r'Viscous coupling $\Lambda$')
    ax.set_title(r'Local structure near the $A_3$ cusp')
    ax.legend(fontsize=9, loc='lower right')
    ax.set_xlim(0.12, 0.28); ax.set_ylim(1.2, 1.9); ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT + 'fig2_cusp_zoom.pdf', dpi=150, bbox_inches='tight')
    plt.close()


# ---------- FIG 3: certified splitting over the certified range ----------
def fig_splitting(fname, title):
    De_ep = 0.15; g = 0.1
    Om_ep, L_ep = find_ep2(De_ep, g, -1.238j, 1.51); Lstar = L_ep.real
    Om_c, Dcusp, Lcusp = find_cusp(0.1)

    def se(eps):
        r = [complex(x) for x in roots_mp(De_ep, g, Lstar * (1 + eps))]
        return sorted(abs(r[i] - r[j]) for i in range(3) for j in range(i + 1, 3))[0]

    def sc(eps):
        r = [complex(x) for x in roots_mp(Dcusp, g, Lcusp * (1 + eps))]
        return sorted(abs(r[i] - r[j]) for i in range(3) for j in range(i + 1, 3))[-1]

    epss = np.logspace(-10, -2, 25)
    yE = np.array([float(se(e)) for e in epss])
    yC = np.array([float(sc(e)) for e in epss])
    mask = epss <= 1e-4
    expE = np.polyfit(np.log(epss[mask]), np.log(yE[mask]), 1)[0]
    expC = np.polyfit(np.log(epss[mask]), np.log(yC[mask]), 1)[0]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.loglog(epss, yE, 'o', color='C0', ms=6,
              label=f'EP2: min gap (slope {expE:.4f}, $\\to 1/2$)')
    ax.loglog(epss, yC, 's', mfc='none', mec='red', ms=7,
              label=f'$A_3$ cusp: max gap (slope {expC:.4f}, $\\to 1/3$)')
    ax.loglog(epss, yE[0] * (epss / epss[0])**0.5, '--', color='C0', alpha=0.5,
              lw=1, label=r'$\epsilon^{1/2}$')
    ax.loglog(epss, yC[0] * (epss / epss[0])**(1 / 3), '--', color='red',
              alpha=0.5, lw=1, label=r'$\epsilon^{1/3}$')
    ax.set_xlabel(r'relative detuning $|\epsilon|$')
    ax.set_ylabel(r'eigenvalue splitting $|\delta\Omega|$')
    ax.set_title(title); ax.legend(fontsize=8.5, loc='lower right')
    ax.grid(alpha=0.3, which='both')
    fig.tight_layout(); fig.savefig(OUT + fname, dpi=150, bbox_inches='tight')
    plt.close()
    # local slopes, decade by decade -- reported so the text can quote the
    # range over which 0.50000 actually holds
    slopes = [(float(epss[i]),
               float((np.log(yE[i]) - np.log(yE[i - 1])) /
                     (np.log(epss[i]) - np.log(epss[i - 1]))))
              for i in range(1, len(epss))]
    return expE, expC, slopes


if __name__ == "__main__":
    print("Generating figures from certified core...")
    fig1(); print("  fig1_ep_curve.pdf   (two branches, ordered by y)")
    fig2(); print("  fig2_cusp_zoom.pdf  (exact discriminant sign)")
    eE, eC, sl = fig_splitting('fig3_splitting.pdf',
                               'Non-analytic scaling of the eigenvalue splitting')
    print(f"  fig3_splitting.pdf  (asymptotic exponents: EP2={eE:.5f}, cusp={eC:.5f})")
    fig_splitting('fig6_certified_scaling.pdf',
                  'Certified non-analytic splitting at the degeneracies')
    print("  fig6_certified_scaling.pdf  (= manuscript figure 3)")
    print("\n  local EP2 slope, top and bottom of the sampled range:")
    print(f"    eps={sl[0][0]:.1e}: {sl[0][1]:.5f}   eps={sl[-1][0]:.1e}: {sl[-1][1]:.5f}")
    print("    -> 0.50000 holds for eps <= 1e-5; quote the range, not 'all decades'.")
    print("\nFigure 4: make_fig4_robustness.py.  Figure 5: make_fig5_postmerger.py.")
