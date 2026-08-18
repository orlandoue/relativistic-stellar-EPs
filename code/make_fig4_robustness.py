"""
make_fig4_robustness.py -- figure 4 of the manuscript.

This file was listed in the README and referenced by make_all_figures.py and by
diagnostic_dense_continuation.py, but was missing from the archived repository,
so figure 4 could not be regenerated.  It is supplied here.

The figure is drawn from the archived continuation data, not recomputed, so
that what is plotted is exactly what is tabulated:

    ../data/robustness_continuation.json   (branch A, complex escape, R2, cusp)
    ../data/branches_Ap_B.json             (branches A' and B; produced by
                                            branches_Ap_B.py)

Run with --published to reproduce the figure exactly as it appears in the
manuscript (branch A plus the three birthplace markers, nothing else).
Run with no flag for the full version, which shows every object the caption
claims is followed: A, A', B, R2, the complex escape and the cusp trajectory.

Usage:  python make_fig4_robustness.py [--published]
"""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = "../figures/"
V0 = 25.25
TAU = 0.03
K = 5.0
ETA_C = 0.1941            # eta_c = 0.1941(1)


def cplx(s):
    return complex(s.replace('(', '').replace(')', '').replace(' ', ''))


def crossing(m, M):
    """Exact defective inter-sector crossing, Eq. (51)."""
    s = -(m + M + 1) / 2.0
    z = (1 + s * TAU) * (s**2 + (2 * m + 1) * s + m * (m + 1) + V0) / (-s)
    return s, z


def load(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def series(rows, i_s=1, i_z=2):
    e = [float(r[0]) for r in rows]
    s = [-cplx(r[i_s]).real for r in rows]
    z = [cplx(r[i_z]).real for r in rows]
    zi = [cplx(r[i_z]).imag for r in rows]
    return np.array(e), np.array(s), np.array(z), np.array(zi)


def main(published=False):
    d = load("../data/robustness_continuation.json")
    if d is None:
        sys.exit("missing ../data/robustness_continuation.json")
    extra = load("../data/branches_Ap_B.json")

    eA, sA, zA, _ = series(d["rama1"])
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 7.5), sharex=True)

    # --- branch A, the fundamental (sector-0) exceptional point
    a1.plot(eA, zA, 'o-', color='C0', ms=5, lw=1.6,
            label='Fundamental EP2 (branch A)')
    a1.plot([0], [zA[0]], 'o', mfc='none', mec='C0', ms=14,
            label=rf'analytic start $k\Lambda_*={zA[0]:.4f}$')
    a2.plot(eA, sA, 'o-', color='C0', ms=5, lw=1.6,
            label='Branch A degenerate eigenvalue')

    # --- exact birthplaces at eta = 0, Eq. (51)
    for (m, M) in [(2, 10), (0, 12), (0, 16)]:
        s0, z0 = crossing(m, M)
        a1.plot(0, z0, '*', color='k', ms=15)
        a1.annotate(f'({m},{M})', (0.002, z0), fontsize=8, va='center')

    if not published:
        # --- branch A', the annihilation partner, born on the (2,10) crossing
        if extra and "branch_Aprime_from_(2,10)" in extra:
            e, s, z, _ = series(extra["branch_Aprime_from_(2,10)"])
            a1.plot(e, z, 's-', color='C3', ms=4, lw=1.3,
                    label="Partner EP2 (branch A'), from (2,10)")
            a2.plot(e, s, 's-', color='C3', ms=4, lw=1.3)
        # --- branch B, born on the (0,12) crossing
        if extra and "branch_B_from_(0,12)" in extra:
            e, s, z, _ = series(extra["branch_B_from_(0,12)"])
            a1.plot(e, z, '^-', color='C2', ms=4, lw=1.3,
                    label='Branch B, from (0,12)')
            a2.plot(e, s, '^-', color='C2', ms=4, lw=1.3)
        # --- R2, the further deformed sectoral branch
        eR, sR, zR, _ = series(d["rama2"])
        o = np.argsort(eR)
        a1.plot(eR[o], zR[o], 'v-', color='C4', ms=4, lw=1.3,
                label='Branch R2 (parent not established)')
        a2.plot(eR[o], sR[o], 'v-', color='C4', ms=4, lw=1.3)
        # --- complex escape of the merged pair, eta > eta_c: plot Re(zeta0*)
        eC, sC, zC, zCi = series(d["rama_compleja"])
        a1.plot(eC, zC, 'd--', color='C1', ms=5, lw=1.3,
                label=r'merged pair, $\mathrm{Re}\,\zeta_0^*$ ($\eta>\eta_c$)')
        a2.plot(eC, sC, 'd--', color='C1', ms=5, lw=1.3)
        # --- cusp trajectory: k*zeta0* of the EP3 branch
        eK = [float(r[0]) for r in d["cusp"]]
        zK = [cplx(r[2]).real for r in d["cusp"]]
        sK = [-cplx(r[3]).real for r in d["cusp"]]
        a1.plot(eK, zK, 'x-', color='0.4', ms=6, lw=1.2,
                label=r'$A_3$ cusp trajectory ($\tau$ free)')
        a2.plot(eK, sK, 'x-', color='0.4', ms=6, lw=1.2)

    for ax in (a1, a2):
        ax.axvline(ETA_C, color='red', ls=':', lw=1.3)
        ax.grid(alpha=0.3)
    a1.annotate(rf'$\eta_c={ETA_C:.4f}$', (ETA_C + 0.006, a1.get_ylim()[0] + 0.15),
                color='red', fontsize=9, ha='left')
    a1.set_ylabel(r'Critical coupling $\zeta_0^*$')
    a2.set_ylabel(r'$-s^* = k\,y^*$')
    a2.set_xlabel(r'Deformation parameter $\eta$')
    a1.legend(fontsize=8, loc='best')
    if published:
        a2.legend(fontsize=8, loc='upper left')
    a1.set_title('Robustness of the fundamental exceptional point\n'
                 'under transport-profile deformation')
    fig.tight_layout()
    name = 'fig5_robustness_branches.pdf' if published \
        else 'fig5_robustness_branches_full.pdf'
    fig.savefig(OUT + name, dpi=150, bbox_inches='tight')
    print(f"wrote {OUT}{name}")
    if not published and extra is None:
        print("  (branches A' and B not drawn: run branches_Ap_B.py first)")


if __name__ == "__main__":
    main(published="--published" in sys.argv)
