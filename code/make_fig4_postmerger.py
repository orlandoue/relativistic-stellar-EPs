"""
Figure 4: post-merger remnant QNM evolution (De=0.11, gamma_GW=0.03).
The sweep in Lambda passes exactly through the two real-axis EP2 coalescences.
Cproj = 5.32e-32 is the fiducial spatial-projection coefficient (see sec 9.5).
Usage: python make_fig4_postmerger.py
"""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from cqg_core import roots_np, find_ep2

De, g, Cproj, f0 = 0.11, 0.03, 5.32e-32, 3.5e3
# exact EP2 of this (De,gamma):
L1 = find_ep2(De, g, -1.15j, 1.71)[1].real   # ~1.712
L2 = find_ep2(De, g, -4.33j, 2.36)[1].real   # ~2.357
z1, z2 = L1/Cproj, L2/Cproj
print(f"Exact EP2: Lambda* = {L1:.6f} (zeta={z1:.3e}), {L2:.6f} (zeta={z2:.3e})")

base = np.linspace(1.5e31, 5.5e31, 500)
dense = np.concatenate([np.linspace(z1*0.98,z1*1.02,400), np.linspace(z2*0.98,z2*1.02,400)])
zetas = np.sort(np.concatenate([base, dense])); Lams = Cproj*zetas
F = np.zeros((len(Lams),3), dtype=complex)
for i,L in enumerate(Lams):
    r = roots_np(De,g,L); F[i] = r[np.argsort(r.imag)]
Ref, Imf = F.real*f0, F.imag*f0
gap = np.minimum(np.abs(F[:,0]-F[:,1]), np.abs(F[:,1]-F[:,2]))*f0

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5))
for j,c in enumerate(['C0','C1','C2']):
    ax1.plot(zetas/1e31, Ref[:,j],'.',color=c,ms=2,label=f'Mode {j+1}')
for ax in (ax1,ax2):
    ax.axvline(z1/1e31,color='red',ls=':',lw=1.3,alpha=0.7)
    ax.axvline(z2/1e31,color='purple',ls=':',lw=1.3,alpha=0.7)
ax1.set_xlabel('Bulk Viscosity ($10^{31}$ g cm$^{-1}$ s$^{-1}$)'); ax1.set_ylabel('Re($f$) (Hz)')
ax1.set_title('Real frequency across the two coalescences'); ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
ax1.twiny().set_xlim(Cproj*1.5e31,Cproj*5.5e31)
ax2.semilogy(zetas/1e31, gap+1,'-',color='darkgreen',lw=1.5,label='min eigenvalue gap')
ax2.annotate(r'EP2 #1 $\Lambda$=1.712',(z1/1e31,5e3),fontsize=8,color='red',ha='center')
ax2.annotate(r'EP2 #2 $\Lambda$=2.357',(z2/1e31,5e3),fontsize=8,color='purple',ha='center')
ax2.set_xlabel('Bulk Viscosity ($10^{31}$ g cm$^{-1}$ s$^{-1}$)'); ax2.set_ylabel(r'Eigenvalue gap (Hz)')
ax2.set_title('Gap collapses at both exact EP2 coalescences'); ax2.legend(fontsize=8); ax2.grid(alpha=0.3,which='both')
ax2.twiny().set_xlim(Cproj*1.5e31,Cproj*5.5e31)
fig.suptitle(r'Post-Merger Remnant QNM Evolution: Exact EP2 Coalescences ($\tau$=0.005 ms, $\gamma_{GW}$=0.03)')
fig.tight_layout(); fig.savefig('../figures/fig4_physical_NS_application.pdf',dpi=150,bbox_inches='tight')
print("fig4 saved")
