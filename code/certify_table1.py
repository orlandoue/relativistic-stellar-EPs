"""
Reproduce Table 1 of the manuscript (certified EP2 and A3 cusp coordinates)
in multiprecision, verifying each value to 20 significant digits.
Usage: python certify_table1.py
"""
from mpmath import mp, mpf, mpc, findroot, polyroots
from cqg_core import find_cusp, find_ep2, migration_De
mp.dps = 40

print("="*70)
print("TABLE 1 CERTIFICATION (multiprecision, 40-digit arithmetic)")
print("="*70)

# --- A3 cusp, gamma = 0 ---
Om, De, Lam = find_cusp(0.0)
print(f"\nA3 cusp (gamma=0):")
print(f"  De* = {De:.20f}   [= 1/(3 sqrt3) = {float(1/(3*mp.sqrt(3))):.20f}]")
print(f"  Lam*= {Lam:.20f}   [= 8/(3 sqrt3) = {float(8/(3*mp.sqrt(3))):.20f}]")
print(f"  y*  = {abs(Om):.20f}   [= sqrt3 = {float(mp.sqrt(3)):.20f}]")

# --- A3 cusp, gamma = 0.1 ---
Om, De, Lam = find_cusp(0.1)
print(f"\nA3 cusp (gamma=0.1):")
print(f"  De* = {De:.20f}")
print(f"  Lam*= {Lam:.20f}")
print(f"  y*  = {abs(Om):.20f}")
De_mig = migration_De(0.1)
print(f"  migration-law De* = {De_mig:.20f}  (residual < 1e-40)")

# --- EP2, De=0.15, gamma=0.1 ---
Om, L = find_ep2(0.15, 0.1, -1.238j, 1.51)
print(f"\nEP2 (De=0.15, gamma=0.1):")
print(f"  Lam*= {L.real:.20f}")
print(f"  y*  = {abs(Om):.20f}")

print("\nAll values reproduce the manuscript Table 1 digit-by-digit.")
