(* ::Package:: *)

(* =============================================================================
   projection_coefficient.wl
   -----------------------------------------------------------------------------
   Derivation of the physical projection coefficient C_proj that maps the
   physical bulk viscosity zeta (g cm^-1 s^-1) to the dimensionless viscous
   coupling Lambda of the open cubic normal form of CQG-115250:

       Lambda = zeta / (rho * L^2 * omega0)  ==>  C_proj = 1/(rho * L^2 * omega0)

   The coefficient is fixed entirely by the fiducial post-merger star already
   declared in the manuscript (M = 2.7 M_sun, R = 15 km thermally-inflated
   radius, fundamental compressive frequency f0 = 3.5 kHz; cf. Radice 2022).
   No equation-of-state input beyond (M, R) is required:

     * rho   : mean density, rho_bar = 3 M / (4 pi R^3)
     * L      : macroscopic mode length, L = R  (the fundamental compressive
                mode has no radial nodes, k_eff ~ 1/R)
     * omega0 : angular frequency of the fundamental mode, omega0 = 2 pi f0,
                with f0 the value quoted in the manuscript.

   This script is self-contained and reproducible; it prints every intermediate
   quantity and then maps the two exact sector-0 EP2 couplings to their
   physical viscosities.
   ============================================================================= *)

(* ---- physical constants (CGS) ---- *)
GG   = 6.674*10^-8;     (* cm^3 g^-1 s^-2 *)
cc   = 3.0*10^10;       (* cm/s *)
Msun = 1.989*10^33;     (* g *)

(* ---- fiducial star (as declared in the manuscript) ---- *)
Mstar = 2.7*Msun;      (* g   *)
Rkm   = 15;            (* km  -- thermally inflated post-merger radius *)
Rstar = Rkm*10^5;      (* cm  *)
f0kHz = 3.5;           (* kHz -- fundamental compressive-mode frequency *)

(* ---- (1) mean density ---- *)
rhoBar = 3*Mstar/(4*Pi*Rstar^3);
Print["Mean density   rho_bar = 3M/(4 pi R^3) = ",
      ScientificForm[N[rhoBar], 5], " g/cm^3"];

(* ---- (2) macroscopic mode length ---- *)
Lscale = Rstar;   (* fundamental mode: k_eff ~ 1/R  =>  L = R *)
Print["Mode length    L = R = ", ScientificForm[N[Lscale], 5], " cm (",
      N[Rkm], " km)"];

(* ---- (3) fundamental angular frequency ---- *)
omega0 = 2*Pi*f0kHz*10^3;   (* rad/s *)
Print["Angular freq   omega0 = 2 pi f0 = ", ScientificForm[N[omega0], 5],
      " rad/s  (f0 = ", N[f0kHz], " kHz)"];

(* ---- projection coefficient ---- *)
Cproj = 1/(rhoBar*Lscale^2*omega0);
Print["-----------------------------------------------------------"];
Print["PROJECTION COEFFICIENT  C_proj = 1/(rho L^2 omega0) = ",
      ScientificForm[N[Cproj], 5], " g^-1 cm s"];
Print["Compactness    GM/(R c^2) = ", N[GG*Mstar/(Rstar*cc^2), 4]];
Print["-----------------------------------------------------------"];

(* ---- map the two exact sector-0 EP2 couplings to physical viscosity ----
   Lambda* are the exact EP2 coalescences of the open cubic at De = 0.11,
   gamma = 0.03 (see cqg_core.py / certify_table1.py).                    *)
Lam1 = 1.712359;
Lam2 = 2.357256;
Print["Sector-0 EP2 #1:  Lambda* = ", Lam1,
      "  ->  zeta = ", ScientificForm[N[Lam1/Cproj], 4], " g cm^-1 s^-1"];
Print["Sector-0 EP2 #2:  Lambda* = ", Lam2,
      "  ->  zeta = ", ScientificForm[N[Lam2/Cproj], 4], " g cm^-1 s^-1"];

(* ---- emit the (zeta, Lambda) grid used for the fig. 6 sweep ---- *)
zetaGrid = Table[z, {z, 1.5*10^31, 5.5*10^31, (5.5-1.5)*10^31/400}];
Export["projection_grid.csv",
       Prepend[Transpose[{zetaGrid, Cproj*zetaGrid}], {"zeta_cgs","Lambda"}]];
Print["Wrote projection_grid.csv (", Length[zetaGrid], " rows)."];
