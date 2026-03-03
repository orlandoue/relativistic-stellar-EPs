# Exceptional points and cusp singularities in relativistic stellar oscillations

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

This repository contains the symbolic derivations, numerical routines, and plotting scripts required to reproduce the geometric and phenomenological results presented in the paper:

> **Exceptional points and cusp singularities in the effective dynamics of relativistic stellar oscillations with causal bulk relaxation** > *Orlando Urbina-Gonzalez* > Pontificia Universidad Católica de Chile (UC)  
> [Link to arXiv / CQG Paper] (Coming soon)

## Abstract
Causal relativistic hydrodynamics (Israel-Stewart theory) elevates the standard quadratic self-adjoint stellar oscillation problem to a cubic non-Hermitian spectral problem. This work demonstrates that this effective macroscopic operator possesses a structured non-Hermitian spectral geometry governed by its discriminant. The spectrum forms a three-sheeted Riemann surface containing a continuous curve of second-order exceptional points (EP2), which terminates at a higher-order $A_2$ (cusp) singularity precisely at the Navier-Stokes boundary ($De \to 0$). We map this topology to phenomenological parameters, demonstrating that the structural transitions are dynamically accessible in post-merger hypermassive neutron stars.

## Repository Structure

All figures generated for the manuscript are fully reproducible using the Python scripts located in the `src/` directory.

* `src/fig1_ep_curve.py`: Solves the cubic discriminant locus to extract the global exceptional-point (EP2) curve in the $(De, \Lambda)$ parameter plane.
* `src/fig2_cusp_zoom.py`: Computes the local geometric structure of the $A_2$ cusp singularity at the Navier-Stokes limit ($De = 0$).
* `src/fig3_splitting.py`: Evaluates the fractional-power branching of the eigenvalues, confirming the non-analytic square-root splitting ($\sim |\epsilon|^{1/2}$) near generic EP2s and the cubic-root scaling ($\sim |\epsilon|^{1/3}$) at the cusp.
* `src/fig4_physical_NS_application.py`: Maps the dimensionless normal form to physical transport scales representing a hypermassive post-merger remnant ($M = 2.7 M_\odot$, $f_0 = 3.5$ kHz, $\tau = 0.005$ ms). It performs a continuous root-tracking algorithm to phenomenologically demonstrate the double-crossing of the EP2 locus driven by Urca-process bulk viscosities ($\zeta \sim 10^{31}$ cgs).

## Installation and Usage

To run the simulations locally, clone this repository and install the standard scientific dependencies:

```bash
git clone [https://github.com/YOUR_USERNAME/relativistic-stellar-EPs.git](https://github.com/YOUR_USERNAME/relativistic-stellar-EPs.git)
cd relativistic-stellar-EPs
pip install -r requirements.txt
