# Grand CFD Portfolio  
### Numerical Methods • Stability Analysis • Convergence Verification • HPC Foundations  

**Author:** Max Minh Le  

---

## Overview

This repository is the first core stage of a structured Computational Fluid Dynamics portfolio focused on research-grade solver development.

The objective is not merely to “solve PDEs,” but to:

- Analyze numerical stability
- Compare time integration schemes
- Verify convergence behavior
- Demonstrate controlled instability
- Implement implicit solvers from first principles
- Build modular, extensible solver architecture
- Prepare for elliptic, hyperbolic, and Navier–Stokes systems

The current stage focuses on **Parabolic PDEs (Heat Equation)** in 1D and 2D.

---

# Pillar 1 — Parabolic Problems (Heat Equation)

## Governing Equation

\[
\frac{\partial u}{\partial t} = \alpha \nabla^2 u
\]

1D:
\[
\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}
\]

2D:
\[
\frac{\partial u}{\partial t} = \alpha \left(
\frac{\partial^2 u}{\partial x^2} +
\frac{\partial^2 u}{\partial y^2}
\right)
\]

---

# 1D Heat Equation

## Implemented Schemes

- Explicit FTCS
- Fully Implicit FTCS
- Dufort–Frankel
- Crank–Nicolson

All schemes are compared against the analytical solution.

---

## What Is Demonstrated

### 1. Stability Behavior

- Explicit FTCS stable under Fourier constraint
- Explicit FTCS catastrophic growth when stability violated
- Dufort–Frankel oscillatory error structure
- Implicit FTCS unconditionally stable but diffusive
- Crank–Nicolson smooth and higher-order behavior

Instability is shown intentionally rather than avoided.

---

### 2. Error Distribution Analysis

For fixed time step:

- Spatial error profiles plotted against exact solution
- Scheme-dependent oscillations identified
- Error sign behavior examined
- Comparison of diffusion vs dispersion characteristics

---

### 3. Time Convergence Study

- Time step successively halved
- L2 and L∞ norms computed
- Log–log convergence plots generated
- Observed order of accuracy verified
- Explicit instability clearly visible at unsafe time steps

This confirms theoretical time accuracy.

---

# 2D Heat Equation

## Implemented Schemes

- Explicit FTCS (with stability monitoring)
- Implicit ADI (Alternating Direction Implicit)

---

## Stability Verification

The 2D explicit scheme enforces:

\[
f_x + f_y \le 0.5
\]

Two cases are shown:

- Stable explicit solution (safe time step)
- Deliberate instability case (checkerboard blow-up pattern)

Instability is documented and analyzed, not hidden.

---

## Implicit ADI Implementation

- ADI splitting method
- Independent tridiagonal Thomas solver
- Fully implicit 2D time advancement
- Stable large time step capability

The ADI implementation serves as reference solution for convergence studies.

---

## 2D Time Convergence Study

- Time step refinement sequence
- L2 and L∞ error norms
- Log–log convergence plots
- Demonstration of expected accuracy behavior

This verifies correctness of the ADI implementation.

---

# Code Architecture

## Languages

- C++ (single-file and modular implementations complete)
- Fortran (in progress)
- Python (in progress)

Each language mirrors the same numerical structure.

---

## Modular Design

C++ implementation structured by responsibility:

- main driver
- simulation driver
- convergence driver
- explicit solver
- implicit ADI solver
- Thomas tridiagonal solver
- error norm computation
- IO utilities
- plot generation
- directory management

No external numerical libraries are used.

All linear algebra routines are implemented from first principles.

---

# Numerical Foundations Established

This stage demonstrates:

- Stability analysis and enforcement
- Controlled instability experiments
- Time integration comparison
- Error norm computation (L2, L∞)
- Log–log convergence verification
- Extension from 1D to 2D
- Independent tridiagonal solver development
- Structured separation between simulation and verification workflows

These foundations prepare for:

- Elliptic Poisson solvers
- Hyperbolic wave systems
- Incompressible Navier–Stokes
- Compressible shock-capturing
- Parallel MPI/OpenMP extensions

---

# Development Philosophy

1. Stability and convergence are demonstrated, not assumed.
2. Instability is shown explicitly when theoretical limits are violated.
3. All solvers are implemented from first principles.
4. Numerical behavior is analyzed quantitatively.
5. Architecture is designed for extensibility toward full CFD systems.

---

# Current Status

- 1D heat solver complete (C++)
- 2D heat solver complete (C++)
- Convergence verification complete
- Fortran and Python ports in progress
- Elliptic solvers planned next

---

This repository evolves incrementally according to a structured CFD roadmap toward full research-grade solver development.
