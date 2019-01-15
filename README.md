## Numerical Optimisation Final Project

The demonstration file is `solution.m` and takes roughly 40 seconds to run.

Guide to support functions:
* `phi_fourier_half_plane.m` - returns real and imaginary components of 2D Fast Fourier Transform on upper half plane of transform domain.
* `phi_t_fourier_half_plane.m` - returns real and imaginary components of 2D Fast Fourier Transform on bottom half plane of transform domain.
* `LineMask.m` - generates k-space frequencies used in the sampling pattern (by Justin Romberg).
* `logbar_L1.M` - log-barrier algorithm for L1 minimisation with quadratic constraints.
* `logbar_TV.m` - log-barrier algorithm for Total Variation minimisation with quadratic constraints.
* `primal_dual_L1.m` - primal-dual algorithm for L1 minimisation with equality constraints.
