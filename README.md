# MATLAB optimisation program to search for time-periodic solutions of the Ostrovsky equation

### - Attempts to minimise an objective that is zero if and only if the solution is time-periodic

### - Uses a quasi-Newton BFGS and adjoint integration to compute a varitional derivative wrt an initial wavefunction

OstrovskyCode and KdVCode contain all of the functions needed to optimise for time-periodic solutions.
There are no dependancies, all functions are included within these two folders. They also contain plots and results scripts 
that can be modified to study particular problems/cases.

Use the function help in each .m file for a more in-depth explanation of the implementation.

KdV code has two BFGS functions that seek solutions on a periodic domain. One implementing a representation of initial 
wave function on evenly spaced Fourier nodes and the other using Fourier mode coefficients. 

Ostrovsky code has several BFGS functions for different representations including Hermite nodes, and allows an optimatation
of wave-packets by including a dampening region in the integration domain. 
