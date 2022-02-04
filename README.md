# Stability-Boundary-Tracking
Matlab codes code to find the stability boundaries for cnoidal wave solutions to the Lugiato-Lefever Equation using dynamical methods.

This package can be used to reproduce Fig. 3 from the below [paper](https://github.com/Computational-Photonics-Laboratory/Stability-Boundary-Tracking/blob/main/PAJ295.pdf):
Z. Qi, S. Wang, J. Jaramillo-Villegas, M. Qi, A. M. Weiner, G. D'Aguanno, T. F. Carruthers, and C. R. Menyuk, "Dissipative cnoidal waves (Turing Rolls) and the soliton limit in microring resonators," Optica 6, 1220-1232 (2019). 

<b>Authors:</b> Logan Courtright, Zhen Qi, Thomas F. Carruthers, and Curtis R. Menyuk

To start: run "main.m" using MATLAB.

Initial Conditions: The program requires four things to be supplied as initial conditions. They are:
	"uout": a highly stable solution in the alpha-F parameter space
	"alpha_0": the location of "uout" in terms of alpha
	"F_0": the location of "uout" in terms of F
	"L": the normalized mode circumference
	"N": the spatial discretization of the solution. This value should not exceed 512.

"uout" is a (2*N,1) vector where the real part is from (1:N,1) and the imaginary part is from (N+1:2*N,1). The rest
of the initial conditions are constants.

The package uses the Lugiato-Lefever Equation (LLE) to obtain cnoidal wave solutions. For details on the LLE
normalizations that were used to obtain "alpha","F", and "L", see Z. Qi, et. al., Optica 6, 1220-1232 (2019).
The periodicity of the cnoidal wave is included in the initial solution "uout" and does not need to be specified.

The user can change three sets of additional "control" variables. They are:
	"alim_l": the lower alpha limit where the program will automatically stop
	"alim_h": the upper alpha limit where the program will automatically stop
	"da": the alpha step size
	"dF": the F step size
	"adir": the initial direction traveled in the alpha direction
	"Fdir": the initial direction traveled in the F direction

Suggested values for these variables are pre-set in the program, but should be changed if necessary.
For best operation of the program, "alim_h" should not be increased past its pre-set value. If this is necessary, 
"N" needs to be increased to at least 1024 to ensure correct operation, but this will drastically increase the 
runtime of the program.

All tests were done on a SAGER NP9873 laptop with 16 GB RAM using L = 50, N = 512, and the pre-set value of control 
variables. The time cost of the program varied from 2.5 - 4 hours depending on starting location and periodicity 
value. Higher periodicity values had lower time costs. 
