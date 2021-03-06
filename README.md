# QuantumSolver

This is a collection of matlab scripts that solve the Schrödinger equation in 1D and 2D in two different ways: by direct evaluation of the differential equation and by finding the eigenvectors of the Hamiltonian and propagating them analytically.

The direct evaluation ones, in quantum_noeigen.m and quantum_noeigen_2d.m, are based on the Crank-Nicolson method (https://en.wikipedia.org/wiki/Crank-Nicolson_method) which is guaranteed to be stable. The setup is a little more complex in two dimensions, as the matricies are N^2 x N^2 and they are band diagonal instead of tridiagonal.

The eigenvalue/eigenvector ones are more like how this would be done analytically. Once the eigenvectors and energy eigenvalues are found, the propagation is just a matter of propagating the eigenvectors analytically.


All code has been tested on Linux Mint 17.3 in Octave 4.0.2. To run it, simply call the file name from Octave.
