## A brief description of the files

### Defination of the model

- A_BdG.m Define the fermionic Hamiltonian in terms of Majorana fermions
- H_BdG.m Define the BdG Hamiltonian in terms of complex fermions
- parameterize.m Define the parameterization of the model 
(*note there are several different vison paths depends on which type of processes we are interested in*)

### Calculating the spectra of the BdG Hamiltonian and BCS states' overlap

- spectra.m Diagonalizing the BaG Hamiltonian, calculating the ground state parity
- overlap.m The matrix element of the (dual of) spin-Z operator between two fermionic BCS states
- overlap_X.m The matrix element of the spin-X operator between two fermionic BCS states
- overlap_Y.m The matrix element of the spin-Y operator between two fermionic BCS states
- hole_overlap* Same matrix elements but with a hole-like BCS ansatz (*suitable for AFM coupled Kitaev model*)
- berry_phase.m Calculate the resulting Berry phase from the calculated vison hopping matrix elements
(*generated by the overlap/hole_overlap files above*)

### Some "auxiliary" packages

- fdet.m Newly defined function for calculating the determinant of a matrix, this is to avoid numerical "0" error for large size matrix
- findex.m Transfer a Wannier state's real space coordinate into its index inside the basis
- fhop.m Give the correct Wannier states' index when defineing the BdG Hamiltonian
- pfaffian_householder.m This is from a package by [M. Wimmer](https://arxiv.org/abs/1102.3440) to calculate the Pfaffian of a skew-symmetric a matrix
