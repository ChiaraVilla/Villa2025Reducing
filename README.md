# Villa2025Reducing
Code to solve a phenotype-structured PDE and the system of ODEs approximating the phenotypic distribution moment dynamics

**Public gitlab repository Villa2025Reducing** <br />
This repository provides Matlab files to perform simulations as described in: <br />
_Reducing phenotype-structured PDE models of cancer evolution to systems of ODEs: a generalised moment dynamics approach_ by _C. Villa, P. Maini, A. Browning, A. Jenner, S. Hamis, T. Cassidy, preprint arXiv:2406.01505 (2025); <br />

For details we refer the interested reader to this preprint (currently under review). <br />
If you use this software in your work then please cite the above named paper.

**How to use** <br />
- Change the definitions of f and V in the V_DEF(x,par) and f_DEF(x,par) functions at the end of the script
- Change the value of "N" moments you wish to track in the approximation (N=2,3 allowed) and of "M" at which to truncate the Taylor series approximation of f and V (M=1,2,3 allowed), in the "Problem set up" section of the script
- Choose between Gaussian ("gaussian") or truncation ("truncate") closure for the ODE system, editing this in the "Problem set up" section of the script
- Edit the date at the top of the script (this is used in the name of the file to save the results)
- Run the scrip

**Notes to users** <br />
_Note 1:_ If your functions are time dependent you many need to edit these to allow for the extra input variable, along with functions dnf_dxn(M,par) and dnV_dxn(M,par), and fn =@(n,mi) and Vn =@(n,mi) in the code. <br />
_Note 2:_ The code uses explicit schemes to solve the PDE and ODE system. This my not be the best choice for custom f and V and for all combinations of N and M, and the scheme may be unstable in certain cases: try reducing the timestep dt (in the "Problem set up" section of the script) or consider using an implicit scheme.

**Copyright notice** <br />
Villa2025Reducing: simulate phenotype-structured PDE and ODE systems approximating the phenotypic distribution moment dynamics <br />
Copyright (C) 2025 C. Villa

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

