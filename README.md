--------------------------------------------------------------------------------
# **The AMUN Code**
## Copyright (C) 2008-2018 Grzegorz Kowal ##
--------------------------------------------------------------------------------

AMUN is a parallel code to perform numerical simulations in fluid approximation
on uniform or non-uniform (adaptive) meshes. The goal in developing this code is
to create a solid framework for simulations with support for number of numerical
methods which can be selected in an easy way through the parameter file. The
following features are already implemented:

* hydrodynamic and magnetohydrodynamic set of equations (HD and MHD),
* both classical and special relativity cases for the above equations,
* Cartesian coordinate system,
* uniform and adaptive mesh generation and update,
* 2nd to 4th order time integration using Strong Stability Preserving
  Runge-Kutta methods,
* 2nd order TVD interpolation with number of limiters and higher order
  reconstructions,
* Riemann solvers of Roe- and HLL-types (HLL, HLLC, and HLLD),
* periodic and open boundary conditions,
* viscous and resistive source terms,
* data stored in the HDF5 format,
* MPI parallelization,
* completely written in Fortran 2003.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.


Developers
==========

 - Grzegorz Kowal <grzegorz@amuncode.org>


Requirements
============

* Fortran 2003 compiler (tested compilers include
  [GNU Fortran](http://gcc.gnu.org/fortran/) version 4.5 or newer,
  [Intel Fortran](https://software.intel.com/en-us/fortran-compilers) compiler
  version 9.0 or newer)
* [HDF5 libraries](http://www.hdfgroup.org/HDF5/) version 1.8 or newer.
* [OpenMPI](https://www.open-mpi.org/) version 1.8 or newer for parallel runs.


Environment Variables
=====================

If the HDF5 libraries are not installed in the default location, i.e. in the
system directory **/usr**, make sure that the environment variable _HDF5DIR_ is
set in your **~/.bashrc** (or **~/.cshrc**) and pointing to the location where
the HDF5 libraries have been installed.


Compilation
===========
1. Clone the AMUN source code: `git clone https://bitbucket.org/amunteam/amun-code.git`,
   or unpack the archive downloaded from page
   [Downloads](https://bitbucket.org/amunteam/amun-code/downloads/).
2. Go to directory **build/hosts/** and copy file **default** to a new file named
   exactly as your host name (name returned by command `hostname`).
3. Customize your compiler and compilation options in your new host file.
4. Go up to directory **build/** and copy file **make.default** to **make.config**.
5. Customize compilation time options in **make.config**.
6. Compile sources by typing `make` in directory **build/**. The executable file
   **amun.x** should be created there.


Usage
=====

In order to run some test problems you can simply copy corresponding parameter
from directory **problems/** to the location when you wish to run your test.
Copy the executable file **amun.x** compiled earlier to the same directory. If
you provide option _-i <parameter_file>_, the code will know that the parameters
have to be read from file _<parameter_file>_. If you don't provide this option,
the code will assume that the parameters are stored in file **params.in** in the
same director.

In order to run serial version, type in your terminal:  `amun.x -i params.in`.

In order to run the parallel version (after compiling the code with MPI
version), type in your terminal: `mpirun -n N ./amun.x -i params.in`, where N is
the number of processors.
