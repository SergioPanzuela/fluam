####Fluam a fluctuating hydrodynamic code
###About this branch
The **Colors** branch allows you to stablish a different interaction for every pair based on an assigned particle type. Currently included only for the Stokes Limit and QuasiNeutrallyBuoyant schemes

You can turn this option on in data.main using the line "colors       1". This way the program will look for the colors in the particle coord file.


Therefore, with colors on, the coords input needs an additional column containing the particle types, starting from 0.

    N
    x1 y1 z1 type1
    x2 y2 z2 type2
    . 
    .
    .

Additionally you have to specify how many different types there are and two matrices with the parameters for both short and long range interactions, this needs to be in a file called LJ.in.
    
    ntypes
    A11 ... A1ntypes
    .   .
    .          .
    .                 .
    Antypes1 ... Antypes_ntypes
    B11 ... B1ntypes
    .   .
    .        .
    .	         .
    Bntypes1 ... Bntypes_ntypes


Both Aij and Bij matrices should be symmetric.


###Contents
0. Introduction
1. Installation instructions
2. Use
3. Contact
4. License and Copyright


##0. Introduction
**fluam** is a code for fluctuating hydrodynamics with immersed structures based
on the immersed boundary method. It offers fluid solvers for the compressible and
incompressible Navier-Stokes equations,
* **Staggered Schemes for Fluctuating Hydrodynamics**, F. Balboa Usabiaga, J. B. Bell, R. Delgado-Buscalioni, A. Donev, T. G. Fai, B. E. Griffith andC. S. Peskin. Multiscale Modeling & Simulation, **10** (3), 1369 (2012). [http://arxiv.org/abs/1108.5188]
 
and particle solvers for regimes ranging from the acoustic time scales to the Stokes limit,
* **Minimal model for acoustic forces on Brownian particles**, F. Balboa Usabiaga and R. Delgado-Buscalioni. Physical Review E, **88**, 063304 (2013). [http://arxiv.org/abs/1307.0702] 
* **Inertial coupling method for particles in an incompressible fluctuating fluid**, F. Balboa Usabiaga, R. Delgado-Buscalioni, B. E. Griffith and A. Donev. Computer Methods in Applied Mechanics and Engineering, **269**, 139 (2014). [http://arxiv.org/abs/1212.6427]
* **Brownian Dynamics without Green's Functions**, S. Delong, F. Balboa Usabiaga, R. Delgado-Buscalioni, B. E. Griffith and A. Donev. The Journal of Chemical Physics, **140** (13), 134110 (2014). [http://arxiv.org/abs/1401.4198]


##1. Installation instructions
a) You will need a NVIDIA GPU with compute capability 1.3
or higher to use **fluam**. You don't need any GPU to compile 
and modify the code.

b) Third-party software: you will need the CUDA compiler
nvcc, CUDA libraries, and the cutil.h files that you can obtain
with the NVIDIA SDK package.

c) Edit the file fluam/bin/MakefileHeader
to include the right path to the NVIDIA SDK files
and the HydroGrid code in case you have it. 
Set the right architecture for your GPU in 
"NVCCFLAGS".

d) Move to fluam/bin/ and type 
make

e) To speed up the compilation process see fluam/work/README

##2. Use
To run **fluam** type
fluam data.main

data.main is a file with the option for the simulation, look
fluam/bin/data.main for the options.


##3. Contact
If you find problems contact the owner of the project
https://github.com/fbusabiaga/fluam


##4. License and Copyright
Source code is available at: https://github.com/fbusabiaga/fluam

**fluam** is released under the terms of the GNU General Public License. See
"COPYING" for more details.

The source files included with this release are copyrighted by their
authors. See each file for more information.

