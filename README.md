
# gmxstrain

gmxstrain is a script for the calculation of finite temperature elastic constants of molecular crystals,
using stress/strain data from MD simulations in the NPT ensemble. The scripts generates input files
(strained simulation cells), executes (as external process) MD simulations with Gromacs, and analyzes
the resulting stress tensors. All time consuming parts of the calculations are performed by Gromacs, which
arguably is one of the fastest MD simulation engines out there, while the rest is book-keeping, and therefore
can be executed in a simple AWK script, as done here. The script is based loosely on a similar script that is part
of [LAMMPS](https://www.sciencedirect.com/science/article/pii/S0010465521002836) (more
details [here](https://github.com/lammps/lammps/tree/develop/examples/ELASTIC), while the accessory script
anisotropy.awk, that calculates moduli and anisotropy, is based on a matlab script
[published by Christopher Kube](https://aip.scitation.org/doi/10.1063/1.4962996)

## Features

- Based on the crystal structure and a classical molecular force field gmxstrain calculates the full elastic tensor for a given crystalline material
- Any force field that is supported by gromacs can be used.
- Based on the elastic tensor, bulk modulus, shear modulus, and the elastic anisotropy can be calculated.
- Since the stress tensors are determined in MD simulations, finite temperature effects are included.
- In our experience (see the paper referenced below) the accuracy of this approach for molecular crystals,
rivals that of DFT based calculations, at a fraction of the computational effort.

## Prerequisites

gmxstrain is a simple AWK script that should run on most Linux distributions, without
compilation, and out of the box. For the MD simulations you need a reasonably recent
version of Gromacs, with its binaries in your PATH.

- [Gromacs](http://www.gromacs.org/)

Other external tools that are not strictly required, but helpful, include:

- [Chimera](https://www.cgl.ucsf.edu/chimera/)
- [openbabel](http://openbabel.org/)
- [acpype](https://github.com/alanwilter/acpype)
- [Ambertools](http://ambermd.org/AmberTools.php)
- [gdis](https://github.com/arohl/gdis)

## Installation

```
git clone https://github.com/mbatgh/gmxstrain.git
```

This will leave a directory gmxstrain in your cwd, with a subdirectory scripts that contains the
script gmxstrain.awk, which can be directly executed on the linux command line.

## Usage

The top section of the script provides brief instructions:

```
# usage:
#
# gmxstrain.awk -v id=ID -v namx=NMAX -v dd=DD -v to=T0
#
# this expects a pdb file id.pdb and the corresponding gmx topology file named id.top to be present
# other variables to set on the command line:
# nmax ... number of different strain values
# dd ..... strain increment, stress tensor will be calculated at strain = -nmax*dd, ..., -dd, 0, dd, ..., nmax*dd
#          for each of 3 compressions and 3 shears.
# t0 ..... time (in ps) to discard from trajectory when calculating average stress tensors
#
# number of time steps to run each MD simulation and temperature need to be set in the mdp files
#
```

For more details see the documented [workflow](examples/README.md) in the examples folder.

## Limitations

- The constant pressure MD simulations use the Berendsen barostate, which does not produce
a (thermodynamically) correct NPT ensemble. One could use the Parrinello-Raman barostate which
is also implemented in Gromacs, but this can lead to large fluctuations of the volume, and numerical
instabilities. In our experience using Berendsen is a minor issue, and the resulting errors are
probably much smaller than those caused by the limitations of the classical force fields. However,
this has not been tested thoroughly, and we'd be happy to hear from you if you are aware of
any systematic study.

- In principle the implemented algorithm should work for any crystal space group. However,
this has not been tested thoroughly. In particular for tri-clinic structures we recommend careful
scrutiny of any results.

## Contact

Please send any questions/bug-reports/suggestions to me. My full name is rare enough, so that,
together with the place where I live (Graz), google will promptly provide you with an email address.

## Reference

If you use this code in published work, please cite the reference given below.

Michael Brunsteiner, Sten Nilsoon-Lill, Lucy M. Morgan, Adrian Davis, and Amrit Paudel
Finite temperature mechanical properties of organic molecular solids from classical molecular simulation.
To be submitted (2022).

## Acknowledgements

The work involved in the design of this tool was done at the [Research Center Pharamaceutical Engineering](http://www.rcpe.at),
in a collaborative K project, towards understanding chemical stability of small molecule drugs in the solid state, with funding
from the [Austrian COMET program](https://www.ffg.at/en/comet/programme), and a number of companies (Astra-Zeneca, Janssen,
Pfizer, and UCB).

