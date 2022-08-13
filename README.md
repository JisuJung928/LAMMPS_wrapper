# LAMMPS_NEB
Any kinds of interatomic potential supported by LAMMPS can be used.  
Neural Network Potentials (NNPs) made by [SIMPLE-NN](https://github.com/MDIL-SNU/SIMPLE-NN_v2) are also possible.  

## Requirement
- CMake >= 2.8.12
- LAMMPS == 29Oct2020

## Installation
1. Build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
```bash
cd /path/to/lammps/src
make yes-replica
make mpi
make mode=shlib mpi
```
2. Set the environment variable. [[link](https://docs.lammps.org/Build_link.html)]
```bash
export LD_LIBRARY_PATH=/path/to/lammps/src:$LD_LIBRARY_PATH
```
3. Modify `LMP_PATH` in `CMakeLists.txt`.
```text
SET ( LMP_PATH /path/to/lammps )
```
4. Make `Makefile` and build
``` bash
cmake CMakeLists.txt
make
```

## Input
```bash
# potential parameter #
NELEMENT        = 2
ATOM_TYPE       = O Pt
PAIR_STYLE      = nn
PAIR_COEFF      = * * potential_saved O Pt
CUTOFF          = 6.0

# neb parameter #
NIMAGES         = 7
MAX_FORCE       = 0.05

# structure parameter #
INIT_CONFIG     = ./POSCAR
INIT_RELAX      = 1
```

|Tag|Description|Units|
|:---|:---|:---|
|NELEMENT|The number of elements||
|ATOM_TYPE|Symbols of atom types in potential||
|PAIR_STYLE|`pair_style` in LAMMPS input file||
|PAIR_COEFF|`pair_coeff` in LAMMPS input file||
|CUTOFF|Cutoff of local energy|Ang|
|NIMAGES|The number of images in diffusion path||
|MAX_FORCE|Force convergence criteria|eV/Ang|
|INIT_CONFIG|Initial configuration in simulation||
|INIT_RELAX|Initial relaxation after preprocess (1: yes, 0: no)||

## Command
```bash
mpirun -np $numproc ./LAMMPS_NEB POSCAR_in POSCAR_fin index
```
$numproc stands for the number of CPU cores in parallel computation.
`index` represents the index of diffusing atom.

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format. 
2. `numproc` in command should be the multiple of nimages. 
