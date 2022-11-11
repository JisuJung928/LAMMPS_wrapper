# LAMMPS_calculator
Any kinds of interatomic potential supported by LAMMPS can be used.  
Neural Network Potentials (NNPs) made by [SIMPLE-NN](https://github.com/MDIL-SNU/SIMPLE-NN_v2) are also possible.  

## Requirement
- CMake >= 3.13
- LAMMPS >= 23Jun2022

## Installation
1. Build LAMMPS as shared library. [[link](https://docs.lammps.org/Build_basics.html)]
```bash
cd lammps
mkdir build; cd build
cmake ../cmake -D BUILD_SHARED_LIBS=yes
cmake --build . --target install
```
2. Build Check compiler type in `CMakeLists.txt`
```bash
cd LAMMPS_calculator
mkdir build; cd build
cmake ../src
cmake --build .
```

## Input
```bash
# potential parameter #
NELEMENT    = 2
ATOM_TYPE   = O Pt
PAIR_STYLE  = nn
PAIR_COEFF  = * * potential_saved O Pt
MAX_FORCE   = 0.05

# calculator #
ONESHOT     = 0
ATOM_RELAX  = 1
CELL_RELAX  = 0
NEB         = 0

# neb parameter #
NIMAGES     = 7
```

|Tag|Description|Units|
|:---|:---|:---|
|NELEMENT|The number of elements||
|ATOM_TYPE|Symbols of atom types in potential||
|PAIR_STYLE|`pair_style` in LAMMPS input file||
|PAIR_COEFF|`pair_coeff` in LAMMPS input file||
|MAX_FORCE|Force convergence criteria|eV/Ang|
|MIN_DIST|Minimum distance in initial images of NEB|Ang|
|ONESHOT|Flag for oneshot calculation||
|ATOM_RELAX|Flag for atomic relaxation||
|CELL_RELAX|Flag for cell relaxation||
|NEB|Flag for nudged elastic band calculation||
|NIMAGES|The number of images in diffusion path||

## Command
```bash
mpirun -np $numproc ./LAMMPS_calculator POSCAR
mpirun -np $numproc ./LAMMPS_calculator POSCAR_in POSCAR_fin
```
`$numproc` stands for the number of CPU cores in parallel computation.

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format. 
2. `numproc` in command should be the multiple of nimages. 
