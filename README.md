# LAMMPS_wrapper
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
cmake ../cmake -D BUILD_SHARED_LIBS=yes -D PKG_REPLICA=yes PKG_PHONON=yes
cmake --build . --target install
```
2. Build Check compiler type in `CMakeLists.txt`
```bash
cd LAMMPS_wrapper
mkdir build; cd build
cmake ../src
cmake --build .
```

## Input
### Potential parameter ###
* **NELEMENT** [integer]
  - *NELEMENT* is the number of elements
* **ATOM_TYPE** [strings]
  - *ATOM_TYPE* is the symbols of elements
* **PAIR_STYLE** [strings]
  - *PAIR_STYLE* stands for the pair style in LAMMPS input.
* **PAIR_COEFF** [strings]
  - *PAIR_COEFF* stands for the pair coeff in LAMMPS input.
### oneshot ###
* **ONESHOT** [0/1]
  - *ONESHOT* determines whether oneshot calculation runs.
### relaxation ###
* **ATOM_RELAX** [0/1]
  - *ATOM_RELAX* determines whether atomic coordinates are optimized.
* **CELL_RELAX** [0/1]
  - *CELL_RELAX* determines whether lattice vectors are also optimized.
* **MAX_FORCE** [real]
  - *MAX_FORCE* sets the force tolerance for relaxation (in eV/Angst).
### neb ###
* **NEB** [0/1]
  - *NEB* determines whether the nudged-elastic band method runs.
* **NIMAGES** [integer]
  - *NIMAGES* sets the number of images during the NEB emthod.
* **MIN_DIST** [real]
  - *MIN_DIST* sets the minimum bond distance between atoms when generating images with interpolation.
### dynamic matrix ###
* **DYNMAT** [0/1]
  - *DYNMAT* determines whether the dynamical matrix is calculated.
* **FINITE_DIFF** [real]
  - *FINITE_DIFF* sets the displacement in finite difference method (in Angst).
### NVT MD ###
* **NVT_MD** [0/1]
  - *NVT_MD* determines whether NVT molecular dynamics run.
* **TIMESTEP** [real]
  - *TIMESTEP* sets the timestep for time integration (in fs).
* **MAX_STEP** [integer]
  - *MAX_STEP* sets the maximum MD step.
* **TEMPERATURE** [real]
  - *TEMPERATURE* sets the temperature of system (in K).

## TARGET (only for DYNMAT)
It contains the atom indices or types to be the targets of dynamical matrix.
```text
I 0 1 2 3
T 1
A
```

* I: Index (starting from 0)
* T: Type (starting from 1)
* A: All

## Command
```bash
mpirun -np $numproc ./LAMMPS_wrapper POSCAR
```

For NEB calculation, two poscar files are needed.
```bash
mpirun -np $numproc ./LAMMPS_wrapper POSCAR_in POSCAR_fin
```
`$numproc` stands for the number of CPU cores in parallel computation.

## Tips  
1. `INIT_CONFIG` should be VASP5 POSCAR format. 
2. `numproc` in command should be the multiple of nimages. 
3. `TARGET` file is needed for calculating dynamical matrix. 
4. Output of dynamical matrix is written as `dynmat.dat`.
