# LAMMPS_wrapper
**LAMMPS_wrapper** makes LAMMPS calculation easier by using POSCAR directly with simple INPUT files.
If you activate the calculation type you want to do, **LAMMPS_wrapper** converts INPUT and POSCAR files internally.
**LAMMPS_wrapper** supports oneshot calculation, relaxation, NEB method, dynamical matrix calculation, and NVT MD simulation.
Be free from writing complex input files and tedious conversion between POSCAR and LAMMPS data format!

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
2. Build `LAMMPS_wrapper`.
Please check compiler type in `CMakeLists.txt`.
```bash
cd LAMMPS_wrapper
mkdir build; cd build
cmake ../src
cmake --build .
```

## Usage
For oneshot, relaxation, NEB, dynamical matrix, and NVT MD simulation, single POSCAR file is needed.
```bash
mpirun -np $numproc ./LAMMPS_wrapper POSCAR
```

For NEB calculation, two POSCAR files are needed.
```bash
mpirun -np $numproc ./LAMMPS_wrapper POSCAR_initial POSCAR_final
```
`$numproc` stands for the number of CPU cores in parallel computation.

### Input
#### Potential parameter ####
* **NELEMENT** [integer]
  - *NELEMENT* is the number of elements
* **ATOM_TYPE** [strings]
  - *ATOM_TYPE* is the symbols of elements
* **PAIR_STYLE** [strings]
  - *PAIR_STYLE* stands for the pair style in LAMMPS input.
* **PAIR_COEFF** [strings]
  - *PAIR_COEFF* stands for the pair coeff in LAMMPS input.
#### oneshot ####
* **ONESHOT** [0/1]
  - *ONESHOT* determines whether oneshot calculation runs.
#### relaxation ####
* **ATOM_RELAX** [0/1]
  - *ATOM_RELAX* determines whether atomic coordinates are optimized.
* **CELL_RELAX** [0/1]
  - *CELL_RELAX* determines whether lattice vectors are also optimized.
* **MAX_FORCE** [real]
  - *MAX_FORCE* sets the force tolerance for relaxation (in eV/Angst).
#### neb ####
* **NEB** [0/1]
  - *NEB* determines whether the nudged-elastic band method runs.
* **NIMAGES** [integer]
  - *NIMAGES* sets the number of images during the NEB emthod.
* **MIN_DIST** [real]
  - *MIN_DIST* sets the minimum bond distance between atoms when generating images with interpolation.
#### dynamic matrix ####
* **DYNMAT** [0/1]
  - *DYNMAT* determines whether the dynamical matrix is calculated.
* **FINITE_DIFF** [real]
  - *FINITE_DIFF* sets the displacement in finite difference method (in Angst).
#### NVT MD ####
* **NVT_MD** [0/1]
  - *NVT_MD* determines whether NVT molecular dynamics run.
* **TIMESTEP** [real]
  - *TIMESTEP* sets the timestep for time integration (in fs).
* **MAX_STEP** [integer]
  - *MAX_STEP* sets the maximum MD step.
* **TEMPERATURE** [real]
  - *TEMPERATURE* sets the temperature of system (in K).

### TARGET (only for DYNMAT)
It contains the atom indices or types to be the targets of dynamical matrix.
```text
I 0 1 2 3
T 1
A
```

* I: Index (starting from 0)
* T: Type (starting from 1)
* A: All

## Tips  
1. **LAMMPS_wrapper** supports only `units metal` currently.
2. `TARGET` file is only needed for calculating dynamical matrix. 
3. Output of dynamical matrix is written as `dynmat.dat`.
