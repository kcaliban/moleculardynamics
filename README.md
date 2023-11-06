# Molecular Dynamics
Code for the [molecular dynamics course](https://pastewka.github.io/MolecularDynamics/) offered at the University of Freiburg. Final report for the course can be found in `report.pdf`.

To build and test:
```bash
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
```

To display all possible command-line arguments for milestone `X`:
```bash
./build/milestones/0X/0X -h
```

## Milestone functionality and commands
4. (Serial) Simple direct-summation Lennard-Jones simulation
```bash
./build/milestones/04/04 -f ./build/milestones/04/lj54.xyz
```
5. (Serial) Lennard-Jones for arbitrary cubic lattice systems, with support for Berendsen thermostat
```bash
./build/milestones/05/05 -a 4
```
6. (Serial) Lennard-Jones with cutoff for arbitrary cubic lattice systems, with support for Berendsen thermostat
```bash
./build/milestones/06/06 -a 4
```
7. (Serial) Embedded-atom potential with support for heating a system.
```bash
# Copy `heat_cap.sh` manually from `./milestones/07` to `./build/milestones/07` (CMake changes the contents of the file).
# Change sz in heat_cap.sh to the wished size
./build/milestones/07/heat_cap.sh
```
8. (Parallel) Embedded-atom potential with MPI support
```bash
mpirun -n 8 ./build/milestones/08/08 -f ./build/milestones/08/cluster_3871_transformed.xyz -o cluster_3871_transformed_out_8.xyz -t 10 -d 100000 -l cluster_3871_transformed_log_8.txt -v --thermostat --t0 300 --tau 1000 
# Or copy `experiments_no_thermo.sh`/`experiments_thermo.sh` from `./milestones/08` to `./build/milestones/08` and run these
```
9. (Parallel) Embedded-atom potential with MPI support and gathering stress-strain data
```bash
# Copy `experiments_small.sh`/`experiments_large` from `./milestones/09` to `./build/milestones/09`
./build/milestones/09/experiments_small.sh
./build/milestones/09/experiments_large.sh
```