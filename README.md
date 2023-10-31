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

## Milestone functionality
4. (Serial) Simple direct-summation Lennard-Jones simulation
5. (Serial) Lennard-Jones for arbitrary cubic lattice systems, with support for Berendsen thermostat
6. (Serial) Lennard-Jones with cutoff for arbitrary cubic lattice systems, with support for Berendsen thermostat
7. (Serial) Embedded-atom potential with support for heating a system
8. (Parallel) Embedded-atom potential with MPI support
9. (Parallel) Embedded-atom potential with MPI support and gathering stress-strain data