# Molecular Dynamics
Code for the [molecular dynamics course](https://pastewka.github.io/MolecularDynamics/) offered at the University of Freiburg. Report for the course can be found in `report.pdf`.

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