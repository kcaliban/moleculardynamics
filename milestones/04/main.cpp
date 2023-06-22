#include "atoms.h"
#include "verlet.h"
#include "kinetic.h"
#include "lj_direct_summation.h"
#include "xyz.h"
#include <iostream>
#include <fstream>

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[]) {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};

    Atoms atoms(positions, velocities);

    double epsilon = 1.0;
    double sigma = 1.0;
    double m = 1.0;
    double total_time = 100 * sqrt(m * sigma * sigma / epsilon);
    double timestep = 0.00001 * sqrt(m * sigma * sigma / epsilon);
    int nb_steps = total_time / timestep;

    // std::ofstream traj("trajectory.xyz");
    std::ofstream output("step_0_00001.txt");
    for (int step = 0; step < nb_steps; step++) {
        auto prev_forces_sum = atoms.forces.sum();

        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        verlet_step2(atoms.velocities, atoms.forces, timestep);

        // Calculate new forces given new positions
        double potential_energy = lj_direct_summation(atoms, epsilon, sigma);

        double time = step * timestep;
        double kin_energy = kinetic_energy(atoms.velocities, m);
        double total_energy = potential_energy + kin_energy;
        output << time << "\t" << total_energy  << "\t" << potential_energy << "\t" << kin_energy << std::endl;
    }
    output.close();

    // traj.close();

    return 0;
}
