#include "atoms.h"
#include "verlet.h"
#include "kinetic.h"
#include "lj.h"
#include "berendsen.h"
#include "temperature.h"
#include "xyz.h"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " Number of atoms per lattice side" << std::endl;
        return 1;
    }

    int atoms_per_side = std::stoi(argv[1]);

    Atoms atoms(pow(atoms_per_side, 3));
    atoms.m = 1.0;
    atoms.velocities.setRandom();

    double cutoff = 0.5;
    double epsilon = 1.0;
    double sigma = 1.0;
    double total_time = 50 * sqrt(atoms.m * sigma * sigma / epsilon);
    double timestep = 0.001 * sqrt(atoms.m * sigma * sigma / epsilon);
    int nb_steps = total_time / timestep;

    double target_temperature = 0;
    // big tau => slow thermostat
    double tau = 2 * sqrt(atoms.m * sigma * sigma / epsilon);

    // Initialize cubic lattice
    double lattice_spacing = 1.0 * sigma;
    for (int i = 0; i < atoms_per_side; i++) {
        for (int j = 0; j < atoms_per_side; j++) {
            for (int k = 0; k < atoms_per_side; k++) {
                int index = i * atoms_per_side * atoms_per_side + j * atoms_per_side + k;
                atoms.positions.col(index)(0) = i * lattice_spacing;
                atoms.positions.col(index)(1) = j * lattice_spacing;
                atoms.positions.col(index)(2) = k * lattice_spacing;
            }
        }
    }

    /*
    std::ofstream lattice("lattice.xyz");
    write_xyz(lattice, atoms);

    std::ofstream traj("trajectory.xyz");
    write_xyz(traj, atoms);
    */

    // std::ofstream output("step_0_001.txt");
    for (int step = 0; step < nb_steps; step++) {
        
        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        verlet_step2(atoms.velocities, atoms.forces, timestep);

        // Apply thermostat
        berendsen_thermostat(atoms, target_temperature, timestep, tau);

        // Calculate new forces given new positions
        double potential_energy = lj(atoms, epsilon, sigma, cutoff);

        double kin_energy = kinetic_energy(atoms.velocities, atoms.m);
        double total_energy = potential_energy + kin_energy;
        double time = step * timestep;
        double temp = temperature(atoms);

        /*
        if (step % 100 == 0)
            write_xyz(traj, atoms);
        */

        // Adjust tau when close to equilibrium
        if (abs(temp - target_temperature) < 0.01)
            tau = 4 * sqrt(atoms.m * sigma * sigma / epsilon);

        // std::cout << "\t" << time << "\t" << temp << "\t" << total_energy  << "\t" << potential_energy << "\t" << kin_energy << std::endl;
    }
    // output.close();
    // traj.close();

    return 0;
}
