#include "atoms.h"
#include "verlet.h"
#include "kinetic.h"
#include "lj_direct_summation.h"
#include "xyz.h"
#include "cxxopts.h"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {
    cxxopts::Options options("Milestone 4", "Molecular dynamics simulation");
    options.add_options()
        ("f,file", "Input file", cxxopts::value<std::string>())
        ("x,trajectory", "Trajectory output file", cxxopts::value<std::string>()->default_value("traj.xyz"))
        ("t,timestep", "Timestep in LJ unit sqrt(m * sigma * sigma / epsilon)", cxxopts::value<double>()->default_value("0.001"))
        ("d,duration", "Duration of simulation in LJ unit sqrt(m * sigma * sigma / epsilon)", cxxopts::value<double>()->default_value("100"))
        ("l,log", "Log file", cxxopts::value<std::string>()->default_value("log.txt"))
        ("h,help", "Print usage");

    auto result{options.parse(argc, argv)};

    if (result.count("help") > 0) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (result.count("file") == 0) {
        std::cout << "No input file given" << std::endl;
        return 1;
    }

    // Read atoms from file
    auto [names, positions, velocities]{read_xyz_with_velocities(result["file"].as<std::string>())};
    Atoms atoms(positions, velocities);

    // Initialize simulation parameters
    const double epsilon{1.0};
    const double sigma{1.0};
    const double m{1.0};
    const double total_time{result["duration"].as<double>() * sqrt(m * sigma * sigma / epsilon)};
    const double timestep{result["timestep"].as<double>() * sqrt(m * sigma * sigma / epsilon)};
    const int nb_steps{static_cast<int>(round(total_time / timestep))};

    // Open output files
    std::ofstream traj(result["trajectory"].as<std::string>());
    std::ofstream output(result["log"].as<std::string>());

    // Simulation loop
    for (int step = 0; step < nb_steps; step++) {
        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);

        // Calculate new forces given new positions
        const double potential_energy{lj_direct_summation(atoms, epsilon, sigma)};

        // Update velocities using new forces
        verlet_step2(atoms.velocities, atoms.forces, timestep);

        // Update time
        const double time{step * timestep};

        // Calculate energies and write to log file
        const double kin_energy{kinetic_energy(atoms.velocities, m)};
        const double total_energy{potential_energy + kin_energy};
        output << time << "\t" << total_energy  << "\t" << potential_energy << "\t" << kin_energy << std::endl;

        // Write trajectory every sqrt(m * sigma * sigma / epsilon)
        if (abs(ceil(time / sqrt(m * sigma * sigma / epsilon)) - (time / sqrt(m * sigma * sigma / epsilon))) < 0.0000001) {
            write_xyz(traj, atoms);
        }
    }

    // Close output files
    output.close();
    traj.close();

    return 0;
}
