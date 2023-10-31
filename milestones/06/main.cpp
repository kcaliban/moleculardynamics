#include "atoms.h"
#include "verlet.h"
#include "kinetic.h"
#include "lj.h"
#include "berendsen.h"
#include "temperature.h"
#include "xyz.h"
#include "cxxopts.h"
#include "cubic_lattice.h"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {
    cxxopts::Options options("Milestone 6", "Molecular dynamics simulation");
    options.add_options()
        ("a,atoms", "Number of atoms per side of the lattice", cxxopts::value<int>())
        ("x,trajectory", "Trajectory output file", cxxopts::value<std::string>()->default_value("traj.xyz"))
        ("c,constant", "Lattice constant, distance in LJ unit sigma", cxxopts::value<double>()->default_value("1"))
        ("b,tau", "Tau for Berendsen thermostat", cxxopts::value<double>()->default_value("2"))
        ("s,steps", "Equilibration steps, in each step, tau is doubled", cxxopts::value<int>()->default_value("5"))
        ("t,timestep", "Timestep in LJ unit sqrt(m * sigma * sigma / epsilon)", cxxopts::value<double>()->default_value("0.001"))
        ("d,duration", "Duration of each equilibration step in LJ unit sqrt(m * sigma * sigma / epsilon)", cxxopts::value<double>()->default_value("100"))
        ("l,log", "Log file", cxxopts::value<std::string>()->default_value("log.txt"))
        ("h,help", "Print usage");

    auto result{options.parse(argc, argv)};

    if (result.count("help") > 0) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (result.count("atoms") == 0) {
        std::cout << "Please provide the number of atoms per side for the lattice" << std::endl;
        return 1;
    }

    // Initialize parameters
    const double epsilon{1.0};
    const double sigma{1.0};

    // Initialize cubic lattice
    const int atoms_per_side{result["atoms"].as<int>()};
    const double lattice_spacing = result["constant"].as<double>() * sigma;
    auto atoms = create_cubic_lattice(lattice_spacing, atoms_per_side);
    atoms.m = 1.0;

    // Initialize time parameters
    const double total_time{result["steps"].as<int>() * result["duration"].as<double>() * sqrt(atoms.m * sigma * sigma / epsilon)};
    const double timestep{result["timestep"].as<double>() * sqrt(atoms.m * sigma * sigma / epsilon)};
    const int nb_steps{static_cast<int>(round(total_time / timestep))};

    // Initialize Berendsen parameters
    const double target_temperature{0};
    double tau{result["tau"].as<double>() * sqrt(atoms.m * sigma * sigma / epsilon)};

    // Open output files
    std::ofstream traj(result["trajectory"].as<std::string>());
    std::ofstream log(result["log"].as<std::string>());

    // Write initial configuration
    write_xyz(traj, atoms);

    // Initialize neighbor list
    const double cutoff{2.5 * sigma};
    NeighborList neighbor_list;

    for (int step = 0; step < nb_steps; step++) {
        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);

        // Update neighbor list (after positions have been updated)
        neighbor_list.update(atoms, cutoff);

        // Calculate new forces given new positions
        const double potential_energy{lj(atoms, neighbor_list, epsilon, sigma)};

        // Update velocities using new forces
        verlet_step2(atoms.velocities, atoms.forces, timestep);

        // Apply thermostat to new velocities
        berendsen_thermostat(atoms, target_temperature, timestep, tau, true);

        // Update time
        const double time{step * timestep};

        // Calculate energies and write to log file
        const double kin_energy{kinetic_energy(atoms.velocities, atoms.m)};
        const double total_energy{potential_energy + kin_energy};
        const double temp{temperature(atoms, true)};
        log << time << "\t" << total_energy  << "\t" << potential_energy << "\t" << kin_energy << "\t" << temp << "\t" << tau << std::endl;

        // Write trajectory every sqrt(m * sigma * sigma / epsilon)
        if (step != 0 && abs(ceil(time / sqrt(atoms.m * sigma * sigma / epsilon)) - (time / sqrt(atoms.m * sigma * sigma / epsilon))) < 0.0000001) {
            write_xyz(traj, atoms);
        }

        // After every duration, double the relaxation constant tau
        if (step != 0 && abs(ceil(time / (result["duration"].as<double>())) - (time / (result["duration"].as<double>()))) < 0.0000001) {
            tau *= 2;
        }
    }

    // Close output files
    traj.close();
    log.close();

    return 0;
}
