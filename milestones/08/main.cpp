#include "atoms.h"
#include "verlet.h"
#include "kinetic.h"
#include "temperature.h"
#include "berendsen.h"
#include "ducastelle.h"
#include "xyz.h"
#include "cxxopts.h"
#include <iostream>
#include <fstream>

#include "mpi.h"
#include "mpi_support.h"
#include "domain.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cxxopts::Options options("Milestone 8", "Molecular dynamics simulation");
    options.add_options()
        ("f,file", "Input file", cxxopts::value<std::string>())
        ("o,output", "Final configuration output file", cxxopts::value<std::string>())
        ("x,trajectory", "Trajectory output file", cxxopts::value<std::string>())
        ("t,timestep", "Timestep in fs", cxxopts::value<double>()->default_value("0.01"))
        ("d,duration", "Duration of simulation in fs", cxxopts::value<double>()->default_value("100"))
        ("l,log", "Log file", cxxopts::value<std::string>()->default_value("log.txt"))
        ("v,velocities", "Ignore velocities", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help") > 0) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (result.count("file") == 0) {
        std::cout << "No input file given" << std::endl;
        return 1;
    }

    if (result.count("output") == 0) {
        std::cout << "No output file given" << std::endl;
        return 1;
    }

    if (result.count("trajectory") == 0) {
        std::cout << "No trajectory file given" << std::endl;
        return 1;
    }

    auto [names, positions, velocities]{read_xyz_with_velocities(result["file"].as<std::string>())};
    Atoms atoms(positions, velocities);
    atoms.m = 196.96657 * 103.6; // Mass of Au im g/mol times factor for time in fs

    Domain domain(
        MPI_COMM_WORLD,
        {40.0, 40.0, 40.0},
        {size, 1, 1},
        {0, 0, 0});

    if (result["velocities"].as<bool>()) {
        std::cout << "Ignoring velocities" << std::endl;
        atoms.velocities = Eigen::Array3Xd::Zero(3, atoms.nb_atoms());
    }

    double timestep = result["timestep"].as<double>();
    double total_time = result["duration"].as<double>();
    int nb_steps = total_time / timestep;

    NeighborList neighbor_list;
    double cutoff = 10; // Standard parameter of Ducastelle

    std::ofstream log(result["log"].as<std::string>(), std::ios_base::app);

    std::ofstream traj(result["trajectory"].as<std::string>());

    // Write initial configuration
    if (rank == 0) {
        log << "Using " << size << " processers" << std::endl;
        log << "Timestep: " << timestep << std::endl;
        log << "Duration: " << total_time << " fs" << std::endl;
        write_xyz(traj, atoms);
    }

    // Switch to decomposed state
    domain.enable(atoms);

    for (int step = 0; step < nb_steps; step++) {
        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.m);

        // Exchange atoms and update ghosts
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, cutoff);

        // Update neighbors list
        neighbor_list.update(atoms, cutoff);

        // Get potential energy
        double potential_energy_total{
            MPI::allreduce(ducastelle_subset(atoms, neighbor_list, domain.nb_local()),
                           MPI_SUM,
                           MPI_COMM_WORLD)};

        // Update positions and velocities using new forces
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.m);

        double time = step * timestep;

        if (ceil(time) - time < 0.00001 && ((int) ceil(time)) % 10 == 0) {
            domain.disable(atoms);

            if (rank == 0) {
                write_xyz(traj, atoms);
                double kin_energy_total = kinetic_energy(atoms.velocities, atoms.m);
                double total_energy_total = potential_energy_total + kin_energy_total;
                log << time << "\t" << total_energy_total  << "\t" << potential_energy_total << "\t" << kin_energy_total << "\t" << std::endl;
            }

            domain.enable(atoms);
        }
    }

    log.close();
    traj.close();

    // Write final configuration
    if (rank == 0) {
        std::ofstream final(result["output"].as<std::string>());
        write_xyz(final, atoms);
        final.close();
    }

    MPI_Finalize();
    return 0;
}
