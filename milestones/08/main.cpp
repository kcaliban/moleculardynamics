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
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cxxopts::Options options("Milestone 8", "Molecular dynamics simulation");
    options.add_options()
        ("f,file", "Input file", cxxopts::value<std::string>())
        ("o,output", "Final configuration output file", cxxopts::value<std::string>()->default_value("out.xyz"))
        ("t,timestep", "Timestep in fs", cxxopts::value<double>()->default_value("0.01"))
        ("d,duration", "Duration of simulation in fs", cxxopts::value<double>()->default_value("100"))
        ("l,log", "Log file", cxxopts::value<std::string>()->default_value("log.txt"))
        ("v,velocities", "Ignore velocities", cxxopts::value<bool>()->default_value("false"))
        ("thermostat", "Enable Berendsen thermostat", cxxopts::value<bool>()->default_value("false"))
        ("t0", "Target temperatures, if using Berendsen", cxxopts::value<double>()->default_value("300"))
        ("tau", "Relaxation time in fs if using Berendsen", cxxopts::value<double>()->default_value("1000"))
        ("gx", "Number of grid cells in x direction (if not provided, decomposition will be automatic if using 1, 2, 4 or 8 processes)", cxxopts::value<int>())
        ("gy", "Number of grid cells in y direction (if not provided, decomposition will be automatic if using 1, 2, 4 or 8 processes)", cxxopts::value<int>())
        ("gz", "Number of grid cells in z direction (if not provided, decomposition will be automatic if using 1, 2, 4 or 8 processes)", cxxopts::value<int>())
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

    if (!result.count("velocities")) {
        std::cout << "Velocities are taken from input file. If program crashes, check that velocities are correct." << std::endl;
    }

    if (result.count("gx") || result.count("gy") || result.count("gz")) {
        if (!result.count("gx") || !result.count("gy") || !result.count("gz")) {
            std::cout << "Please either provide number of grid cells in all directions, or none." << std::endl;
            return 1;
        }
    }

    // Read atoms from file
    auto [names, positions, velocities]{read_xyz_with_velocities(result["file"].as<std::string>())};
    Atoms atoms(positions, velocities);
    atoms.m = 196.96657 * 103.6; // Mass of Au im g/mol times factor for time in fs

    // Store total number of atoms for calculation of thermodynamic properties
    const size_t total_number_atoms = atoms.nb_atoms();

    // Calculate domain size
    const auto domain_buffer = 20; // spacing between system and domain boundaries
    Eigen::Vector3d max = atoms.positions.rowwise().maxCoeff();
    Eigen::Array3d domain_size = {max[0] + domain_buffer, max[1] + domain_buffer, max[2] + domain_buffer};

    // Determine domain decomposition depending on number of processes or use provided values
    Eigen::Array3i domain_decomposition;

    if (result.count("gx")) {
        const auto gx = result["gx"].as<int>();
        const auto gy = result["gy"].as<int>();
        const auto gz = result["gz"].as<int>();

        if (gx * gy * gz != size) {
            std::cout << "Please use more processes or change grid size! " << gx * gy * gz << " != " << size << std::endl;
            return 1;
        }

        domain_decomposition = {gx, gy, gz};
    } else {
        if (size == 1) {
            domain_decomposition = {1, 1, 1};
        } else if (size == 2) {
            domain_decomposition = {2, 1, 1};
        } else if (size == 4) {
            domain_decomposition = {2, 2, 1};
        } else if (size == 8) {
            domain_decomposition = {2, 2, 2};
        }
    }

    // Set domain periodicity
    Eigen::Array3i domain_periodicity = {0, 0, 0};
    
    // Initialize domain
    Domain domain(
        MPI_COMM_WORLD,
        domain_size,
        domain_decomposition,
        domain_periodicity);

    if (result["velocities"].as<bool>()) {
        std::cout << "Ignoring velocities" << std::endl;
        // Set velocities to zero
        atoms.velocities = Eigen::Array3Xd::Zero(3, atoms.nb_atoms());
    }

    // Initialize time parameters
    const double timestep = result["timestep"].as<double>();
    const double total_time = result["duration"].as<double>();
    const int nb_steps = total_time / timestep;

    // Initialize Berendsen parameters
    const bool berendsen = result["thermostat"].as<bool>();
    const double target_temperature = result["t0"].as<double>();
    const double tau = result["tau"].as<double>();

    // Open log file, log timestep
    std::ofstream log(result["log"].as<std::string>(), std::ios_base::app);

    // Write simulation info
    if (rank == 0) {
        log << "Using " << size << " processers" << std::endl;
        log << "Timestep: " << timestep << std::endl;
        log << "Duration: " << total_time << " fs" << std::endl;
    }

    // Initialize neighbor list
    NeighborList neighbor_list;
    const double cutoff = 10; // Standard parameter of Ducastelle

    // Switch to decomposed state
    domain.enable(atoms);

    // Simulation loop
    for (int step = 0; step < nb_steps; step++) {
        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.m);

        // Exchange atoms and update ghosts
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * cutoff);

        // Update neighbors list
        neighbor_list.update(atoms, cutoff);

        // Get potential energy
        const double potential_energy_total{
            MPI::allreduce(ducastelle_subset(atoms, neighbor_list, domain.nb_local()),
                           MPI_SUM,
                           MPI_COMM_WORLD)};

        // Update velocities using new forces
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.m);

        // Get kinetic energy
        const double kinetic_energy_total{
            MPI::allreduce(kinetic_energy_subset(atoms.velocities, atoms.m, domain.nb_local()),
                           MPI_SUM,
                           MPI_COMM_WORLD)};

        if (berendsen) {
            berendsen_thermostat(atoms, kinetic_energy_total, target_temperature, timestep, tau, total_number_atoms);
        }
        
        // Calculate current time for logging
        const double time = step * timestep;

        // Write log
        if (rank == 0) {
            const double total_energy_total{potential_energy_total + kinetic_energy_total};
            const double temp{temperature(kinetic_energy_total, total_number_atoms)};
            log << time << "\t" << total_energy_total  << "\t" << potential_energy_total << "\t" << kinetic_energy_total << "\t" << temp << std::endl;
        }

        // Every 10000fs, output number of local atoms per rank (as a sanity check)
        if (ceil(time) - time < 0.00001 && ((int) ceil(time)) % 10000 == 0) {
            std::cout << rank << ": " << domain.nb_local() << std::endl;
        }
    }

    // Close log file
    log.close();

    // Write final configuration and close file
    domain.disable(atoms);
    if (rank == 0) {
        std::ofstream final(result["output"].as<std::string>());
        write_xyz(final, atoms);
        final.close();
    }

    MPI_Finalize();
    return 0;
}
