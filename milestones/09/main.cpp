#include "atoms.h"
#include "verlet.h"
#include "kinetic.h"
#include "temperature.h"
#include "berendsen.h"
#include "ducastelle.h"
#include "virial.h"
#include "xyz.h"
#include "cxxopts.h"
#include "volume.h"
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

    cxxopts::Options options("Milestone 9", "Molecular dynamics simulation");
    options.add_options()
        ("f,file", "Input file", cxxopts::value<std::string>())
        ("o,output", "Final configuration output file", cxxopts::value<std::string>()->default_value("out.xyz"))
        ("x,trajectory", "Trajectory output file", cxxopts::value<std::string>()->default_value("traj.xyz"))
        ("t,timestep", "Timestep in fs", cxxopts::value<double>()->default_value("0.01"))
        ("d,duration", "Duration of simulation in fs", cxxopts::value<double>()->default_value("100"))
        ("l,log", "Log file", cxxopts::value<std::string>()->default_value("log.txt"))
        ("v,velocities", "Ignore velocities", cxxopts::value<bool>()->default_value("false"))
        ("s", "Strain-stress output", cxxopts::value<std::string>()->default_value("strain-stress.txt"))
        ("tau", "Tau parameter for Berendsen thermostat", cxxopts::value<double>()->default_value("500"))
        ("temp", "Target temperature for Berendsen thermostat", cxxopts::value<double>()->default_value("300"))
        ("relax", "Relaxation time before straining", cxxopts::value<double>()->default_value("10000"))
        ("rate", "Power of strain rate, i.e. 10^{rate} * s^{-1}", cxxopts::value<double>()->default_value("8"))
        ("factor", "Factor of strain rate, i.e. factor * 10^{rate} * s^{-1}", cxxopts::value<double>()->default_value("1"))
        ("dx", "Size of domain in x-direction", cxxopts::value<double>()->default_value("40"))
        ("dy", "Size of domain in y-direction", cxxopts::value<double>()->default_value("40"))
        ("dz", "Size of domain in z-direction", cxxopts::value<double>()->default_value("145"))
        ("gx", "Number of grid cells in x direction (if not provided, use max. decomposition along z axis)", cxxopts::value<int>())
        ("gy", "Number of grid cells in y direction (if not provided, use max. decomposition along z axis)", cxxopts::value<int>())
        ("gz", "Number of grid cells in z direction (if not provided, use max. decomposition along z axis)", cxxopts::value<int>())
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
    // and to print warnings in case atoms fly out of the domain
    size_t total_number_atoms = atoms.nb_atoms();

    // Length of whisker_small: 144.24978336
    // X, Y of whisker_small: 40
    // Length of whisker_large: 288.49956672
    // X, Y of whisker_large: 84
    const double initial_length = result["dz"].as<double>();
    Eigen::Array3d domain_length = {result["dx"].as<double>(), result["dy"].as<double>(), initial_length};

    // Calculate domain decomposition, either decomposing along z axis
    // or using user-provided decomposition
    Eigen::Array3i domain_decomposition = {1, 1, size};
    if (result.count("gx")) {
        const auto gx = result["gx"].as<int>();
        const auto gy = result["gy"].as<int>();
        const auto gz = result["gz"].as<int>();

        if (gx * gy * gz != size) {
            std::cout << "Please use more processes or change grid size! " << gx * gy * gz << " != " << size << std::endl;
            return 1;
        }

        domain_decomposition = {gx, gy, gz};
    }

    // Set domain periodicity
    Eigen::Array3i domain_periodicity = {0, 0, 1};

    // Initialize domain
    Domain domain(
        MPI_COMM_WORLD,
        domain_length,
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
    const double target_temperature = result["temp"].as<double>();
    const double tau = result["tau"].as<double>();
    const double relaxation_time = result["relax"].as<double>(); // in fs

    // Initialize straining parameters
    const double scaling_fac = result["factor"].as<double>();
    const double scaling_rate = result["rate"].as<double>();
    // scaling rate is per second, one second = 10^{15}fs, 
    const double scaling_factor = scaling_fac * pow(10, scaling_rate - 15);

    // Open output files
    std::ofstream log(result["log"].as<std::string>(), std::ios_base::app);
    std::ofstream strain_stress_output(result["s"].as<std::string>(), std::ios_base::app);
    std::ofstream traj(result["trajectory"].as<std::string>());

    // Write initial configuration
    if (rank == 0) {
        log << "Using " << size << " processers" << std::endl;
        log << "Timestep: " << timestep << std::endl;
        log << "Duration: " << total_time << " fs" << std::endl;
        write_xyz(traj, atoms);
    }

    // Initialize neighbor list
    NeighborList neighbor_list;
    const double cutoff = 10.0; // Standard parameter of Ducastelle

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

        // Get potential energy and embeddings (the latter required for stress calculations) for local atoms
        const auto [potential_local, embeddings] = ducastelle_subset_embedding(atoms, neighbor_list, domain.nb_local());
        const double potential_energy_total{
            MPI::allreduce(potential_local,
                           MPI_SUM,
                           MPI_COMM_WORLD)};

        // Update velocities using new forces
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.m);

        // Get kinetic energy for local atoms
        const auto kinetic_local = kinetic_energy_subset(atoms.velocities, atoms.m, domain.nb_local());
        const double kinetic_energy_total{
            MPI::allreduce(kinetic_local,
                           MPI_SUM,
                           MPI_COMM_WORLD)};

        // Apply berendsen 
        berendsen_thermostat(atoms, kinetic_energy_total, target_temperature, timestep, tau, total_number_atoms);

        // Calculate current time, get current domain length
        const double time = step * timestep;
        const auto current_domain_length = domain.domain_length();

        // Write log
        if (rank == 0) {
            const double temp{temperature(kinetic_energy_total, total_number_atoms)};
            const double total_energy_total{potential_energy_total + kinetic_energy_total};
            log << time << "\t" << total_energy_total  << "\t" << potential_energy_total << "\t" << kinetic_energy_total << "\t" << temp << std::endl;
        }

        // Write strain-stress every step after relaxation time
        if (ceil(time) - time < 0.00001 && ((int) ceil(time)) > relaxation_time) {
            const double volume = current_domain_length[0] * current_domain_length[1] * current_domain_length[2];
            const auto [local_stress, local_stress_ghosts] = stress_z(atoms, volume, neighbor_list, embeddings, domain.nb_local());
            const double stress_local{
                    MPI::allreduce(local_stress,
                                   MPI_SUM,
                                   MPI_COMM_WORLD)
            };
            const double stress_ghost{
                    MPI::allreduce(local_stress_ghosts,
                                   MPI_SUM,
                                   MPI_COMM_WORLD)
            };
            // Ghost stresses are always double
            const double stress = stress_local + stress_ghost / 2;
            if (rank == 0) {
                double strain = (current_domain_length[2] - initial_length) / initial_length;
                strain_stress_output << time << "\t" << strain << "\t" << stress << std::endl;
            }
        }

        // Every 10000fs, output number of local atoms & ghost atoms per rank (as a sanity check)
        if (ceil(time) - time < 0.00001 && ((int) ceil(time)) % 10000 == 0) {
            std::cout << time << "\t" << rank  << ": " << domain.nb_local() << " (Ghost atoms: " << atoms.nb_atoms() - domain.nb_local() << ")" << std::endl;
        }

        // Write trajectory every 500fs
        if (ceil(time) - time < 0.00001 && ((int) ceil(time)) % 500 == 0) {
            domain.disable(atoms);
            if (rank == 0) {
                write_xyz(traj, atoms);
            }

            if (atoms.nb_atoms() != total_number_atoms) {
                std::cout << "Warning: Lost " << total_number_atoms - atoms.nb_atoms() <<  " atom(s)" << std::endl;
                total_number_atoms = atoms.nb_atoms();
            }
            domain.enable(atoms);
        }

        // Start scaling after relaxation time
        if (ceil(time) - time < 0.00001 && ((int) ceil(time)) >= relaxation_time) {
            const auto scaling = scaling_factor * (time - relaxation_time);
            domain.scale(atoms, {current_domain_length[0], current_domain_length[1], initial_length * (1 + scaling)});
        }
    }

    // Close output files
    log.close();
    traj.close();
    strain_stress_output.close();

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
