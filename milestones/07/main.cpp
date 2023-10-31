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

int main(int argc, char *argv[]) {
    cxxopts::Options options("Milestone 7", "Molecular dynamics simulation");
    options.add_options()
        ("f,file", "Input file", cxxopts::value<std::string>())
        ("o,output", "Final configuration output file", cxxopts::value<std::string>()->default_value("out.xyz"))
        ("x,trajectory", "Trajectory output file", cxxopts::value<std::string>()->default_value("traj.xyz"))
        ("t,timestep", "Timestep in fs", cxxopts::value<double>()->default_value("1"))
        ("d,duration", "Duration of simulation in fs", cxxopts::value<double>()->default_value("1000"))
        ("b,thermostat", "Use Berendsen thermostat (target temperature 0K)", cxxopts::value<bool>()->default_value("false"))
        ("tau", "Tau for Berendsen thermostat", cxxopts::value<double>()->default_value("2"))
        ("q,addQ", "Add heat to the system", cxxopts::value<bool>()->default_value("false"))
        ("r,relax", "Relax time for heat capacity calculation in fs", cxxopts::value<double>()->default_value("2500"))
        ("m,measure", "Measurement time for heat capacity calculation in fs", cxxopts::value<double>()->default_value("2500"))
        ("k,kelvin", "How much kelvin to add during add_q", cxxopts::value<double>()->default_value("100"))
        ("l,log", "Log file", cxxopts::value<std::string>()->default_value("log.txt"))
        ("et", "Energy-temperature file", cxxopts::value<std::string>()->default_value("et.txt"))
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

    if (!result.count("velocities")) {
        std::cout << "Velocities are taken from input file. If program crashes, check that velocities are correct." << std::endl;
    }

    const bool thermostat = result["thermostat"].as<bool>();
    const bool add_Q = result["addQ"].as<bool>();

    if (thermostat && add_Q) {
        std::cout << "Cannot use thermostat and add_Q at the same time" << std::endl;
        return 1;
    }

    if (!thermostat && result.count("tau") > 0) {
        std::cout << "Tau provided with disabled thermostat, ignoring" << std::endl;
    }

    if (!add_Q && result.count("relax") > 0) {
        std::cout << "Relax time provided with disabled add_Q, ignoring" << std::endl;
    }

    if (!add_Q && result.count("measure") > 0) {
        std::cout << "Measurement time provided with disabled add_Q, ignoring" << std::endl;
    }

    if (!add_Q && result.count("kelvin") > 0) {
        std::cout << "Kelvin provided with disabled add_Q, ignoring" << std::endl;
    }

    // Read atoms from file
    auto [names, positions, velocities]{read_xyz_with_velocities(result["file"].as<std::string>())};
    Atoms atoms(positions, velocities);
    atoms.m = 196.96657 * 103.6; // Mass of Au im g/mol times factor for time in fs

    if (result["velocities"].as<bool>()) {
        std::cout << "Ignoring velocities" << std::endl;
        // Set velocities to zero
        atoms.velocities = Eigen::Array3Xd::Zero(3, atoms.nb_atoms());
    }

    // Initialize time parameters (non-const, change for add_Q)
    double timestep = result["timestep"].as<double>();
    double total_time = result["duration"].as<double>();
    int nb_steps = total_time / timestep;
    std::cout << "Using timestep: " << timestep << " fs" << std::endl;

    // Initialize Berendsen parameters
    const double target_temperature = 0; // 0K as target temperature
    const double tau = result["tau"].as<double>();

    // Open log file, log timestep
    std::ofstream log(result["log"].as<std::string>(), std::ios_base::app);
    log << "Timestep: " << timestep << std::endl;

    // Initialize parameters for addQ
    const double relax_time = result["relax"].as<double>();
    const double measurement_time = result["measure"].as<double>();
    std::ofstream et;
    // For heat capacity calculation, we increase the temperature and let the system relax
    if (add_Q) {
        std::cout << "Adding heat to the system" << std::endl;
        const auto kinetic_energy_before = kinetic_energy(atoms.velocities, atoms.m);
        const auto temp_before = temperature(atoms);
        log << "Temperature before: " << temp_before << std::endl;
        const auto targert_temperature_after = temp_before + result["kelvin"].as<double>();
        atoms.velocities *= sqrt(targert_temperature_after / temp_before);
        log << "Temperature after: " << temperature(atoms) << std::endl;
        log << "Difference in temperature: " << temperature(atoms) - temp_before << std::endl << std::endl;
        log << "Kinetic energy before: " << kinetic_energy_before << std::endl;
        log << "Kinetic energy after: " << kinetic_energy(atoms.velocities, atoms.m) << std::endl;
        log << "Difference in kinetic energy: " << kinetic_energy(atoms.velocities, atoms.m) - kinetic_energy_before << std::endl << std::endl;
        log << "Relax time: " << relax_time << std::endl;
        log << "Measurement time: " << measurement_time << std::endl;
        total_time = relax_time + measurement_time;
        nb_steps = total_time / timestep;
    
        et.open(result["et"].as<std::string>(), std::ios_base::app);
    }

    // Log calculated duration
    log << "Duration: " << total_time << " fs" << std::endl;
    
    // Initialize variables used for averaging
    double average_temperature = 0;
    double average_total_energy = 0;
    double average_kinetic_energy = 0;
    double average_potential_energy = 0;
    double average_temperature_count = 0;

    // Open trajectory file
    std::ofstream traj(result["trajectory"].as<std::string>());

    // Write initial configuration
    write_xyz(traj, atoms);

    // Initialize neighbor list
    NeighborList neighbor_list;
    const double cutoff = 10; // Standard parameter of Ducastelle

    // Simulation loop
    for (int step = 0; step < nb_steps; step++) {
        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.m);

        // Update neighbors list
        neighbor_list.update(atoms, cutoff);

        // Calculate new forces given new positions
        double potential_energy = ducastelle(atoms, neighbor_list);

        // Update velocities using new forces
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.m);

        // Apply thermostat
        if (thermostat)
            berendsen_thermostat(atoms, target_temperature, timestep, tau);

        // Calculate energies and temperature
        double kin_energy = kinetic_energy(atoms.velocities, atoms.m);
        double total_energy = potential_energy + kin_energy;
        double temp = temperature(atoms);
        double time = step * timestep;

        // Write trajectory every 500 fs
        if (ceil(time) - time < 0.00001 && ((int) ceil(time)) % 500 == 0) {
            write_xyz(traj, atoms);
        }

        // If we are adding heat and we are in the measurement phase, we add the temperature to the average
        if (add_Q && time > relax_time) {
            average_temperature += temp;
            average_total_energy += total_energy;
            average_kinetic_energy += kin_energy;
            average_potential_energy += potential_energy;
            average_temperature_count++;
        }
        
        // Write data to log file
        log << time << "\t" << total_energy  << "\t" << potential_energy << "\t" << kin_energy << "\t" << temp << std::endl;

        // Print progress in percentage
        auto percent = (step + 1) * 100 / nb_steps;
        std::cout << "\r" << percent << "%";
        std::cout.flush();
    }
    // Newline after printing percentage
    std::cout << std::endl;

    // If we are adding heat, we calculate the average temperature
    if (add_Q) {
        average_temperature /= average_temperature_count;
        log << "Average temperature: " << average_temperature << std::endl;
        log << "Average kinetic energy: " << average_kinetic_energy / average_temperature_count << std::endl;
        log << "Average potential energy: " << average_potential_energy / average_temperature_count << std::endl;
        log << "Average total energy: " << average_total_energy / average_temperature_count << std::endl;
        et << average_total_energy / average_temperature_count << "\t" << average_temperature << std::endl;
    }

    // Close files
    log.close();
    traj.close();
    if (et.is_open()) et.close();

    // Write final configuration and close file
    std::ofstream final(result["output"].as<std::string>());
    write_xyz(final, atoms);
    final.close();

    return 0;
}
