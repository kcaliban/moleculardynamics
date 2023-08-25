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
        ("o,output", "Final configuration output file", cxxopts::value<std::string>())
        ("x,trajectory", "Trajectory output file", cxxopts::value<std::string>())
        ("t,timestep", "Timestep in fs", cxxopts::value<double>()->default_value("0.01"))
        ("d,duration", "Duration of simulation in fs", cxxopts::value<double>()->default_value("100"))
        ("b,thermostat", "Use Berendsen thermostat (target temperature 0K)", cxxopts::value<bool>()->default_value("false"))
        ("tau", "Tau for Berendsen thermostat", cxxopts::value<double>()->default_value("2"))
        ("q,add_Q", "Add heat to the system", cxxopts::value<bool>()->default_value("false"))
        ("r,relax", "Relax time for heat capacity calculation in fs", cxxopts::value<double>()->default_value("1000"))
        ("m,measure", "Measurement time for heat capacity calculation in fs", cxxopts::value<double>()->default_value("500"))
        ("k,kelvin", "How much kelvin to add during add_q", cxxopts::value<double>()->default_value("100"))
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

    bool thermostat = result["thermostat"].as<bool>();
    bool add_Q = result["add_Q"].as<bool>();

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

    auto [names, positions, velocities]{read_xyz_with_velocities(result["file"].as<std::string>())};
    Atoms atoms(positions, velocities);
    atoms.m = 196.96657 * 103.6; // Mass of Au im g/mol times factor for time in fs

    if (result["velocities"].as<bool>()) {
        std::cout << "Ignoring velocities" << std::endl;
        atoms.velocities = Eigen::Array3Xd::Zero(3, atoms.nb_atoms());
    }

    double timestep = result["timestep"].as<double>();
    std::cout << "Using timestep: " << timestep << " fs" << std::endl;
    double total_time = result["duration"].as<double>();
    int nb_steps = total_time / timestep;

    NeighborList neighbor_list;
    double cutoff = 10; // Standard parameter of Ducastelle

    double target_temperature = 0; // 0K as target temperature
    double tau = result["tau"].as<double>();

    std::ofstream log(result["log"].as<std::string>(), std::ios_base::app);

    log << "Timestep: " << timestep << std::endl;

    double relax_time = result["relax"].as<double>();
    double measurement_time = result["measure"].as<double>();
    // For heat capacity calculation, we increase the temperature and let the system relax
    if (add_Q) {
        std::cout << "Adding heat to the system" << std::endl;
        auto kinetic_energy_before = kinetic_energy(atoms.velocities, atoms.m);
        auto temp_before = temperature(atoms);
        log << "Temperature before: " << temp_before << std::endl;
        auto targert_temperature_after = temp_before + result["kelvin"].as<double>();
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
    }

    log << "Duration: " << total_time << " fs" << std::endl;
    
    double average_temperature = 0;
    double average_total_energy = 0;
    double average_kinetic_energy = 0;
    double average_potential_energy = 0;
    double average_temperature_count = 0;

    std::ofstream traj(result["trajectory"].as<std::string>());

    // Write initial configuration
    write_xyz(traj, atoms);

    for (int step = 0; step < nb_steps; step++) {
        // Apply thermostat
        if (thermostat)
            berendsen_thermostat(atoms, target_temperature, timestep, tau);

        // Update positions and velocities using previous forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.m);

        // Update neighbors list
        neighbor_list.update(atoms, cutoff);

        // Calculate new forces given new positions
        double potential_energy = ducastelle(atoms, neighbor_list);

        // Update positions and velocities using new forces
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.m);

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
        // log << time << "\t" << total_energy  << "\t" << potential_energy << "\t" << kin_energy << "\t" << temp << std::endl;

        // Print progress
        auto percent = (step + 1) * 100 / nb_steps;
        std::cout << "\r" << percent << "%";
        std::cout.flush();
    }
    std::cout << std::endl;

    // If we are adding heat, we calculate the average temperature
    if (add_Q) {
        average_temperature /= average_temperature_count;
        log << "Average temperature: " << average_temperature << std::endl;
        log << "Average kinetic energy: " << average_kinetic_energy / average_temperature_count << std::endl;
        log << "Average potential energy: " << average_potential_energy / average_temperature_count << std::endl;
        log << "Average total energy: " << average_total_energy / average_temperature_count << std::endl;
    }

    log.close();
    traj.close();

    // Write final configuration
    std::ofstream final(result["output"].as<std::string>());
    write_xyz(final, atoms);
    final.close();

    return 0;
}
