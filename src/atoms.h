#ifndef ATOMS_H
#define ATOMS_H

#include "types.h"
 
class Atoms { 
public: 
    Positions_t positions; 
    Velocities_t velocities; 
    Forces_t forces; 
 
    Atoms(const Positions_t &p) : 
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} { 
        velocities.setZero(); 
        forces.setZero(); 
    } 
 
    Atoms(const Positions_t &p, const Velocities_t &v) : 
            positions{p}, velocities{v}, forces{3, p.cols()} { 
        assert(p.cols() == v.cols());
        forces.setZero(); 
    } 

    Atoms(const int nb_atoms) :
        positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms} {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
    }

    size_t nb_atoms() const { 
        return positions.cols(); 
    }
};

#endif  // ATOMS_H