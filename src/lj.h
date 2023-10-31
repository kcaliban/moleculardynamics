#ifndef LJ_H
#define LJ_H

#include "atoms.h"
#include "neighbors.h"

/**
 * Calculate the Lennard-Jones potential using neighbor lists
 * 
 * @param atoms Atoms to calculate potential for
 * @param neighbor_list Neighbor list
 * @param epsilon Lennard-Jones epsilon
 * @param sigma Lennard-Jones sigma
 * @return Potential energy (dimensionless)
*/
double lj(Atoms &atoms, NeighborList& neighbor_list, double epsilon, double sigma);

#endif // LJ_H