#ifndef LJ_DIRECT_SUMMATION_H
#define LJ_DIRECT_SUMMATION_H

#include "atoms.h"

/**
 * Calculate the Lennard-Jones potential using direct summation
 * 
 * @param atoms Atoms to calculate potential for
 * @param epsilon Lennard-Jones epsilon
 * @param sigma Lennard-Jones sigma
 * @return Potential energy (dimensionless)
 */
double lj_direct_summation(Atoms &atoms, double epsilon, double sigma);

#endif // LJ_DIRECT_SUMMATION_H