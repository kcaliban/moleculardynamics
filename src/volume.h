#ifndef VOLUME_H
#define VOLUME_H

#include <Eigen/Dense>
#include "atoms.h"

/** 
 * Calculate the volume of the simulation box using minimum and maximum coordinates
 * 
 * @param atoms Atoms to calculate the volume for
 * @return Volume of the simulation box
*/
double volume(const Atoms &atoms);

#endif // VOLUME_H