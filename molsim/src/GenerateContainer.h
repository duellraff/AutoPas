#ifndef SRC_GENERATECONTAINER_H_
#define SRC_GENERATECONTAINER_H_

#include "ParticleMS.h"
#include "../mdutils.h"
#include "ContainerProperties.h"

#include <assert.h>
#include <omp.h>
#include "MaxwellBoltzmannDistribution.h"
#include <list>
#include "utils/Vector.h"
#include <vector>
#include <array>
#include <cmath>

using namespace autopas;

	/** Generate a container
	 * @param autopasMS ParticleContainer to be used
	 * @param timestep time betweech each step
	 * @param BrownianFactor
	 * @param functor to be used
	 * @param Rcutoff only calculate forces between two particles if the distance is smaller than Rcutoff
	 * @param array with the dimensions of the box
	 * @param boundary boundary conditions to be applied
	 * @param threads number of threads used for calculation (not used)
	 * @return the properties of the generated Container
	 */

template<typename autopasClass, typename Functor>
ContProperties<autopasClass, Functor> Generate(autopasClass* autopasMS, double timestep,
					   double BrownianFactor, Functor* fm,
					   double Rcutoff, std::array<double, 3> Boxlimits,
					   std::array<int, 6> boundary,
					   int threads){

	std::array<double, 3> BoxMin({0,0,0});
	std::array<double, 3> BoxMax({Boxlimits});

	return ContProperties<autopasClass, Functor>(autopasMS, fm, boundary, timestep, BrownianFactor);

}

/**
 * generates a cuboid of particles
 * @param Container to which the particle shall be put
 * @param properties of the container
 * @param type type of each particle
 * @param pos coordinates of lower left front-side corner
 * @param dim number of particles in each dimension
 * @param h mesh width
 * @Param properties of the container
 * @param v velocity of each particle
 */
template<class autopasClass, class ContProp >
unsigned long generateCuboid(autopasClass* autopasMS, ContProp &prop, int type,
                    std::array<double,3> pos,
					std::array<int,3> dim,
					double h,
					std::array<double,3> v, unsigned long index) {

	// for each step in each dimension
	for (int x = 0; x < dim[0]; x++){
		for (int y = 0; y < dim[1]; y++) {
			for (int z = 0; z < dim[2]; z++) {

				// create a particle and add to container
				ParticleMS p = ParticleMS(type, index);

				std::array<double ,3> r  = {pos[0] + x*h, pos[1] + y*h, pos[2] + z*h};
				p.setR(r);

				p.setV(v);


				MaxwellBoltzmannDistribution(p, prop.getBF(), 3);
				autopasMS->addParticle(p);
				++index;
			}
		}
	}
	return index;
}


/**
 * generates a cuboid of membrane particles
 * @param autopasClass to which the particle shall be put
 * @param properties of the container
 * @param stiffness stiffness k of the membrane
 * @param r_zero average bond length of a molecule pair
 * @param t_end time when the force up should disappear
 * @param force force acting on certain particles in the membrane
 * @param force_coord list of mesh coordinates of the particles supporting an external force
 * @param type type of each particle
 * @param pos coordinates of lower left front-side corner
 * @param dim number of particles in each dimension
 * @param h mesh width
 * @Param properties of the container
 * @param v velocity of each particle
 */

template<class autopasClass, class ContProp >
unsigned long generateMembrane(autopasClass* autopasMS, ContProp &prop, int type, double stiffness,
				   double r_zero, double t_end,
				   std::array<double, 3> force,
				   std::list<std::array<int,3>> coord_force,
				   std::array<double,3> pos,
				   std::array<int,3> dim,
				   double h,
				   std::array<double,3> v, unsigned long index){

	//set the membrane parameter
	prop.setMembrane(stiffness, r_zero, true, t_end);

	// for each step in each dimension
	for (int z = 0; z < dim[2]; z++){

		int i = 0;

		std::vector<ParticleMS*> PartIndex;
		for (int y = 0; y < dim[1]; y++) {
			for (int x = 0; x < dim[0]; x++) {

				//get actual mesh coordinates in vector

				std::array<int, 3> meshcoord;
				meshcoord[0] = x;
				meshcoord[1] = y;
				meshcoord[2] = z;

				// create a particle and add to container
				ParticleMS p = ParticleMS(type, index);

				std::array<double ,3> r  = {pos[0] + x*h, pos[1] + y*h, pos[2] + z*h};
				p.setR(r);

				p.setV(v);

				MaxwellBoltzmannDistribution(p, prop.getBF(), 3);

				//add external force on certain particles
				for(std::list<std::array<int,3>>::iterator c = coord_force.begin(); c != coord_force.end(); ++c){
					if (*c == meshcoord){
						p.getConstF() = force;
					}
				}

				//add to the container
				autopasMS->addParticle(p);
				++index;
				PartIndex.push_back(&p);

				if(PartIndex[i] == nullptr){
					continue;
					i++;
				}
				ParticleMS& actual = *PartIndex[i];


				//complete neighbour vectors
				//all except the lower row
				if (i >= dim[0]){
					actual.getDirectNeighbours()[3] = PartIndex[i - dim[0]];
					actual.getDirectNeighbours()[3]->getDirectNeighbours()[1] = &actual;

					//ignore the left side
					if (i % dim[0] != 0){
						actual.getDirectNeighbours()[0] = PartIndex[i - 1];
						actual.getDiagonalNeighbours()[3] = PartIndex[i - dim[0] - 1];
						actual.getDirectNeighbours()[0]->getDirectNeighbours()[2] = &actual;
						actual.getDiagonalNeighbours()[3]->getDiagonalNeighbours()[1] = &actual;
					}
					//ignore the right side
					if ((i+1) % dim[0] != 0){
						actual.getDiagonalNeighbours()[2] = PartIndex[i - dim[0] + 1];
						actual.getDiagonalNeighbours()[2]->getDiagonalNeighbours()[0] = &actual;
					}
				}
				//fill the lower row
				else if (i > 0){
					actual.getDirectNeighbours()[0] = PartIndex[i - 1];
					actual.getDirectNeighbours()[0]->getDirectNeighbours()[2] = &actual;
				}

				++i;
			}
		}
	}
	return index;
}

/**
 * generates a sphere of particles
 * @param Container to which the particle shall be put
 * @param properties of the container
 * @param type type of each particle
 * @param pos coordinates of the center of sphere
 * @param r number of particles as radius
 * @param h mesh width
 * @Param properties of the container
 * @param m mass of each particle
 * @param v velocity of each particle
 */

template<typename autopasClass, class ContProp >
unsigned long generateSphere(autopasClass* autopasMS, ContProp &prop, int type,
                    std::array<double,3> pos,
					int r, double h,
					std::array<double,3> v, unsigned long index) {
	int r_sqr = r*r;

	// for each step in each dimension
	for (int x = -r; x < r; x++){
		int x_sqr = x*x;
		for (int y = -r; y < r; y++) {
			int y_sqr = y*y;
			for (int z = 0; z <= r; z++) {
				int z_sqr = z*z;

				if (x_sqr + y_sqr + z_sqr < r_sqr) {
					// create a particle and add to container
					ParticleMS p = ParticleMS(type, index);

					std::array<double ,3> r  = {pos[0] + x*h, pos[1] + y*h, pos[2] + z*h};
					p.setR(r);

					p.setV(v);

					MaxwellBoltzmannDistribution(p, prop.getBF(), 3);
					++index;
					autopasMS->addParticle(p);


				}
			}
		}
	}
	return index;
}


/**
 * add a thermostat which regulate the temperature to T
 * @param Container to which the thermostat shall be applied
 * @param T new temperature
 * @param ignoreY true if y-component is ignored
 */

template<typename autopasClass>
void setTemperature(autopasClass* autopasMS, float T_target, bool ignoreY){


	if (ignoreY) {
		double T = 0;
		int count = 0;

		// get sum of velocity squared in xz and count non-fix particles
		for (auto pp = autopasMS->begin(); pp.isValid(); ++pp ){

			ParticleMS& p = *pp;

			if (p.getFixed()) {
				continue;
			}

			std::array<double, 3> v = p.getV();
			double v_xz_sqr = v[0]*v[0] + v[2]*v[2];

			T += p.getM() * v_xz_sqr;
			count++;
		}

		T /= 3 * count;
		double beta = sqrt(T_target / T);

		// scale xz-velocity of all non-fix particles
		for (auto pp = autopasMS->begin(); pp.isValid(); ++pp) {
			ParticleMS& p = *pp;

			if (p.getFixed()) {
				continue;
			}

			std::array<double, 3> v = p.getV();
			v[0] *= beta;
			v[2] *= beta;
			p.setV(v);

		}
	}
	else
	{
		double T = 0;
		int count = 0;

		// get sum of velocity squared and count non-fix particles
		for (auto pp = autopasMS->begin(); pp.isValid(); ++pp) {
			ParticleMS& p = *pp;
			if (p.getFixed()) {
				continue;
			}

			std::array<double, 3> vel = p.getV();
			double Vnormsqr = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
			T += p.getM() * Vnormsqr;
			count++;
		}

		T /= 3 * count;
		double beta = sqrt(T_target / T);

		// scale velocity of all non-fix particles
		for (auto pp = autopasMS->begin(); pp.isValid(); ++pp) {
			ParticleMS& p = *pp;

			if (p.getFixed()) {
				continue;
			}

			p.getV() = arrayMath::mulScalar(p.getV(), beta);
		}
	}
}




#endif
