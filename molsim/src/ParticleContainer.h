/*
 * ParticleContainer.h
 *
 *  Created on: Oct 31, 2017
 *      Author: raffi
 */

#ifndef SRC_PARTICLECONTAINER_H_
#define SRC_PARTICLECONTAINER_H_

#include "Particle.h"

#include <list>
#include "utils/Vector.h"


/** @brief Class to hold a set of particles
 *
 */
class ParticleContainer {
public:
	enum ForceMethod{
		Gravity,
		LennardJones
	};

protected:
	//list containing all particles
	std::list<Particle> partlist;
	double delta_t;
	ForceMethod fm;
	utils::Vector<double,3> (ParticleContainer::*calcFFunction)(Particle&, Particle&);

	double BrownianFactor;

	//membrane data
	double stiffness = 0;
	double r0 = 0;
	bool membrane_active = false;
	double t_end_force = 0;


public:

	/** Constructor
	 * @param f Function function to calculate the force between two particles
	 * @param timestep step between each step
	 * @param BrownianFactor the brownian Factor
	 */
//	ParticleContainer(ForceMethod fm, double timestep=0.014, double BrownianFactor=0.1);
//	virtual ~ParticleContainer();

	/** set the timestep
	 * @param timestep time between steps
	 */
//	void setDeltaT(double timestep);

	/** set the gravitational force
	 * @param g the gravitational acceleration
	 */
//	virtual void setGravity(utils::Vector<double,3> g);
	/** sets the membrane parameters
	 * @param stiff membrane stiffness
	 * @param r_zero average bond lenght
	 * @param status whether the membrane is active or not
	 * @param tend force t_end
	 */
//	virtual void setMembrane(double stiff, double r_zero, bool status, double tend);

	/** define all particles in the container
	 * @param plist the list to fill container
	 */
//	virtual void setList(std::list<Particle> plist);
	/** add one particle to the container
	 * @param p the particle to add
	 * @return a pointer to the newly created particle
	*/
//	virtual Particle* addParticle(Particle& p);
	/** remove a particle from the container
	 * @param p the particle to remove
	 */
//	virtual void deleteParticle(Particle& p);

	/** tests if the force up should be switched of
	 * @param time current simulation time
	 */
//	virtual void testForceEnd(double time);


	/** get a list of all particle in the container
	 * @return list of all particles
	 */
//	virtual std::list<Particle>& getAll();
	/** get a list of pointer to all particle in the container
	 * @return list of all particles
	 */
//	virtual std::list<Particle*> getAllPointer();
	/** set the force to all particles according to gravity force
	 */
	virtual void fReset();

	/**
	 * calculate the force for all particles
	 */
//	virtual void calculateF();

	/**
	 * calculate the position for all particles
	 * using the Velocity-Stoermer-Verlet-Algorithm
	 */
	virtual void calculateX();

	/**
	 * calculate the position for all particles
	 * using the Velocity-Stoermer-Verlet-Algorithm
	 */
//	virtual void calculateV();

	/**
	 * calculate forces using selected function
	 * @param p1 first particle
	 * @param p2 second particle
	 * @return vector containing forces
	 */
//	utils::Vector<double,3> calcF(Particle& p1, Particle& p2);

	/**
	 * Calculate forces acting in a membrane
	 * @param p particle
	 * @param k Stiffness of the membrane
	 * @param r0 average bond length
	 */
//	void HarmonicPotential(double k, double r_zero);
	/**
	 * calculate forces using gravity
	 * @param p1 first particle
	 * @param p2 second particle
	 * @return vector containing forces
	 */
//	utils::Vector<double,3> calcGravity(Particle& p1, Particle& p2);

	/**
	 * calculate forces using Lennard-Jones Potential
	 * @param p1 first particle
	 * @param p2 second particle
	 * @return vector containing forces
	 */
//	utils::Vector<double,3> calcLennardJones(Particle& p1, Particle& p2);

	/**
	 * generates a cuboid of particles
	 * @param type type of each particle
	 * @param pos coordinates of lower left front-side corner
	 * @param dim number of particles in each dimension
	 * @param h mesh width
	 * @param v velocity of each particle
	 */
//	void generateCuboid(int type,
	                    utils::Vector<double,3> pos,
						utils::Vector<int,3> dim,
						double h,
						utils::Vector<double,3> v);

	/**
	 * generates a cuboid of membrane particles
	 * @param stiffness stiffness k of the membrane
	 * @param r_zero average bond length of a molecule pair
	 * @param t_end time when the force up should disappear
	 * @param force force acting on certain particles in the membrane
	 * @param force_coord list of mesh coordinates of the particles supporting an external force
	 * @param type type of each particle
	 * @param pos coordinates of lower left front-side corner
	 * @param dim number of particles in each dimension
	 * @param h mesh width
	 * @param v velocity of each particle
	 */
//	void generateMembrane(int type, double stiffness,
					   double r_zero, double t_end,
					   utils::Vector<double, 3> force,
					   std::list<utils::Vector<int,3>> coord_force,
					   utils::Vector<double,3> pos,
					   utils::Vector<int,3> dim,
					   double h,
					   utils::Vector<double,3> v);

	/**
	 * generates a sphere of particles
	 * @param type type of each particle
	 * @param pos coordinates of the center of sphere
	 * @param r number of particles as radius
	 * @param h mesh width
	 * @param m mass of each particle
	 * @param v velocity of each particle
	 */
//	void generateSphere(int type,
	                    utils::Vector<double,3> pos,
						int r, double h,
						utils::Vector<double,3> v);

	/**
	 * add a thermostat which regulate the temperature to T
	 * @param T new temperature
	 * @param ignoreY true if y-component is ignored
	 */
//	void setTemperature(float T, bool ignoreY);
};

//inline utils::Vector<double,3> ParticleContainer::calcF(Particle& p1, Particle& p2){
	return (this->*calcFFunction)(p1, p2);
}


#endif /* SRC_PARTICLECONTAINER_H_ */
