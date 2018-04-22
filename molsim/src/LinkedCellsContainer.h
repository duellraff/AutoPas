/*
 * LinkedCellsContainer.h
 *
 *  Created on: Nov 27, 2017
 *      Author: raffi
 */

#ifndef SRC_LINKEDCELLSCONTAINER_H_
#define SRC_LINKEDCELLSCONTAINER_H_
#include "Particle.h"
#include "ParticleContainer.h"

#include <list>
#include "utils/Vector.h"
#include <vector>


/** @brief Class to hold cells containing a set of particles
 *
 */
class LinkedCellsContainer:public ParticleContainer {

private:


	//number of cells in each dimension
	int NcellsX;
	int NcellsY;
	int NcellsZ;
	int NcellsXY;
	int NcellsTot;

	// 0=x; 1=y; 2=z
	int MaxDim;
	int NumThreads;

	//mesh width in each dimension
	double HcellsX;
	double HcellsY;
	double HcellsZ;

	//domain limits (3D)
	utils::Vector<double,3> lowerleft;
	utils::Vector<double,3> upperright;

	//vector containing all cells
	std::vector<ParticleContainer> celllist;

	//boundary conditions 1 for reflecting cells, 0 (or anything else) for outflow, order is left, top, right, bottom
	utils::Vector<int,6> BoundCond;


public:

	/** Constructor
	 * @param f Function function to calculate the force between two particles
	 * @param timestep time betweech each step
	 * @param Rcutoff only calculate forces between two particles if the distance is smaller than Rcutoff
	 * @param domainX x dimension of domain
	 * @param domainY y dimension of domain
	 * @param domainZ z dimension of domain
	 * @param boundary boundary conditions to be applied
	 * @param threads number of threads used for calculation
	 */
	LinkedCellsContainer(ForceMethod fm, double timestep, double BrownianFactor,
						double Rcutoff, double domainX, double domainY, double domainZ,
						 utils::Vector<int,6> boundary,
						 int threads=1);
	virtual ~LinkedCellsContainer();

	/** set the gravitational force along the z-axis
	 * @param g the gravitational acceleration
	 */
	void setGravity(utils::Vector<double,3> g);

	/**set domain
	 * @param x position vector of the lower left front-side corner
	 * @param domainX X-size of the domain
	 * @param domainY Y-size of the domain
	 */
	void setDomain(double domainX, double domainY, double domainZ);
	/** set the cell mesh
	 * @param Rcutoff the cut-off radius in the simulation
	 * @param Domain the domain for our simulation
	 */
	void setMesh(double Rcut, double domainX, double domainY, double domainZ);
	/** sets he membrane parameters to all cells
	 * @param stiff membrane stiffness
	 * @param r_ero average bond lenght
	 * @param status whether the membrane is active or not
	 * @param tend membrane t_end
	 */
	void setMembrane(double stiff, double r_zero, bool status, double tend);

	/** define all particles in the container
	 * @param plist the list to fill container
	 */
	void setList(std::list<Particle> plist);


	/** add one particle to the container
	 * @param p the particle to add
	 * @return a pointer to the particle
	*/
	Particle* addParticle(Particle& p);
	/** remove a particle from the container
	 * @param p the particle to remove
	 */
	void deleteParticle(Particle& p);
	/** get a list of all particle in the container
	 * @return list of all particles
	 */
	std::list<Particle>& getAll();
	/** get a list of pointer to all particle in the container
	 * @return list of all particles
	 */
	virtual std::list<Particle*> getAllPointer();
	/** get a list of all particle in the boundary
	 * @return list of all particles
	 */
	std::list<Particle>& getBoundary();
	/** get a list of all particle in the halo
	 * @return list of all particles
	 */
	std::list<Particle>& getHalo();
	/** delete all particles in halo
	 */
	void clearHalo();
	/** set the force to all particles according to gravity force
	 */
	void fReset();

	/** tests if the force up should be switched of
	 * @param time current simulation time
	 */
	void testForceEnd(double time);
	/** computes the forces in a membrane
	 * */
	void HarmonicPotential();

	/**applies the reflecting boundary condition
	 * @param p the particle near the border
	 * @param i the cell index of the current cell
	 */
	void ReflectingBoundary(Particle& p, int i);

	/** creates halo particle for particles on the boundary
	 */
	void PeriodicBoundary();

	/**
	 * calculate the force for all particles
	 */
	void calculateF();

	/** calculate the force for given cell
	 * @param x x-coord of cell
	 * @param y y-coord of cell
	 * @param z z-coord of cell
	 */
	void calculateForceForCell(int x, int y, int z);

	/**
	 * calculate the position for all particles
	 * using the Velocity-Stoermer-Verlet-Algorithm
	 */
	void calculateX();

	/**
	 * calculate the position for all particles
	 * using the Velocity-Stoermer-Verlet-Algorithm
	 */
	void calculateV();

private:
	/**
	 * Calculate the index of a cell
	 * @param x x-coord of cell
	 * @param y y-coord of cell
	 * @param z z-coord of cell
	 * @return index of cell
	 */
	int INDEX(int x, int y, int z);
	/**
	 * Calculate the x-coord of a cell
	 * @param i index of cell
	 * @return x-coord of cell
	 */
	int INDEX_X(int i);
	/**
	 * Calculate the y-coord of a cell
	 * @param i index of cell
	 * @return y-coord of cell
	 */
	int INDEX_Y(int i);
	/**
	 * Calculate the z-coord of a cell
	 * @param i index of cell
	 * @return z-coord of cell
	 */
	int INDEX_Z(int i);

	/**
	 * Calculate the x-coord of the cell which contains a particle with given x-coord
	 * @param x x-coord of a particle
	 * @return x-coord of cells which contain given particle
	 */
	int GET_X(double x);
	/**
	 * Calculate the y-coord of the cell which contains a particle with given y-coord
	 * @param y y-coord of a particle
	 * @return y-coord of cells which contain given particle
	 */
	int GET_Y(double y);
	/**
	 * Calculate the z-coord of the cell which contains a particle with given z-coord
	 * @param z z-coord of a particle
	 * @return z-coord of cells which contain given particle
	 */
	int GET_Z(double z);
};

inline int LinkedCellsContainer::INDEX(int x, int y, int z) { return (z)*NcellsXY + (y)*NcellsX + (x); }
inline int LinkedCellsContainer::INDEX_X(int i) { return i % NcellsX; }
inline int LinkedCellsContainer::INDEX_Y(int i) { return (i / NcellsX) % NcellsY; }
inline int LinkedCellsContainer::INDEX_Z(int i) { return i / NcellsXY; }

inline int LinkedCellsContainer::GET_X(double x) { return std::floor((x - lowerleft[0]) / HcellsX) + 1; }
inline int LinkedCellsContainer::GET_Y(double y) { return std::floor((y - lowerleft[1]) / HcellsY) + 1; }
inline int LinkedCellsContainer::GET_Z(double z) { return std::floor((z - lowerleft[2]) / HcellsZ) + 1; }


#endif /* SRC_LINKEDCELLSCONTAINER_H_ */
