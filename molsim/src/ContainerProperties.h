#ifndef SRC_CONTAINERPROPETIES_H_
#define SRC_CONTAINERPROPETIES_H_

#include "../mdutils.h"

#include <assert.h>
#include <omp.h>
#include "MaxwellBoltzmannDistribution.h"
#include <list>
#include "utils/Vector.h"
#include <vector>

using namespace autopas;

class ContProperties{
private:

	double delta_t;
	Functor<ParticleMS, FullParticleCell<ParticleMS>>* fm;

	double BrownianFactor;

	std::array<int, 6> boundary;

	ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>* container;

	//membrane data
	double stiffness = 0;
	double r0 = 0;
	bool membrane_active = false;
	double t_end_force = 0;

public:
	/** Constructor
	 * @param pointer to the assoiated container
	 * @param f Functor function to calculate the force between two particles
	 * @param boundary conditions on the sides of the container
	 * @param timestep step between each step
	 * @param Brownian the brownian Factor
	 */
	ContProperties(ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>* cont,
			Functor<ParticleMS, FullParticleCell<ParticleMS>>* f,
			std::array<int, 6> bound, double timestep=0.014, double Brownian=0.1){
		container = cont; fm = f; delta_t = timestep; BrownianFactor = Brownian; boundary = bound;}

	//virtual ~ContProperties();

	/** sets the membrane parameters
	 * @param stiff membrane stiffness
	 * @param r_zero average bond lenght
	 * @param status whether the membrane is active or not
	 * @param tend force t_end
	 */
	virtual void setMembrane(double stiff, double r_zero, bool status, double tend){
		stiffness = stiff;
		r0 = r_zero;
		membrane_active = status;
		t_end_force = tend;
	}

	double getDeltaT(){return delta_t;}

	/** set the timestep
	 * @param timestep time between steps
	 */
	void setDeltaT(double timestep){delta_t = timestep;}

	Functor<ParticleMS, FullParticleCell<ParticleMS>>* getFunctor(){return fm;}

	double getBF(){return BrownianFactor;}

	double getStiffness(){return stiffness;}

	double getTEndForce(){return t_end_force;}

	bool getMembraneStatus(){return membrane_active;}

	ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>* getPointerContainer(){return container;}
//typename ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>::iterator
	void setGravity(std::array<double, 3> g){
		for (auto it = container->begin(); it.isValid(); ++it){
			it->getGrav() = g;
			it->getConstF() = arrayMath::add(it->getConstF(), arrayMath::mulScalar(g, it->getM()));
		}
	}

	long unsigned getNumParticles(){
		long unsigned count = 0;
		//typename ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>::iterator
		for (auto it = container->begin(); it.isValid(); ++it){
			++count;
		}
		return count;
	}

	};


#endif /* SRC_CONTAINERPROPETIES_H_ */
