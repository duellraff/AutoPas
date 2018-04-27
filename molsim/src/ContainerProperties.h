#ifndef SRC_CONTAINERPROPETIES_H_
#define SRC_CONTAINERPROPETIES_H_

#include "autopasIncludes.h"

#include <assert.h>
#include <omp.h>
#include "MaxwellBoltzmannDistribution.h"
#include <list>
#include "utils/Vector.h"
#include <vector>

using namespace autopas;

template<class autopasClass, class Functor>
class ContProperties{
private:

	double delta_t;
	Functor* fm;

	double BrownianFactor;

	std::array<int, 6> boundary;

	autopasClass* AutoPasMS;

	//membrane data
	double stiffness = 0;
	double r0 = 0;
	bool membrane_active = false;
	double t_end_force = 0;

public:
	/** Constructor
	 * @param pointer to the associated container
	 * @param f Functor function to calculate the force between two particles
	 * @param boundary conditions on the sides of the container
	 * @param timestep step between each step
	 * @param Brownian the brownian Factor
	 */

	ContProperties(autopasClass* autopasMS, Functor* f,
			std::array<int, 6> bound, double timestep=0.014, double Brownian=0.1){
		AutoPasMS = autopasMS; fm = f; delta_t = timestep; BrownianFactor = Brownian; boundary = bound;}


	/** sets the membrane parameters
	 * @param stiff membrane stiffness
	 * @param r_zero average bond lenght
	 * @param status whether the membrane is active or not
	 * @param tend force t_end
	 */
	void setMembrane(double stiff, double r_zero, bool status, double tend){
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

	Functor* getFunctor(){return fm;}

	double getBF(){return BrownianFactor;}

	double getStiffness(){return stiffness;}

	double getTEndForce(){return t_end_force;}

	bool getMembraneStatus(){return membrane_active;}

	autopasClass* getPointerContainer(){return AutoPasMS;}

	void setGravity(std::array<double, 3> g){
		for (auto it = AutoPasMS->begin(); it.isValid(); ++it){
			it->getGrav() = g;
			it->getConstF() = arrayMath::add(it->getConstF(), arrayMath::mulScalar(g, it->getM()));
		}
	}

	long unsigned getNumParticles(){
		long unsigned count = 0;
		//typename ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>::iterator
		for (auto it = AutoPasMS->begin(); it.isValid(); ++it){
			++count;
		}
		return count;
	}

	};


#endif /* SRC_CONTAINERPROPETIES_H_ */
