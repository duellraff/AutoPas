/*
 * GravityFunctor.cpp
 *
 *  Created on: Apr 6, 2018
 *      Author: raffi
 */

#include "GravityFunctor.h"

using namespace autopas;

/*template<typename Particle, typename ParticleCell>
void GravityFunctor<Particle, ParticleCell>::AoSFunctor(Particle &i, Particle &j) {
	std::array<double, 3> r = arrayMath::sub(i.getR(), j.getR());
	double rabs = sqrt(arrayMath::dot(r,r));

	if (rabs < r_cutoff){
		assert(!(rabs==0.0));
		double rabs3 = rabs * rabs * rabs;

		std::array<double, 3> addforce = arrayMath::mulScalar(r, i.getM()*j.getM() / rabs3);

		i.addF(addforce);
		j.subF(addforce);
	}

}*/





