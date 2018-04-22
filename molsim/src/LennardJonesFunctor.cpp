/*
 * GravityFunctor.cpp
 *
 *  Created on: Apr 6, 2018
 *      Author: raffi
 */

#include "LennardJonesFunctor.h"

using namespace autopas;

/*template<typename Particle, typename ParticleCell>
void LennardJonesFunctor<Particle, ParticleCell>::AoSFunctor(Particle &i, Particle &j) {
	std::array<double, 3> r = arrayMath::sub(i.getR(), j.getR());

	double rabs2 = arrayMath::dot(r,r);
	assert(!(rabs2==0.0));

	if (rabs2<r_cutoff_square){
		double sigma = (i.getSigma() + j.getSigma()) / 2;
		double epsilon = ParticleType::getEpsilon(i.getTypeID(), j.getTypeID());

		double fraction = sigma * sigma / rabs2;
		double fraction3 = fraction * fraction * fraction;
		double fraction6 = fraction3 * fraction3;

		std::array<double, 3> addforce = arrayMath::mulScalar(r, 24 * epsilon / rabs2 * (fraction3 - 2 * fraction6));

		i.addF(addforce);
		j.subF(addforce);
	}

}*/



