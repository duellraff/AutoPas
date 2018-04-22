
#include "MembraneFunctor.h"

#define Sqrt2 1.41421356

using namespace autopas;


template<typename Particle, typename ParticleCell>
void MembraneFunctor<Particle, ParticleCell>::AoSFunctor(Particle &p, double k, double r_zero) {
	//initialize membrane force and direction vector
	std::array<double,3> f({0.0, 0.0, 0.0});
	std::array<double,3> dir({0.0, 0.0, 0.0});


	// apply potential to direct neighbours
	for(int n=0; n<2; n++){
		if(p.getDirectNeighbours()[n] != NULL){
			Particle& pother = *p.getDirectNeighbours()[n];

			//compute direct potential
			dir = arrayMath::sub(p.getR(), pother.getR());
			double norm2 = sqrt(arrayMath::dot(dir, dir));
			f = arrayMath::mulScalar(dir, k*(1 - r_zero / norm2));

			//add the new force
			p.subF(f);
			pother.addF(f);
		}
	}

	// apply potential to diagonal neighbours
	for(int n=0; n<2; n++){
		if(p.getDiagonalNeighbours()[n] != NULL){
			Particle& pother = *p.getDiagonalNeighbours()[n];

			//compute diagonal potential
			dir = arrayMath::sub(p.getR(), pother.getR());
			double norm2 = sqrt(arrayMath::dot(dir, dir));
			f = arrayMath::mulScalar(dir, k*(1 - Sqrt2 * r_zero / norm2));

			//add the new force
			p.subF(f);
			pother.addF(f);
		}
	}
}




