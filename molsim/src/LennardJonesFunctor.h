/* GravityFunctor.h
 *
 *  Created on: Apr 6, 2018
 *      Author: raffi
 */

#ifndef SRC_LENNARDJONESFUNCTOR_H_
#define SRC_LENNARDJONESFUNCTOR_H_

#include "../mdutils.h"
#include "ParticleMS.h"

#include <array>
#include <assert.h>
#include <cmath>


using namespace autopas;

template <class Particle, class ParticleCell>
class LennardJonesFunctor : public Functor<Particle, ParticleCell>  {
public:
	/**@brief Constructor of the functor
	 * @param cut-off radius
	 */
	LennardJonesFunctor<Particle, ParticleCell>(double rcut){r_cutoff = rcut; r_cutoff_square = rcut*rcut;}

	  /**
	   * @brief Functor for arrays of structures (AoS).
	   *
	   * This functor should calculate the forces or any other pair-wise interaction
	   * between two particles.
	   * This should include a cutoff check if needed!
	   */
	  void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override{
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
				j.addF(addforce);
				i.subF(addforce);
			}

		}

	  /**TODO
	   * @brief Functor for structure of arrays (SoA)
	   *
	   * This functor should calculate the forces or any other pair-wise interaction
	   * between all particles in soa.
	   * This should include a cutoff check if needed!
	   *
	   * @param soa Structure of arrays
	   */
	  void SoAFunctor(SoA &soa, bool newton3 = true) override {}

	  /**TODO
	   * @brief Functor for structure of arrays (SoA)
	   *
	   * This functor should calculate the forces or any other pair-wise interaction
	   * between all particles of soa1 and soa2.
	   * This should include a cutoff check if needed!
	   *
	   * @param soa1 First structure of arrays.
	   * @param soa2 Second structure of arrays.
	   */
	  void SoAFunctor(SoA &soa1, SoA &soa2, bool newton3 = true) override {}

	  /**TODO
	   * @brief Copies the AoS data of the given cell in the given soa.
	   *
	   * @param cell Cell from where the data is loaded.
	   * @param soa  Structure of arrays where the data is copied to.
	   */
	  void SoALoader(ParticleCell &cell, SoA *soa) override {}

	  /**TODO
	   * @brief Copies the data stored in the soa in the cell.
	   *
	   * @param cell Cell where the data should be stored.
	   * @param soa  Structure of arrays from where the data is loaded.
	   */
	  void SoAExtractor(ParticleCell *cell, SoA *soa) {}

private:
	  double r_cutoff;
	  double r_cutoff_square;
};





#endif /* SRC_LENNARDJONESFUNCTOR_H_ */
