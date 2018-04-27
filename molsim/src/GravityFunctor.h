/*
 * GravityFunctor.h
 *
 *  Created on: Apr 6, 2018
 *      Author: raffi
 */

#ifndef SRC_GRAVITYFUNCTOR_H_
#define SRC_GRAVITYFUNCTOR_H_

#include "../mdutils.h"
#include "ParticleMS.h"

#include <array>
#include <assert.h>
#include <cmath>


using namespace autopas;

template <class Particle, class ParticleCell>
class GravityFunctor : public Functor<Particle, ParticleCell>  {
public:

	/**@brief Constructor of the functor
	 * @param cut-off radius
	 */
	GravityFunctor<Particle, ParticleCell>(double rcut){r_cutoff = rcut; r_cutoff_square = rcut*rcut;}

	  /**
	   * @brief Functor for arrays of structures (AoS).
	   *
	   * This functor should calculate the forces or any other pair-wise interaction
	   * between two particles.
	   * This should include a cutoff check if needed!
	   */
	  void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {

		  //get the distance between both particles
		  std::array<double, 3> r = arrayMath::sub(j.getR(), i.getR());
			double rabs = sqrt(arrayMath::dot(r,r));

			assert(!(rabs==0.0));
			double rabs3 = rabs * rabs * rabs;

			//compute the new force
			std::array<double, 3> addforce = arrayMath::mulScalar(r, i.getM()*j.getM() / rabs3);

			//add it to the particles
			i.addF(addforce);
			j.subF(addforce);
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

	  bool allowsNewton3() override { return true; }

	  bool allowsNonNewton3() override { return false; }
private:

	  double r_cutoff;
	  double r_cutoff_square;
};




#endif /* SRC_GRAVITYFUNCTOR_H_ */
