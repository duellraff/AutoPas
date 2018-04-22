/* MembraneFunctor.h
 *
 *  Created on: Apr 6, 2018
 *      Author: raffi
 */

#ifndef SRC_MEMBRANEFUNCTOR_H_
#define SRC_MEMBRANEFUNCTOR_H_

#include "../mdutils.h"
#include "ParticleMS.h"

#include <array>
#include <assert.h>
#include <cmath>


using namespace autopas;

template <class Particle, class ParticleCell>
class MembraneFunctor : public Functor<Particle, ParticleCell>  {
public:
	MembraneFunctor<Particle, ParticleCell>(){}
	  /**
	   * @brief Functor for arrays of structures (AoS).
	   *
	   * This functor should calculate the membrane forces
	   * between a molecule and its direct and diagonal neighbours
	   * This should include a cutoff check if needed!
	   *
	   * @param i the considered Particle
	   * @param k stiffness of the membrane
	   * @param r_zero unstressed distance between molecules
	   */
	  void AoSFunctor(Particle &i, double k, double r_zero) override;

	  /**TODO
	   * @brief Functor for structure of arrays (SoA)
	   *
	   * This functor should calculate the forces or any other pair-wise interaction
	   * between all particles in soa.
	   * This should include a cutoff check if needed!
	   *
	   * @param soa Structure of arrays
	   */
	  void SoAFunctor(SoA &soa) override {}

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
	  void SoAFunctor(SoA &soa1, SoA &soa2) override {}

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
	};


#endif /* SRC_MEMBRANEJONESFUNCTOR_H_ */
