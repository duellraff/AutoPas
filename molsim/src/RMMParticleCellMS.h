/*
 * RMMParticleCellMS.h
 *
 *  Created on: Apr 25, 2018
 *      Author: raffi
 */

#ifndef SRC_RMMPARTICLECELLMS_H_
#define SRC_RMMPARTICLECELLMS_H_

#include "autopasIncludes.h"
#include "ParticleMS.h"

using namespace autopas;

template <class Particle, class Iterator>
class RMMParticleCellMS : public ParticleCell<Particle, Iterator,
			RMMParticleCellMS<Particle, Iterator>> {
public:
	RMMParticleCellMS(){
	    _particleSoABuffer.initArrays(
	        {Particle::AttributeNames::id, Particle::AttributeNames::posX,
	    	 Particle::AttributeNames::posY, Particle::AttributeNames::posZ,
			 Particle::AttributeNames::velX, Particle::AttributeNames::velY,
			 Particle::AttributeNames::velZ, Particle::AttributeNames::forceX,
	         Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ,
			 Particle::AttributeNames::oldforceX, Particle::AttributeNames::oldforceY,
			 Particle::AttributeNames::oldforceZ, Particle::AttributeNames::constforceX,
	         Particle::AttributeNames::constforceY, Particle::AttributeNames::constforceZ,
			 Particle::AttributeNames::ggravX, Particle::AttributeNames::ggravY,
			 Particle::AttributeNames::ggravZ, Particle::AttributeNames::directNeigh1,
			 Particle::AttributeNames::directNeigh2, Particle::AttributeNames::directNeigh3,
			 Particle::AttributeNames::directNeigh4, Particle::AttributeNames::diagonalNeigh1,
			 Particle::AttributeNames::diagonalNeigh2, Particle::AttributeNames::diagonalNeigh3,
			 Particle::AttributeNames::diagonalNeigh4, Particle::AttributeNames::typeID
	        });

	}

	void addParticle(Particle &m) override {
		_particleSoABuffer.push(Particle::AttributeNames::id, static_cast<double>(m.getID()));
		_particleSoABuffer.push(Particle::AttributeNames::posX, m.getR()[0]);
	    _particleSoABuffer.push(Particle::AttributeNames::posY, m.getR()[1]);
	    _particleSoABuffer.push(Particle::AttributeNames::posZ, m.getR()[2]);
		_particleSoABuffer.push(Particle::AttributeNames::velX, m.getV()[0]);
	    _particleSoABuffer.push(Particle::AttributeNames::velY, m.getV()[1]);
	    _particleSoABuffer.push(Particle::AttributeNames::velZ, m.getV()[2]);
	    _particleSoABuffer.push(Particle::AttributeNames::forceX, m.getF()[0]);
	    _particleSoABuffer.push(Particle::AttributeNames::forceY, m.getF()[1]);
	    _particleSoABuffer.push(Particle::AttributeNames::forceZ, m.getF()[2]);
	    _particleSoABuffer.push(Particle::AttributeNames::oldforceX, m.getOldF()[0]);
	    _particleSoABuffer.push(Particle::AttributeNames::oldforceY, m.getOldF()[1]);
	    _particleSoABuffer.push(Particle::AttributeNames::oldforceZ, m.getOldF()[2]);
	    _particleSoABuffer.push(Particle::AttributeNames::constforceX, m.getConstF()[0]);
	    _particleSoABuffer.push(Particle::AttributeNames::constforceY, m.getConstF()[1]);
	    _particleSoABuffer.push(Particle::AttributeNames::constforceZ, m.getConstF()[2]);
	    _particleSoABuffer.push(Particle::AttributeNames::ggravX, m.getGrav()[0]);
	    _particleSoABuffer.push(Particle::AttributeNames::ggravY, m.getGrav()[1]);
	    _particleSoABuffer.push(Particle::AttributeNames::ggravZ, m.getGrav()[2]);
/*	    _particleSoABuffer.push(Particle::AttributeNames::directNeigh1, static_cast<double>(m.getDirectNeighbours()[0]->getID()));
	    _particleSoABuffer.push(Particle::AttributeNames::directNeigh2, static_cast<double>(m.getDirectNeighbours()[1]->getID()));
	    _particleSoABuffer.push(Particle::AttributeNames::directNeigh3, static_cast<double>(m.getDirectNeighbours()[2]->getID()));
	    _particleSoABuffer.push(Particle::AttributeNames::directNeigh4, static_cast<double>(m.getDirectNeighbours()[3]->getID()));
	    _particleSoABuffer.push(Particle::AttributeNames::diagonalNeigh1, static_cast<double>(m.getDiagonalNeighbours()[0]->getID()));
	    _particleSoABuffer.push(Particle::AttributeNames::diagonalNeigh2, static_cast<double>(m.getDiagonalNeighbours()[1]->getID()));
	    _particleSoABuffer.push(Particle::AttributeNames::diagonalNeigh3, static_cast<double>(m.getDiagonalNeighbours()[2]->getID()));
	    _particleSoABuffer.push(Particle::AttributeNames::diagonalNeigh4, static_cast<double>(m.getDiagonalNeighbours()[3]->getID()));
*/		_particleSoABuffer.push(Particle::AttributeNames::typeID, static_cast<double>(m.getTypeID()));
	}

	  unsigned long numParticles() const override {
	    return _particleSoABuffer.getNumParticles();
	  }
	  bool isNotEmpty() const override { return numParticles() > 0; }

	  void clear() override { _particleSoABuffer.clear(); }

	  void deleteByIndex(int index) override {
	    assert(index >= 0 and index < numParticles());
	    assert(numParticles() > 0);
	    if (index < numParticles() - 1) {
	      _particleSoABuffer.swap(index, numParticles() - 1);
	    }
	    _particleSoABuffer.pop_back();
	  }

	  /**
	   * the soa buffer of the particle, all information is stored here.
	   */
	  SoA _particleSoABuffer;

	  /**
	   * iterator to iterate through ParticleCell
	   * If you need to explicitly store this iterator use
	   * typename RMMParticleCell<ParticleType>::iterator iter;
	   */
	  typedef Iterator iterator;

private:
	  void buildParticleFromSoA(size_t i, ParticleMS &rmm_or_not) {
	    rmm_or_not.getR() = _particleSoABuffer.read<3>(
	        {Particle::AttributeNames::posX, Particle::AttributeNames::posY,
	         Particle::AttributeNames::posZ},
	        i);
	    rmm_or_not.getV() =_particleSoABuffer.read<3>(
	        {Particle::AttributeNames::velX, Particle::AttributeNames::velY,
	         Particle::AttributeNames::velZ},
	        i);
	    rmm_or_not.getF() = _particleSoABuffer.read<3>(
	        {Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
	         Particle::AttributeNames::forceZ},
	        i);
	    rmm_or_not.getOldF() =_particleSoABuffer.read<3>(
	    	{Particle::AttributeNames::oldforceX, Particle::AttributeNames::oldforceY,
	         Particle::AttributeNames::oldforceZ},
	        i);
	    rmm_or_not.getConstF() =_particleSoABuffer.read<3>(
	    	{Particle::AttributeNames::constforceX, Particle::AttributeNames::constforceY,
	         Particle::AttributeNames::constforceZ},
	        i);
	    rmm_or_not.getGrav() =_particleSoABuffer.read<3>(
	    	{Particle::AttributeNames::ggravX, Particle::AttributeNames::ggravY,
	         Particle::AttributeNames::ggravZ},
	        i);
	    rmm_or_not.getTypeID() = static_cast<int>(_particleSoABuffer.read(ParticleMS::AttributeNames::typeID, i));
	  }

	  void writeParticleToSoA(size_t index, Particle &particle) {
	    _particleSoABuffer.write<3>(
	        {Particle::AttributeNames::posX, Particle::AttributeNames::posY,
	         Particle::AttributeNames::posZ},
	        index, particle.getR());
	    _particleSoABuffer.write<3>(
	        {Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
	         Particle::AttributeNames::forceZ},
	        index, particle.getF());
	    _particleSoABuffer.write<3>(
	        {Particle::AttributeNames::velX, Particle::AttributeNames::velY,
	         Particle::AttributeNames::velZ},
	        index, particle.getV());
	    _particleSoABuffer.write<3>(
	    	{Particle::AttributeNames::oldforceX, Particle::AttributeNames::oldforceY,
	         Particle::AttributeNames::oldforceZ},
	        index, particle.getOldF());
	    _particleSoABuffer.write<3>(
	    	{Particle::AttributeNames::constforceX, Particle::AttributeNames::constforceY,
	         Particle::AttributeNames::constforceZ},
	        index, particle.getConstF());
	    _particleSoABuffer.write<3>(
	    	{Particle::AttributeNames::ggravX, Particle::AttributeNames::ggravY,
	         Particle::AttributeNames::ggravZ},
	        index, particle.getGrav());
	    _particleSoABuffer.write<1>({ParticleMS::AttributeNames::typeID}, index, {static_cast<double>(particle.getTypeID())});
	  }


	  template <class ParticleType>
	  friend class RMMParticleCellMSIterator;

};

/**
 * SingleCellIterator for the RMMParticleCell
 * @tparam Particle
 */
template <class Particle>
class RMMParticleCellMSIterator {
 public:
  /**
   * default constructor of SingleCellIterator
   */
  RMMParticleCellMSIterator() : _cell(nullptr), _index(0), _deleted(false) {}

  /**
   * constructor of SingleCellIterator
   * @param cell_arg pointer to the cell of particles
   * @param ind index of the first particle
   */
  RMMParticleCellMSIterator(
      RMMParticleCellMS<Particle, RMMParticleCellMSIterator<Particle>> *cell_arg,
      int ind = 0)
      : _cell(cell_arg), _index(ind), _deleted(false) {}

  //  SingleCellIterator(const SingleCellIterator &cellIterator) {
  //    _cell = cellIterator._cell;
  //    _AoSReservoir = cellIterator._AoSReservoir;
  //    _index = cellIterator._index;
  //  }

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  Particle &operator*() {
    // Particle * ptr = nullptr;
    // ptr = const_cast<Particle *>(& _AoSReservoir);
    Particle& p = _AoSReservoir;
    _cell->buildParticleFromSoA(_index, p);
    //_cell->particleAt(_index, ptr);
    return p;
  }

  /**
   * access particle using "iterator->"
   *
   * this is the member of pointer operator
   * @return current particle
   */
  Particle *operator->() { return &(this->operator*()); }

  /**
   * equality operator.
   * if both iterators are invalid or if they point to the same particle, this
   * returns true
   * @param rhs
   * @return
   */
  bool operator==(const RMMParticleCellMSIterator &rhs) const {
    return (not this->isValid() and not rhs.isValid()) or
           (_cell == rhs._cell && _index == rhs._index);
  }

  /**
   * inequality operator
   * descrition see operator==
   * @param rhs
   * @return
   */
  bool operator!=(const RMMParticleCellMSIterator &rhs) const {
    return !(rhs == *this);
  }

  /**
   * increment operator to get the next particle
   * @return the next particle, usually ignored
   */
  RMMParticleCellMSIterator &operator++() {
    if (not _deleted) {
		AutoPasLogger->info(_AoSReservoir.getR()[0]);
      _cell->writeParticleToSoA(_index, _AoSReservoir);
      ++_index;
    }
    _deleted = false;
    return *this;
  }

  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  bool isValid() const {
    return _cell != nullptr and _index < _cell->numParticles();
  }

  /**
   * Get the index of the particle in the cell
   * @return index of the current particle
   */
  int getIndex() const { return _index; }

  /**
   * Deletes the current particle
   */
  void deleteCurrentParticle() {
    _cell->deleteByIndex(_index);
    _deleted = true;
  }

 private:
  RMMParticleCellMS<Particle, RMMParticleCellMSIterator<Particle>> *_cell;
  Particle _AoSReservoir;
  int _index;
  bool _deleted;
};

// provide a simpler template for RMMParticleCell, i.e.
// RMMParticleCell<Particle>
template <class Particle>
using RMMParticleCell =
    RMMParticleCellMS<Particle, RMMParticleCellMSIterator<Particle>>;

#endif /* AUTOPAS_SRC_CELLS_RMMPARTICLECELLMS_H_ */
