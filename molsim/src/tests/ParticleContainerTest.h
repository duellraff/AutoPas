/*
 * ParticleContainerTest.h
 *
 *  Created on: Nov 13, 2017
 *      Author: raffi
 */

#ifndef SRC_PARTICLECONTAINERTEST_H_
#define SRC_PARTICLECONTAINERTEST_H_


#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>
#include <../mdutils.h>
#include <ParticleMS.h>


/**
 * @brief TestSuite for ParticleContainer
 */
class ParticleContainerTest : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(ParticleContainerTest);
	CPPUNIT_TEST(addParticleTest);
	CPPUNIT_TEST(deleteParticleTest);
	CPPUNIT_TEST(fclearTest);
	CPPUNIT_TEST(calculateFTest);
	CPPUNIT_TEST(GenerateCuboidTest);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp();
	/**
	 * Tests if particles are added,
	 * but copies are ignored
	 */
	void addParticleTest();
	/**
	 * Tests if particles are deleted if in container
	 */
	void deleteParticleTest();
	/**
	 * Tests if forces are set to zero
	 * and the old-force variable is set accordingly
	 */
	void fclearTest();
	/**
	 * Test if the forces are calculated correctly
	 */
	void calculateFTest();
	/**
	 * Test if the right amount of particles are generated
	 */
	void GenerateCuboidTest();


	ParticleContainerTest();
	virtual ~ParticleContainerTest();




};
#endif /* SRC_PARTICLECONTAINERTEST_H_ */
