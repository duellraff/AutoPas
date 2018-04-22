/*
 * ParticleContainerTest.h
 *
 *  Created on: Nov 13, 2017
 *      Author: raffi
 */

#ifndef SRC_LINKEDCONTAINERTEST_H_
#define SRC_LINKEDCONTAINERTEST_H_


#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../mdutils.h"
#include "ParticleMS.h"


/**
 * @brief TestSuite for LinkedContainer
 */
class LinkedContainerTest : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(LinkedContainerTest);
	CPPUNIT_TEST(addParticleTest);
	CPPUNIT_TEST(deleteParticleTest);
	CPPUNIT_TEST(outOfDomainTest);
	CPPUNIT_TEST(boundaryForceTest);
	CPPUNIT_TEST(boundaryTest);
	CPPUNIT_TEST(calculateFTest);
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
	 * Tests if particles are deleted if out of domain
	 */
	void outOfDomainTest();

	/**
	 * Tests if particles are deleted if out of domain
	 */
	void boundaryForceTest();

	/**
	 * Tests if particles are recognized in the boundary
	 */
	void boundaryTest();
	/**
	 * Test if the forces are calculated correctly
	 */
	void calculateFTest();


	LinkedContainerTest();
	virtual ~LinkedContainerTest();

};
#endif /* SRC_LINKEDCONTAINERTEST_H_ */
