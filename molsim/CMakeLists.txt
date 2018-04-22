project(MolSim)

cmake_minimum_required(VERSION 3.3)

#CMAKE policies, etc.
CMAKE_POLICY(SET CMP0057 NEW)
set(CMAKE_CXX_STANDARD 11)
# cmake < 3.6 ignores "set(CMAKE_CXX_STANDARD 11)" when using icpc
if (CMAKE_VERSION VERSION_LESS "3.6" AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

set(LINKER_FLAGS "-lxerces-c  -lcppunit -lapr-1  -laprutil-1")

ADD_EXECUTABLE(MolSim
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(MolSim
	${LINKER_FLAGS})

TARGET_LINK_LIBRARIES(MolSim
	autopas)
