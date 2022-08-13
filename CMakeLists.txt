CMAKE_MINIMUM_REQUIRED ( VERSION 2.8.12 )

SET ( LMP_PATH /path/to/lammps )
PROJECT ( "LAMMPS_NEB" )
SET ( PROJECT_VERSION_MAJOR 1 )
SET ( PROJECT_VERSION_MINOR 0 )

SET ( CMAKE_VERBOSE_MAKEFILE true )
SET ( SRC_FILES 
        main.c 
        calculator.c
        config.c
        input.c
        neb.c
        utils.c )
SET ( CMAKE_C_COMPILER "mpiicc" )

ADD_COMPILE_OPTIONS ( -O3 -g -std=c11 -Wall -fPIC )
LINK_LIBRARIES( lammps m )
LINK_DIRECTORIES( ${LMP_PATH}/src )
INCLUDE_DIRECTORIES ( ${LMP_PATH}/src )

ADD_EXECUTABLE (LAMMPS_NEB ${SRC_FILES})
