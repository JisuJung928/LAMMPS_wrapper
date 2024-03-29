CMAKE_MINIMUM_REQUIRED ( VERSION 3.13 )
PROJECT ( LAMMPS_wrapper )
ENABLE_LANGUAGE ( C )

SET ( LAMMPS_SRC_DIRECTORY "" CACHE PATH "Path for lammps source" )
IF ( NOT LAMMPS_SRC_DIRECTORY STREQUAL "" AND EXISTS ${LAMMPS_SRC_DIRECTORY}/cmake/CMakeLists.txt )
  OPTION( BUILD_LIB "Build LAMMPS library" ON )
  ADD_SUBDIRECTORY ( ${LAMMPS_SRC_DIRECTORY}/cmake lammps )
ELSE()
  FIND_PACKAGE( LAMMPS REQUIRED )
ENDIF()

SET ( CMAKE_VERBOSE_MAKEFILE true )
SET ( SRC_FILES 
        main.c 
        calculator.c
        config.c
        input.c
        target.c
        utils.c )
SET ( CMAKE_C_COMPILER "mpiicc" )
ADD_COMPILE_OPTIONS ( -O3 -g -Wall )
LINK_LIBRARIES ( m )
ADD_EXECUTABLE ( LAMMPS_wrapper ${SRC_FILES} )

TARGET_LINK_LIBRARIES ( LAMMPS_wrapper LAMMPS::lammps )
