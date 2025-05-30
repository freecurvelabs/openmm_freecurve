#---------------------------------------------------
# OpenMM ARROW Plugin
#
# Creates OpenMMArrow plugin library.
#
# Windows:
#   OpenMMArrow.dll
#   OpenMMArrow.lib
#   OpenMMArrow_static.lib
# Unix:
#   libOpenMMArrow.so
#   libOpenMMArrow_static.a
#----------------------------------------------------

#INCLUDE(Dart)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_ARROW_PLUGIN_SOURCE_SUBDIRS . openmmapi )

SET(OPENMM_ARROW_LIBRARY_NAME OpenMMArrow)
SET(SHARED_ARROW_TARGET ${OPENMM_ARROW_LIBRARY_NAME})
SET(STATIC_ARROW_TARGET ${OPENMM_ARROW_LIBRARY_NAME}_static)

# These are all the places to search for header files which are
# to be part of the API.
SET(API_ARROW_INCLUDE_DIRS) # start empty
FOREACH(subdir ${OPENMM_ARROW_PLUGIN_SOURCE_SUBDIRS})
    # append
    SET(API_ARROW_INCLUDE_DIRS ${API_ARROW_INCLUDE_DIRS}
                                     ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include
                                     ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_ARROW_REL_INCLUDE_FILES)   # start these out empty
SET(API_ARROW_ABS_INCLUDE_FILES)

FOREACH(dir ${API_ARROW_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ARROW_ABS_INCLUDE_FILES ${API_ARROW_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_ARROW_REL_INCLUDE_FILES ${API_ARROW_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_ARROW_FILES) # empty
SET(SOURCE_ARROW_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_ARROW_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_ARROW_FILES         ${SOURCE_ARROW_FILES}         ${src_files})   #append
    SET(SOURCE_ARROW_INCLUDE_FILES ${SOURCE_ARROW_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# ----------------------------------------------------------------------------

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

ADD_LIBRARY(${SHARED_ARROW_TARGET} SHARED ${SOURCE_ARROW_FILES} ${SOURCE_ARROW_INCLUDE_FILES} ${API_ARROW_ABS_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_ARROW_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_ARROW_BUILDING_SHARED_LIBRARY -DLEPTON_BUILDING_SHARED_LIBRARY")

FILE(GLOB serialization_files  ${CMAKE_CURRENT_SOURCE_DIR}/serialization/src/*.cpp)
SET_SOURCE_FILES_PROPERTIES(${serialization_files} PROPERTIES COMPILE_FLAGS "-DOPENMM_ARROW_BUILDING_SHARED_LIBRARY -DTIXML_USE_STL")

IF(OPENMM_BUILD_STATIC_LIB)
  ADD_LIBRARY(${STATIC_ARROW_TARGET} STATIC ${SOURCE_ARROW_FILES} ${SOURCE_ARROW_INCLUDE_FILES} ${API_ARROW_ABS_INCLUDE_FILES})
  SET_TARGET_PROPERTIES(${STATIC_ARROW_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_ARROW_BUILDING_STATIC_LIBRARY -DOPENMM_USE_STATIC_LIBRARIES -DOPENMM_BUILDING_STATIC_LIBRARY -DLEPTON_USE_STATIC_LIBRARIES -DLEPTON_BUILDING_STATIC_LIBRARY")
ENDIF(OPENMM_BUILD_STATIC_LIB)

# ----------------------------------------------------------------------------

# On Linux need to link to libdl
IF(NOT APPLE)
    FIND_LIBRARY(DL_LIBRARY dl)
    IF(DL_LIBRARY)
      TARGET_LINK_LIBRARIES(${SHARED_ARROW_TARGET} ${DL_LIBRARY})
      IF( CREATE_SERIALIZABLE_OPENMM_ARROW )
          TARGET_LINK_LIBRARIES(${SHARED_ARROW_SERIALIZABLE_TARGET} ${DL_LIBRARY})
      ENDIF( CREATE_SERIALIZABLE_OPENMM_ARROW )
      IF(OPENMM_BUILD_STATIC_LIB)
        TARGET_LINK_LIBRARIES(${STATIC_ARROW_TARGET} ${DL_LIBRARY})
      ENDIF(OPENMM_BUILD_STATIC_LIB)
    ENDIF(DL_LIBRARY)
ENDIF()

SET( OpenMMLib OpenMM )

TARGET_LINK_LIBRARIES( ${SHARED_ARROW_TARGET} ${SHARED_TARGET} )
IF(OPENMM_BUILD_STATIC_LIB)
    TARGET_LINK_LIBRARIES( ${STATIC_ARROW_TARGET} ${STATIC_TARGET} )
ENDIF(OPENMM_BUILD_STATIC_LIB)

# Which hardware platforms to build

ADD_SUBDIRECTORY(platforms/reference)
ADD_SUBDIRECTORY(platforms/common)
ADD_SUBDIRECTORY(platforms/cuda)

# turned off ARROW OPENCL library building for now - set next line ON when OPENCL version is available 
IF(OPENMM_BUILD_OPENCL_LIB)
    SET(OPENMM_BUILD_ARROW_OPENCL_LIB OFF CACHE BOOL "Build Arrow implementation for OpenCL")
ELSE(OPENMM_BUILD_OPENCL_LIB)
    SET(OPENMM_BUILD_ARROW_OPENCL_LIB OFF CACHE BOOL "Build Arrow implementation for OpenCL")
ENDIF(OPENMM_BUILD_OPENCL_LIB)
IF(OPENMM_BUILD_ARROW_OPENCL_LIB)
    ADD_SUBDIRECTORY(platforms/opencl)
ENDIF(OPENMM_BUILD_ARROW_OPENCL_LIB)

# turned off ARROW CUDA library building for now - set next line ON when OPENCL version is available 
IF(OPENMM_BUILD_CUDA_LIB)
    SET(OPENMM_BUILD_ARROW_CUDA_LIB ON CACHE BOOL "Build Arrow implementation for CUDA")
ELSE(OPENMM_BUILD_CUDA_LIB)
    SET(OPENMM_BUILD_ARROW_CUDA_LIB OFF CACHE BOOL "Build Arrow implementation for CUDA")
ENDIF(OPENMM_BUILD_CUDA_LIB)
#IF(OPENMM_BUILD_ARROW_CUDA_LIB)
#    ADD_SUBDIRECTORY(platforms/cuda)
#ENDIF(OPENMM_BUILD_ARROW_CUDA_LIB)

INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_ARROW_TARGET})
IF(OPENMM_BUILD_STATIC_LIB)
  INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${STATIC_ARROW_TARGET})
ENDIF(OPENMM_BUILD_STATIC_LIB)
FILE(GLOB CORE_HEADERS     include/*.h          */include/*.h)
FILE(GLOB TOP_HEADERS      include/openmm/*.h          */include/openmm/*.h)
FILE(GLOB INTERNAL_HEADERS include/openmm/internal/*.h */include/openmm/internal/*.h )
INSTALL_FILES(/include                 FILES ${CORE_HEADERS})
INSTALL_FILES(/include/openmm          FILES ${TOP_HEADERS})
INSTALL_FILES(/include/openmm/internal FILES ${INTERNAL_HEADERS})

#
# Testing  
#

ENABLE_TESTING()

IF (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
ELSE (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH .)
ENDIF (EXECUTABLE_OUTPUT_PATH)

#INCLUDE(ApiDoxygen.cmake)

IF(BUILD_TESTING AND OPENMM_BUILD_SERIALIZATION_TESTS)
    ADD_SUBDIRECTORY(serialization/tests)
ENDIF(BUILD_TESTING AND OPENMM_BUILD_SERIALIZATION_TESTS)
