#---------------------------------------------------
# OpenMM Reference Arrow Force Library
#
# Creates OpenMMArrowReference library.
#
# Windows:
#   OpenMMArrowReference.dll
#   OpenMMArrowReference.lib
# Unix:
#   libOpenMMArrowReference.so
#----------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS . )
SET(ARBALEST_ROOT ${CMAKE_SOURCE_DIR}/../ARBALEST)
SET(FFTW_DIR  ${CMAKE_SOURCE_DIR}/../FFTW_LIBS_SHARED)

# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMARROWREFERENCE_LIBRARY_NAME OpenMMArrowReference)

SET(SHARED_TARGET ${OPENMMARROWREFERENCE_LIBRARY_NAME})

# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS) # start empty
FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_REL_INCLUDE_FILES)   # start these out empty
SET(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.c)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/reference/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/reference/src)

#Add Arbalest include directories
INCLUDE_DIRECTORIES(BEFORE ${ARBALEST_ROOT})
INCLUDE_DIRECTORIES(BEFORE ${ARBALEST_ROOT}/ARBALEST)

# Adding CUDA includes for Arblest libraries that require cuda.h
INCLUDE_DIRECTORIES(BEFORE ${CUDAToolkit_INCLUDE_DIRS})

# Create the library

INCLUDE_DIRECTORIES(${REFERENCE_INCLUDE_DIR})

ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${OPENMM_LIBRARY_NAME})
TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${SHARED_ARROW_TARGET} optimized ${SHARED_ARROW_TARGET})
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_BUILDING_SHARED_LIBRARY")

# link directories for Arbalest libraries 
target_link_directories(${SHARED_TARGET} PUBLIC ${ARBALEST_ROOT}/FFT/lib/x64)
IF(WIN32)
# WINDOWS
target_link_directories(${SHARED_TARGET} PUBLIC
    $<$<CONFIG:Debug>:${ARBALEST_ROOT}/x64/Bin/Debug>
    $<$<CONFIG:Release>:${ARBALEST_ROOT}/x64/Bin/Release>
)

find_package(CUDAToolkit REQUIRED)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/ARBALEST.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/CommonLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/ConfigurationFileLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/FileIOLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/ForceFieldLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/SimulationCoreLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/SimulationToolsLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/SystemLoadingLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/TopologyLib.lib)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/GpuLib.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_LIB_DIR}/SBDescLib.lib)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/FFT/lib/x64/libfftw3-3.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/FFT/lib/x64/libfftw3f-3.lib)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/FFT/lib/x64/libfftw3l-3.lib)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/TFLib/libtensorflow-cpu-windows-x86_64-2.7.0/lib/tensorflow.lib)

ELSE()
# Linux 
  
find_package(CUDAToolkit REQUIRED)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/ConfigurationFileLib/libConfigurationFileLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/ForceFieldLib/libForceFieldLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/SimulationCoreLib/libSimulationCoreLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/SimulationToolsLib/libSimulationToolsLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/SystemLoadingLib/libSystemLoadingLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/TopologyLib/libTopologyLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/SimulationCoreLib/libSimulationCoreLibGpu.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/FileIOLib/libFileIOLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/CommonLib/libCommonLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/SimulationCoreLib/libSimulationCoreLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/ConfigurationFileLib/libConfigurationFileLib.a)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/GpuLib/libGpuLib.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/build/SBDescLib/libSBDescLibShared.so)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} -fopenmp)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} /srv/data/permanent/arbalest/TFLIB/libtensorflow-cpu-linux-x86_64-2.7.0/lib/libtensorflow.so)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${FFTW_DIR}/libfftw3.so)
#TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/FFT/lib/gcc-7.4/libfftw3f.a)
#TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/FFT/lib/gcc-7.4/libfftw3l.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${FFTW_DIR}/libfftw3_omp.so)
#TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/FFT/lib/gcc-7.4/libfftw3f_omp.a)
#TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${ARBALEST_ROOT}/FFT/lib/gcc-7.4/libfftw3l_omp.a)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lfftw3)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lfftw3f)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lfftw3l)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lfftw3_omp)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lfftw3f_omp)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lfftw3l_omp)

TARGET_LINK_LIBRARIES(${SHARED_TARGET}  )
TARGET_LINK_LIBRARIES(${SHARED_TARGET} CUDA::cudart CUDA::cufft CUDA::cuda_driver)

TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lboost_program_options )
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lboost_filesystem )
TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lboost_regex )

TARGET_LINK_LIBRARIES(${SHARED_TARGET} -lopenblas )

set(COORDSPEC "")
set(DOUBSPEC "")
set(PRECISION "-D__FLOATDOUBLEFIXED__")
#set(PRECISION "-D__DOUBLE__")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-128 -finline-functions -Wno-inline -pthread -fPIC -std=c++0x ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARBALEST_CUDA -DARBALEST_NNADJ -DUSE_FFTIMPL -DUSE_OPENMP -D__NN_DOUBLE__ -D__DOUBLE__ -march=native -MD -MT ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARBALEST_CUDA -DARBALEST_NNADJ -DUSE_FFTIMPL -DUSE_OPENMP  -march=native -MD -MT ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-128 -finline-functions -Wno-inline -pthread -fPIC -std=c++11 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PRECISION} ${COORDSPEC} ${DOUBSPEC} ")

ENDIF()

INSTALL(TARGETS ${SHARED_TARGET} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/plugins)

IF(BUILD_TESTING AND OPENMM_BUILD_REFERENCE_TESTS)
    SUBDIRS (tests)
ENDIF(BUILD_TESTING AND OPENMM_BUILD_REFERENCE_TESTS)
