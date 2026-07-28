# CMakeFiles.txt for the EBTKS library
#
# Vladimir S. FONOV - vladimir.fonov@gmail.com

CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
SET(EBTKS_PACKAGE_VERSION_MAJOR 1)
SET(EBTKS_PACKAGE_VERSION_MINOR 6)
SET(EBTKS_PACKAGE_VERSION_PATCH 40)

# WARNING: compiling library as shared library produces wrong code for some reason
# Untill this problem is solved the code should be compiled as static

SET(EBTKS_BUILD_SHARED_LIBS OFF)
# Packaging defines
IF(NOT MINC_TOOLKIT_BUILD)
  PROJECT(EBTKS)
  SET(CPACK_GENERATOR TGZ)
  SET(CPACK_PACKAGE_VERSION_MAJOR ${EBTKS_PACKAGE_VERSION_MAJOR})
  SET(CPACK_PACKAGE_VERSION_MINOR ${EBTKS_PACKAGE_VERSION_MINOR})
  SET(CPACK_PACKAGE_VERSION_PATCH ${EBTKS_PACKAGE_VERSION_PATCH})
  #OPTION(EBTKS_BUILD_SHARED_LIBS "Build shared libraries." OFF)
  INCLUDE(CPack)
ENDIF(NOT MINC_TOOLKIT_BUILD)

ENABLE_TESTING()
INCLUDE(CTest)

SET(PACKAGE "EBTKS")
SET(PACKAGE_BUGREPORT "a.janke@gmail.com")
SET(PACKAGE_NAME "EBTKS")
SET(PACKAGE_VERSION "${EBTKS_PACKAGE_VERSION_MAJOR}.${EBTKS_PACKAGE_VERSION_MINOR}.${EBTKS_PACKAGE_VERSION_PATCH}")
SET(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_C_STANDARD 99)


IF(NOT MINC_TOOLKIT_BUILD)
	FIND_PACKAGE( LIBMINC REQUIRED )
ENDIF(NOT MINC_TOOLKIT_BUILD)

INCLUDE( ${LIBMINC_USE_FILE} )

ADD_DEFINITIONS(-DHAVE_CONFIG_H)
ADD_DEFINITIONS(-DUSE_COMPMAT)
ADD_DEFINITIONS(-DUSE_FCOMPMAT)
ADD_DEFINITIONS(-DUSE_DBLMAT)

# add for building relocatable library
IF(UNIX)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF(UNIX)

# Disable floating-point contraction (fused multiply-add) for this library.
#
# N3's NU estimation is an iterative feedback loop (~31 iterations) whose inner
# numerics -- the Bunch-Kaufman factorisation in dsysv_ and the FFT histogram
# deconvolution -- are exquisitely sensitive to rounding. Pivot selection in
# dsysv_ is a discrete magnitude comparison, so a 1-ULP difference flips the
# factorisation path and the trajectories diverge.
#
# ARM64 has FMA in the baseline ISA and compilers contract a*b+c into a fused
# op by default; baseline x86-64 has no FMA instruction and must round the
# intermediate. That alone makes the two platforms disagree: N3's
# nu_reference_1 regression test drifts ~8.5e-4 relative RMS on Apple Silicon
# against an x86-64-generated reference, over the test's 1e-4 tolerance.
#
# Verified on x86-64: rebuilding only this library and N3 with -mfma reproduces
# the failure (4.8e-3 RMS, 48x tolerance); adding -ffp-contract=off emits zero
# FMA instructions and restores a bit-identical result. Both N3 and EBTKS need
# the flag -- either one alone still contracts enough to diverge.
#
# -ffp-contract=off is GCC/Clang spelling; MSVC uses /fp: and rejects it.
IF(NOT MSVC)
  SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -ffp-contract=off")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-contract=off")
ENDIF()

# aliases
SET(VERSION "${PACKAGE_VERSION}")

# check for prereqs
INCLUDE(CheckFunctionExists)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckSymbolExists)

set(CMAKE_REQUIRED_LIBRARIES "m")

CHECK_SYMBOL_EXISTS(finite      "math.h" HAVE_FINITE)
CHECK_SYMBOL_EXISTS(isfinite    "math.h" HAVE_ISFINITE)

CHECK_FUNCTION_EXISTS(mkstemp  HAVE_MKSTEMP)

CHECK_INCLUDE_FILES(float.h     HAVE_FLOAT_H)
CHECK_INCLUDE_FILES(sys/dir.h   HAVE_SYS_DIR_H)
CHECK_INCLUDE_FILES(sys/ndir.h  HAVE_SYS_NDIR_H)
CHECK_INCLUDE_FILES(sys/stat.h  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(sys/wait.h  HAVE_SYS_WAIT_H)
CHECK_INCLUDE_FILES(values.h    HAVE_VALUES_H)
CHECK_INCLUDE_FILES(unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(dirent.h    HAVE_DIRENT_H)
CHECK_INCLUDE_FILES(memory.h    HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(malloc.h    HAVE_MALLOC_H)
CHECK_INCLUDE_FILES(stdlib.h    HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(fcntl.h     HAVE_FCNTL_H)
CHECK_INCLUDE_FILES(vfork.h     HAVE_VFORK_H)
CHECK_INCLUDE_FILES(inttypes.h  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(string.h    HAVE_STRING_H)
CHECK_INCLUDE_FILES(strings.h   HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(dlfcn.h     HAVE_DLFCN_H)
CHECK_INCLUDE_FILES(complex     HAVE_COMPLEX)
CHECK_INCLUDE_FILES(complex.h   HAVE_COMPLEX_H)

# HACKY HACK
IF(APPLE)
  SET(HAVE_ISFINITE ON)
ENDIF(APPLE)


message("HAVE_ISFINITE=${HAVE_ISFINITE}")
# config files for build
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)

# others
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/epm-header.in ${CMAKE_CURRENT_BINARY_DIR}/epm-header)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include  ${CMAKE_CURRENT_SOURCE_DIR}/templates )

SET(EBTKS_HEADERS
        include/EBTKS/amoeba.h
        include/EBTKS/assert.h
        include/EBTKS/backProp.h
        include/EBTKS/Complex.h
        include/EBTKS/dcomplex.h
        include/EBTKS/fcomplex.h
        include/EBTKS/FileIO.h
        include/EBTKS/Histogram.h
        include/EBTKS/matlabSupport.h
        include/EBTKS/Minc.h
        include/EBTKS/MPoint.h
        include/EBTKS/MString.h
        include/EBTKS/MTypes.h
        include/EBTKS/OpTimer.h
        include/EBTKS/OrderedCltn.h
        include/EBTKS/Path.h
        include/EBTKS/Polynomial.h
        include/EBTKS/popen.h
        include/EBTKS/TrainingSet.h
        include/EBTKS/trivials.h

        templates/EBTKS/Array.h
        templates/EBTKS/CachedArray.h
        templates/EBTKS/Dictionary.h
        templates/EBTKS/Matrix3D.h
        templates/EBTKS/Matrix.h
        templates/EBTKS/MatrixSupport.h
        templates/EBTKS/MatrixTest.h
        templates/EBTKS/miscTemplateFunc.h
        templates/EBTKS/Pool.h
        templates/EBTKS/SimpleArray.h
        templates/EBTKS/Stack.h
        templates/EBTKS/ValueMap.h
)

# some variables
SET(EBTKS_LIB_SRCS
  ${EBTKS_HEADERS}
  version.cc
  src/FileIO.cc
  src/Histogram.cc
  src/MPoint.cc
  src/MString.cc
  src/OpTimer.cc
  src/OrderedCltn.cc
  src/Path.cc
  src/Polynomial.cc
  src/TrainingSet.cc
  src/amoeba.cc
  src/backProp.cc
  src/fcomplex.cc
  src/dcomplex.cc

  templates/Array.cc
  templates/CachedArray.cc
  templates/Dictionary.cc
  templates/Matrix.cc
  templates/MatrixSupport.cc
  templates/Pool.cc
  templates/SimpleArray.cc
  templates/ValueMap.cc
  templates/miscTemplateFunc.cc
  templates/MatrixConv.cc

  clapack/dcopy.c
  clapack/dgemm.c
  clapack/dgemv.c
  clapack/dger.c
  clapack/dlasyf.c
  clapack/dscal.c
  clapack/dswap.c
  clapack/dsyr.c
  clapack/dsysv.c
  clapack/dsytf2.c
  clapack/dsytrf.c
  clapack/dsytrs.c
  clapack/idamax.c
  clapack/ieeeck.c
  clapack/ilaenv.c
  clapack/lsame.c
  clapack/s_cmp.c
  clapack/s_copy.c
  clapack/xerbla.c
)

IF(EBTKS_BUILD_SHARED_LIBS)
  SET(LIBRARY_TYPE SHARED)
  SET(LIBRARY_INSTALL LIBRARY)
ELSE(EBTKS_BUILD_SHARED_LIBS)
  SET(LIBRARY_TYPE STATIC)
  SET(LIBRARY_INSTALL ARCHIVE)
ENDIF(EBTKS_BUILD_SHARED_LIBS)

# to fix a bug with diverging N3
#STRING(REGEX REPLACE "-O3" "-O2" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
#STRING(REGEX REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})

ADD_LIBRARY(EBTKS ${LIBRARY_TYPE} ${EBTKS_LIB_SRCS})
SET(EBTKS_LIBRARIES EBTKS)

INSTALL(FILES ${EBTKS_HEADERS}
   DESTINATION include/EBTKS
   COMPONENT  Development)

INSTALL(TARGETS EBTKS ${LIBRARY_INSTALL} DESTINATION lib${LIB_SUFFIX})

set(EBTKS_USE_FILE_CONFIG      "${CMAKE_CURRENT_BINARY_DIR}/UseEBTKS.cmake")
set(EBTKS_INCLUDE_DIRS_CONFIG   ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/templates)
set(EBTKS_LIBRARY_DIRS_CONFIG   ${CMAKE_CURRENT_BINARY_DIR})
set(EBTKS_LIBRARIES_CONFIG      EBTKS)

configure_file(EBTKSConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/EBTKSConfig.cmake @ONLY )

configure_file(UseEBTKS.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/UseEBTKS.cmake @ONLY)

# config for install dir
set(EBTKS_USE_FILE_CONFIG     "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}/UseEBTKS.cmake")
set(EBTKS_INCLUDE_DIRS_CONFIG  ${CMAKE_INSTALL_PREFIX}/include)
set(EBTKS_LIBRARY_DIRS_CONFIG  ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})

IF(SUPERBUILD_STAGING_PREFIX)
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" EBTKS_INCLUDE_DIRS_CONFIG  "${EBTKS_INCLUDE_DIRS_CONFIG}")
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" EBTKS_LIBRARY_DIRS_CONFIG  "${EBTKS_LIBRARY_DIRS_CONFIG}")
ENDIF(SUPERBUILD_STAGING_PREFIX)


configure_file(EBTKSConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/EBTKSConfig.cmake @ONLY )

configure_file(UseEBTKS.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/UseEBTKS.cmake @ONLY)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/UseEBTKS.cmake ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/EBTKSConfig.cmake
  DESTINATION lib${LIB_SUFFIX}
  COMPONENT
   Development)

# testing
IF(BUILD_TESTING)
   ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
