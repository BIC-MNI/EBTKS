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
)

# Which BLAS/LAPACK legacy/EBTKS and legacy/N3 link. The choice is explicit and never
# inferred: nothing here is searched for, and a named backend that cannot be satisfied is a
# configure-time error, not a fallback to bundled. See /app/PLAN.md §10, /app/TODO_BLAS.md.
SET(EBTKS_BLAS_BACKEND "bundled" CACHE STRING
    "Which LAPACK/BLAS legacy/EBTKS and legacy/N3 link: bundled (default, compiles \
clapack/), lapack (link a named Fortran-ABI library directly, no shim), lapacke (link via \
a LAPACKE C shim), or cblas (bundled LAPACK over a named CBLAS shim).")
SET_PROPERTY(CACHE EBTKS_BLAS_BACKEND PROPERTY STRINGS bundled lapack lapacke cblas)

SET(EBTKS_LAPACK_LIBRARIES "" CACHE STRING
    "External LAPACK/BLAS to link when EBTKS_BLAS_BACKEND != bundled, e.g. 'openblas' or \
an absolute path to libmkl_rt.so. Absolute paths under a superbuild.")
SET(EBTKS_LAPACK_LIBRARY_DIRS "" CACHE STRING
    "Extra -L directories for EBTKS_LAPACK_LIBRARIES, when its entries are bare names.")
SET(EBTKS_LAPACK_INCLUDE_DIRS "" CACHE STRING
    "Extra -I directories for lapacke.h/cblas.h, needed by the lapacke and cblas backends.")

IF(NOT EBTKS_BLAS_BACKEND MATCHES "^(bundled|lapack|lapacke|cblas)$")
  MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=${EBTKS_BLAS_BACKEND} is not one of: bundled, lapack, lapacke, cblas")
ENDIF()
IF(NOT EBTKS_BLAS_BACKEND STREQUAL "bundled" AND NOT EBTKS_LAPACK_LIBRARIES)
  MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=${EBTKS_BLAS_BACKEND} requires EBTKS_LAPACK_LIBRARIES to be set")
ENDIF()

# Validation, not detection: nothing is searched for here. Only the symbol a named backend
# depends on is checked, and only to fail at configure time with a better message than the
# linker's -- this is deliberately not FIND_PACKAGE(BLAS)/FIND_PACKAGE(LAPACK), which locate
# Fortran libraries and would find neither C header. This assumes the named library already
# exists when EBTKS is configured, which a superbuild gives with
# ExternalProject_Add(... DEPENDS OpenBLAS); a superbuild that configures EBTKS first fails
# this check correctly rather than silently.
IF(NOT EBTKS_BLAS_BACKEND STREQUAL "bundled")
  IF(POLICY CMP0075)
    CMAKE_POLICY(SET CMP0075 NEW)
  ENDIF()
  INCLUDE(CheckFunctionExists)
  INCLUDE(CheckIncludeFile)

  SET(CMAKE_REQUIRED_INCLUDES ${EBTKS_LAPACK_INCLUDE_DIRS})
  SET(_ebtks_blas_required_libs "")
  FOREACH(_ebtks_blas_dir ${EBTKS_LAPACK_LIBRARY_DIRS})
    LIST(APPEND _ebtks_blas_required_libs "-L${_ebtks_blas_dir}")
  ENDFOREACH()
  LIST(APPEND _ebtks_blas_required_libs ${EBTKS_LAPACK_LIBRARIES})
  SET(CMAKE_REQUIRED_LIBRARIES ${_ebtks_blas_required_libs})

  IF(EBTKS_BLAS_BACKEND STREQUAL "lapack")
    UNSET(EBTKS_HAVE_DSYSV CACHE)
    CHECK_FUNCTION_EXISTS(dsysv_ EBTKS_HAVE_DSYSV)
    IF(NOT EBTKS_HAVE_DSYSV)
      MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=lapack: dsysv_ not found in "
        "EBTKS_LAPACK_LIBRARIES=${EBTKS_LAPACK_LIBRARIES} "
        "EBTKS_LAPACK_LIBRARY_DIRS=${EBTKS_LAPACK_LIBRARY_DIRS}")
    ENDIF()
    # TBSpline.cc:70 declares dsysv_'s integer arguments as the f2c 8-byte "long int", which
    # is wrong against an ILP64 (64-bit Fortran INTEGER) build -- LP64 only. Reject what can
    # be detected: OpenBLAS's own openblas_config.h, if EBTKS_LAPACK_INCLUDE_DIRS names the
    # directory holding it, records how it was built. Silence, not a positive LP64 check: no
    # portable API reports the interface width of an arbitrary named LAPACK, so anything
    # other than a confirmed ILP64 OpenBLAS passes uncontested and stays undiagnosed.
    IF(EBTKS_LAPACK_INCLUDE_DIRS)
      UNSET(_ebtks_openblas_config CACHE)
      FIND_FILE(_ebtks_openblas_config openblas_config.h PATHS ${EBTKS_LAPACK_INCLUDE_DIRS} NO_DEFAULT_PATH)
      IF(_ebtks_openblas_config)
        FILE(STRINGS "${_ebtks_openblas_config}" _ebtks_openblas_ilp64 REGEX "#define OPENBLAS_USE64BITINT")
        IF(_ebtks_openblas_ilp64)
          MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=lapack: ${_ebtks_openblas_config} defines "
            "OPENBLAS_USE64BITINT -- this OpenBLAS is ILP64, which TBSpline.cc's 8-byte "
            "f2c integer declaration does not support. Use backend lapacke instead, whose "
            "shim converts explicitly.")
        ENDIF()
      ENDIF()
      UNSET(_ebtks_openblas_config CACHE)
    ENDIF()
  ELSEIF(EBTKS_BLAS_BACKEND STREQUAL "lapacke")
    UNSET(EBTKS_HAVE_LAPACKE_H CACHE)
    CHECK_INCLUDE_FILE(lapacke.h EBTKS_HAVE_LAPACKE_H)
    IF(NOT EBTKS_HAVE_LAPACKE_H)
      MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=lapacke: lapacke.h not found in "
        "EBTKS_LAPACK_INCLUDE_DIRS=${EBTKS_LAPACK_INCLUDE_DIRS}")
    ENDIF()
    UNSET(EBTKS_HAVE_LAPACKE_DSYSV_WORK CACHE)
    CHECK_FUNCTION_EXISTS(LAPACKE_dsysv_work EBTKS_HAVE_LAPACKE_DSYSV_WORK)
    IF(NOT EBTKS_HAVE_LAPACKE_DSYSV_WORK)
      MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=lapacke: LAPACKE_dsysv_work not found in "
        "EBTKS_LAPACK_LIBRARIES=${EBTKS_LAPACK_LIBRARIES}")
    ENDIF()
  ELSEIF(EBTKS_BLAS_BACKEND STREQUAL "cblas")
    UNSET(EBTKS_HAVE_CBLAS_H CACHE)
    CHECK_INCLUDE_FILE(cblas.h EBTKS_HAVE_CBLAS_H)
    IF(NOT EBTKS_HAVE_CBLAS_H)
      MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=cblas: cblas.h not found in "
        "EBTKS_LAPACK_INCLUDE_DIRS=${EBTKS_LAPACK_INCLUDE_DIRS}")
    ENDIF()
    UNSET(EBTKS_HAVE_CBLAS_DGEMM CACHE)
    CHECK_FUNCTION_EXISTS(cblas_dgemm EBTKS_HAVE_CBLAS_DGEMM)
    IF(NOT EBTKS_HAVE_CBLAS_DGEMM)
      MESSAGE(FATAL_ERROR "EBTKS_BLAS_BACKEND=cblas: cblas_dgemm not found in "
        "EBTKS_LAPACK_LIBRARIES=${EBTKS_LAPACK_LIBRARIES}")
    ENDIF()
  ENDIF()

  UNSET(CMAKE_REQUIRED_INCLUDES)
  UNSET(CMAKE_REQUIRED_LIBRARIES)
  UNSET(_ebtks_blas_required_libs)
ENDIF()

# The 19 clapack/*.c files. EBTKS_CLAPACK_BUNDLED_SRCS keeps the exact interleaved order the
# unconditional list had before this switch existed, so the bundled archive's member order,
# and therefore its bytes, do not move. EBTKS_CLAPACK_BLAS_SRCS/_LAPACK_SRCS are the same 19
# files split by role, each sublist in its own original relative order, for the backends
# that need only one half -- their concatenation is not the bundled order and must not be
# used in place of EBTKS_CLAPACK_BUNDLED_SRCS.
SET(EBTKS_CLAPACK_BUNDLED_SRCS
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
SET(EBTKS_CLAPACK_BLAS_SRCS
  clapack/dcopy.c
  clapack/dgemm.c
  clapack/dgemv.c
  clapack/dger.c
  clapack/dscal.c
  clapack/dswap.c
  clapack/dsyr.c
  clapack/idamax.c
)
SET(EBTKS_CLAPACK_LAPACK_SRCS
  clapack/dlasyf.c
  clapack/dsysv.c
  clapack/dsytf2.c
  clapack/dsytrf.c
  clapack/dsytrs.c
  clapack/ieeeck.c
  clapack/ilaenv.c
  clapack/lsame.c
  clapack/s_cmp.c
  clapack/s_copy.c
  clapack/xerbla.c
)

IF(EBTKS_BLAS_BACKEND STREQUAL "bundled")
  # dsysv_ resolved from clapack/dsysv.c, which calls the rest of clapack/ internally.
  LIST(APPEND EBTKS_LIB_SRCS ${EBTKS_CLAPACK_BUNDLED_SRCS})
ELSEIF(EBTKS_BLAS_BACKEND STREQUAL "cblas")
  # dsysv_ still resolved from clapack/dsysv.c, but the EBTKS_dcopy/EBTKS_dgemm/etc. it calls
  # internally (clapack/blaswrap.h renames the plain f2c names before dsytrf.c/dsytf2.c/
  # dlasyf.c/dsytrs.c reference them, precisely to avoid clashing with a real linked BLAS) now
  # come from shim/blas_cblas.c, forwarding to the named CBLAS, instead of the 8 bundled
  # clapack/*.c BLAS sources.
  LIST(APPEND EBTKS_LIB_SRCS ${EBTKS_CLAPACK_LAPACK_SRCS} shim/blas_cblas.c)
ELSEIF(EBTKS_BLAS_BACKEND STREQUAL "lapacke")
  # dsysv_ resolved from shim/dsysv_lapacke.c, which calls LAPACKE_dsysv_work.
  LIST(APPEND EBTKS_LIB_SRCS shim/dsysv_lapacke.c)
ENDIF()
# lapack: nothing appended -- dsysv_ resolved directly from EBTKS_LAPACK_LIBRARIES.

IF(NOT EBTKS_BLAS_BACKEND STREQUAL "bundled")
  IF(EBTKS_LAPACK_LIBRARY_DIRS)
    LINK_DIRECTORIES(${EBTKS_LAPACK_LIBRARY_DIRS})
  ENDIF()
  IF(EBTKS_LAPACK_INCLUDE_DIRS)
    INCLUDE_DIRECTORIES(${EBTKS_LAPACK_INCLUDE_DIRS})
  ENDIF()
ENDIF()

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
IF(NOT EBTKS_BLAS_BACKEND STREQUAL "bundled")
  TARGET_LINK_LIBRARIES(EBTKS PUBLIC ${EBTKS_LAPACK_LIBRARIES})
ENDIF()
SET(EBTKS_LIBRARIES EBTKS)

MESSAGE(STATUS "EBTKS_BLAS_BACKEND: ${EBTKS_BLAS_BACKEND}")
IF(NOT EBTKS_BLAS_BACKEND STREQUAL "bundled")
  MESSAGE(STATUS "EBTKS_LAPACK_LIBRARIES: ${EBTKS_LAPACK_LIBRARIES}")
  MESSAGE(STATUS "EBTKS_LAPACK_LIBRARY_DIRS: ${EBTKS_LAPACK_LIBRARY_DIRS}")
  MESSAGE(STATUS "EBTKS_LAPACK_INCLUDE_DIRS: ${EBTKS_LAPACK_INCLUDE_DIRS}")
ENDIF()

INSTALL(FILES ${EBTKS_HEADERS}
   DESTINATION include/EBTKS
   COMPONENT  Development)

INSTALL(TARGETS EBTKS ${LIBRARY_INSTALL} DESTINATION lib${LIB_SUFFIX})

set(EBTKS_USE_FILE_CONFIG      "${CMAKE_CURRENT_BINARY_DIR}/UseEBTKS.cmake")
set(EBTKS_INCLUDE_DIRS_CONFIG   ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/templates)
set(EBTKS_LIBRARY_DIRS_CONFIG   ${CMAKE_CURRENT_BINARY_DIR})
set(EBTKS_LIBRARIES_CONFIG      EBTKS)
# Empty under bundled: dsysv_ is inside libEBTKS.a itself and nothing further is needed.
# EBTKS_LAPACK_INCLUDE_DIRS is deliberately not exported here -- shim/dsysv_lapacke.c and
# shim/blas_cblas.c are the only translation units that need lapacke.h/cblas.h, and both are
# compiled inside this build, so no consumer needs the headers.
set(EBTKS_LAPACK_LIBRARIES_CONFIG    ${EBTKS_LAPACK_LIBRARIES})
set(EBTKS_LAPACK_LIBRARY_DIRS_CONFIG ${EBTKS_LAPACK_LIBRARY_DIRS})

configure_file(EBTKSConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/EBTKSConfig.cmake @ONLY )

configure_file(UseEBTKS.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/UseEBTKS.cmake @ONLY)

# config for install dir
set(EBTKS_USE_FILE_CONFIG     "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}/UseEBTKS.cmake")
set(EBTKS_INCLUDE_DIRS_CONFIG  ${CMAKE_INSTALL_PREFIX}/include)
set(EBTKS_LIBRARY_DIRS_CONFIG  ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
set(EBTKS_LAPACK_LIBRARIES_CONFIG    ${EBTKS_LAPACK_LIBRARIES})
set(EBTKS_LAPACK_LIBRARY_DIRS_CONFIG ${EBTKS_LAPACK_LIBRARY_DIRS})

IF(SUPERBUILD_STAGING_PREFIX)
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" EBTKS_INCLUDE_DIRS_CONFIG  "${EBTKS_INCLUDE_DIRS_CONFIG}")
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" EBTKS_LIBRARY_DIRS_CONFIG  "${EBTKS_LIBRARY_DIRS_CONFIG}")
  # A named EBTKS_LAPACK_LIBRARIES/_LIBRARY_DIRS can itself be an absolute path into the
  # staging tree (the "absolute paths under a superbuild" case in the CACHE STRING help
  # text above); strip it the same way or the installed EBTKSConfig.cmake names a directory
  # that no longer exists once the superbuild finishes.
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" EBTKS_LAPACK_LIBRARIES_CONFIG    "${EBTKS_LAPACK_LIBRARIES_CONFIG}")
  STRING(REPLACE "${SUPERBUILD_STAGING_PREFIX}/" "" EBTKS_LAPACK_LIBRARY_DIRS_CONFIG "${EBTKS_LAPACK_LIBRARY_DIRS_CONFIG}")
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
