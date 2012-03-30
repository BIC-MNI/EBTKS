AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/templates -DUSE_COMPMAT -DUSE_DBLMAT

lib_LTLIBRARIES = libEBTKS.la


EXTRA_DIST = $(m4_files) \
	epm-header.in

# Some #include'd C++ sources which are not independently compiled.
#
noinst_HEADERS = \
	templates/SimpleArraySpec.cc \
	templates/MatrixSpec.cc \
	clapack/f2c.h \
	clapack/blaswrap.h

# Headers which are installed.
#
include_ebtksdir = $(includedir)/EBTKS

include_ebtks_HEADERS = \
	include/EBTKS/amoeba.h \
	include/EBTKS/assert.h \
	include/EBTKS/backProp.h \
	include/EBTKS/Complex.h \
	include/EBTKS/dcomplex.h \
	include/EBTKS/fcomplex.h \
	include/EBTKS/FileIO.h \
	include/EBTKS/Histogram.h \
	include/EBTKS/matlabSupport.h \
	include/EBTKS/Minc.h \
	include/EBTKS/MPoint.h \
	include/EBTKS/MString.h \
	include/EBTKS/MTypes.h \
	include/EBTKS/OpTimer.h \
	include/EBTKS/OrderedCltn.h \
	include/EBTKS/Path.h \
	include/EBTKS/Polynomial.h \
	include/EBTKS/popen.h \
	include/EBTKS/TrainingSet.h \
	include/EBTKS/trivials.h \
	templates/EBTKS/Array.h \
	templates/EBTKS/CachedArray.h \
	templates/EBTKS/Dictionary.h \
	templates/EBTKS/Matrix3D.h \
	templates/EBTKS/Matrix.h \
	templates/EBTKS/MatrixSupport.h \
	templates/EBTKS/MatrixTest.h \
	templates/EBTKS/miscTemplateFunc.h \
	templates/EBTKS/Pool.h \
	templates/EBTKS/SimpleArray.h \
	templates/EBTKS/Stack.h \
	templates/EBTKS/ValueMap.h

# Source files
#
libEBTKS_la_SOURCES = version.cc \
	src/FileIO.cc \
	src/Histogram.cc \
	src/MPoint.cc \
	src/MString.cc \
	src/OpTimer.cc \
	src/OrderedCltn.cc \
	src/Path.cc \
	src/Polynomial.cc \
	src/TrainingSet.cc \
	src/amoeba.cc \
	src/backProp.cc \
	src/fcomplex.cc \
	src/dcomplex.cc \
        templates/Array.cc \
        templates/CachedArray.cc \
        templates/Dictionary.cc \
        templates/Matrix.cc \
        templates/MatrixSupport.cc \
        templates/Pool.cc \
        templates/SimpleArray.cc \
        templates/ValueMap.cc \
        templates/miscTemplateFunc.cc \
        templates/MatrixConv.cc \
	clapack/dcopy.c \
	clapack/dgemm.c \
	clapack/dgemv.c \
	clapack/dger.c \
	clapack/dlasyf.c \
	clapack/dscal.c \
	clapack/dswap.c \
	clapack/dsyr.c \
	clapack/dsysv.c \
	clapack/dsytf2.c \
	clapack/dsytrf.c \
	clapack/dsytrs.c \
	clapack/idamax.c \
	clapack/ieeeck.c \
	clapack/ilaenv.c \
	clapack/lsame.c \
	clapack/s_cmp.c \
	clapack/s_copy.c \
	clapack/xerbla.c


m4_files = m4/mni_REQUIRE_LIB.m4		\
	m4/mni_REQUIRE_MNILIBS.m4		\
	m4/mni_REQUIRE_OPENINVENTOR.m4		\
	m4/mni_cxx_have_koenig_lookup.m4	\
	m4/smr_CGAL_MAKEFILE.m4			\
	m4/smr_OPTIONAL_LIB.m4			\
	m4/smr_REQUIRED_LIB.m4			\
	m4/smr_WITH_BUILD_PATH.m4

