AUTOMAKE_OPTIONS = check-news

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/templates -DUSE_COMPMAT -DUSE_DBLMAT

lib_LTLIBRARIES = libEBTKS.la

# Some #include'd C++ sources which are not independently compiled.
#
noinst_HEADERS = \
	templates/SimpleArraySpec.cc \
	templates/MatrixSpec.cc

# Headers which are installed.
#
include_ebtksdir = $(includedir)/EBTKS

include_ebtks_HEADERS = \
	include/amoeba.h \
	include/assert.h \
	include/backProp.h \
	include/Complex.h \
	include/dcomplex.h \
	include/fcomplex.h \
	include/FileIO.h \
	include/Histogram.h \
	include/matlabSupport.h \
	include/Minc.h \
	include/MPoint.h \
	include/MString.h \
	include/MTypes.h \
	include/OpTimer.h \
	include/OrderedCltn.h \
	include/Path.h \
	include/Polynomial.h \
	include/popen.h \
	include/Spline.h \
	include/TBSpline.h \
	include/TrainingSet.h \
	include/trivials.h \
	templates/Array.h \
	templates/CachedArray.h \
	templates/Dictionary.h \
	templates/Matrix3D.h \
	templates/Matrix.h \
	templates/MatrixSupport.h \
	templates/MatrixTest.h \
	templates/miscTemplateFunc.h \
	templates/Pool.h \
	templates/SimpleArray.h \
	templates/Stack.h \
	templates/ValueMap.h

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
	src/Spline.cc \
	src/TBSpline.cc \
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
        templates/MatrixConv.cc



