dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([src/FileIO.cc])
AM_INIT_AUTOMAKE(EBTKS, 0.99.1)
AM_CONFIG_HEADER(config.h)

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

dnl Build only static libs by default
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

smr_WITH_BUILD_PATH

# Check for required C libraries
mni_REQUIRE_VOLUMEIO
AC_CHECK_FUNCS(finite isfinite mkstemp)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

