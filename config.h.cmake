/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H 1

/* Define to 1 if you have the `finite' function. */
#cmakedefine HAVE_FINITE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H 1

/* Define to 1 if you have the `isfinite' function. */
#cmakedefine HAVE_ISFINITE 1

/* Define to 1 if you have the <malloc.h> header file. */
#cmakedefine HAVE_MALLOC_H 1

/* Define to 1 if you have the <memory.h> header file. */
#cmakedefine HAVE_MEMORY_H 1

/* Define to 1 if you have the `mkstemp' function. */
#cmakedefine HAVE_MKSTEMP 1

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#cmakedefine HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H 1

/* Define to 1 if you have the <values.h> header file. */
#cmakedefine HAVE_VALUES_H 1

/* Define to 1 if you have the ANSI C header files. */
#cmakedefine STDC_HEADERS 1

#define PACKAGE "@PACKAGE@"
#define PACKAGE_NAME "@PACKAGE_NAME@"
#define PACKAGE_BUGREPORT "@PACKAGE_BUGREPORT@"
#define PACKAGE_VERSION "@PACKAGE_VERSION@"
#define PACKAGE_STRING "@PACKAGE_STRING@"
#define VERSION PACKAGE_VERSION


#if defined(HAVE_ISFINITE)
#define FINITE(x) isfinite(x)
#elif defined(HAVE_FINITE)
#define FINITE(x) finite(x)
#else
#error "Need isfinite or finite"
#endif /* HAVE_ISFINITE */
