SUBDIRS = src templates include

AUTOMAKE_OPTIONS = check-news

lib_LTLIBRARIES = libEBTKS.la

libEBTKS_la_SOURCES = version.cc

libEBTKS_la_LIBADD = \
	$(top_builddir)/src/libAZgen_s.la \
	$(top_builddir)/templates/libAZgen_t.la


