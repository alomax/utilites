# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

subpackage = query
include $(top_srcdir)/pkgs/subpackage.am

lib_LTLIBRARIES = libcencalvmquery.la

libcencalvmquery_la_SOURCES = \
	VMQuery.cc

libcencalvmquery_la_LDFLAGS =

INCLUDES +=

libcencalvmquery_la_LIBADD = \
	$(top_builddir)/pkgs/storage/lib/libcencalvmstorage.la \
	-letree

subpkginclude_HEADERS = \
	VMQuery.h \
	VMQuery.icc

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# version
# $Id: Makefile.am,v 1.2 2005/10/25 15:46:44 brad Exp $

# End of file 
