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

subpackage = create
include $(top_srcdir)/pkgs/subpackage.am

lib_LTLIBRARIES = libcencalvmcreate.la

libcencalvmcreate_la_SOURCES = \
	VMCreator.cc \
	GridIngester.cc

libcencalvmcreate_la_LDFLAGS =

INCLUDES +=

libcencalvmcreate_la_LIBADD = \
	-letree \
	$(top_builddir)/pkgs/storage/lib/libcencalvmstorage.la

subpkginclude_HEADERS = \
	VMCreator.h \
	VMCreator.icc

noinst_HEADERS = \
	GridIngester.h

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# version
# $Id: Makefile.am,v 1.1 2005/10/10 17:30:18 brad Exp $

# End of file 
