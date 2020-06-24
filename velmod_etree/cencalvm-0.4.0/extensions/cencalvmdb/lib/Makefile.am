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

subpackage = extensions/cencalvmdb
include $(top_srcdir)/extensions/subpackage.am

lib_LTLIBRARIES = libcencalvmdb.la

libcencalvmdb_la_SOURCES = \
	CenCalVMDB.cc

libcencalvmdb_la_LDFLAGS =

INCLUDES +=

libcencalvmdb_la_LIBADD = \
	$(top_builddir)/pkgs/query/lib/libcencalvmquery.la \
	-lspatialdata -lspatialdatageocoords

subpkginclude_HEADERS = \
	CenCalVMDB.h \
	CenCalVMDB.icc

# version
# $Id: Makefile.am,v 1.2 2005/10/26 23:26:45 brad Exp $

# End of file 
