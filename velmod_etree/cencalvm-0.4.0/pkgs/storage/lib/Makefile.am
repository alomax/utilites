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

subpackage = storage
include $(top_srcdir)/pkgs/subpackage.am

lib_LTLIBRARIES = libcencalvmstorage.la

libcencalvmstorage_la_SOURCES = \
	ErrorHandler.cc \
	Geometry.cc \
	Payload.cc \
	Projector.cc

libcencalvmstorage_la_LDFLAGS =

INCLUDES +=

libcencalvmstorage_la_LIBADD = \
	-lproj \
	-letree

subpkginclude_HEADERS = \
	ErrorHandler.h \
	ErrorHandler.icc \
	Geometry.h \
	Geometry.icc \
	Payload.h \
	etreefwd.h

noinst_HEADERS = \
	Projector.h \
	proj4fwd.h

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# version
# $Id: Makefile.am,v 1.2 2005/10/11 18:27:42 brad Exp $

# End of file 
