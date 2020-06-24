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

SUBDIRS = data

TESTS = testcreate

check_PROGRAMS = testcreate

testcreate_SOURCES = \
	TestVMCreator.cc \
	testcreate.cc

noinst_HEADERS = \
	TestVMCreator.h

testcreate_LDFLAGS =

INCLUDES +=

testcreate_LDADD = \
	-lcppunit -ldl \
	-letree \
	$(top_builddir)/pkgs/create/lib/libcencalvmcreate.la \
	$(top_builddir)/pkgs/storage/lib/libcencalvmstorage.la

# version
# $Id: Makefile.am,v 1.1 2005/10/10 17:30:18 brad Exp $

# End of file 
