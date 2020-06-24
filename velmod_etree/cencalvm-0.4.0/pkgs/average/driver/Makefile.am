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

subpackage = average
include $(top_srcdir)/pkgs/subpackage.am

bin_PROGRAMS = avgcencalvm

avgcencalvm_SOURCES = \
	avgcencalvm.cc

INCLUDES += 

avgcencalvm_LDADD = \
	$(top_builddir)/pkgs/average/lib/libcencalvmaverage.la \
	$(top_builddir)/pkgs/storage/lib/libcencalvmstorage.la

# version
# $Id: Makefile.am,v 1.1 2005/10/10 17:28:13 brad Exp $

# End of file 
