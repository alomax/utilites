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

# Link using C++ compiler to make sure we satisfy all of the C++
# symbols in the C++ query library which this program uses
F77LD = @CXX@

subpackage = query
include $(top_srcdir)/pkgs/subpackage.am

bin_PROGRAMS = cencalvmf77query

cencalvmf77query_SOURCES = \
	cencalvmf77query.f

cencalvmf77query_LDFLAGS =

INCLUDES +=

cencalvmf77query_LDADD = \
	$(top_builddir)/pkgs/query/liblang/libcencalvmquerylang.la \
	$(FCLIBS)

CLEANFILES = \
	test_maxres.out \
	test_fixedres.out \
	test_waveres.out \
	test.log \
	testext_maxres.out \
	testext_fixedres.out \
	testext_waveres.out


# version
# $Id: Makefile.am,v 1.4 2005/10/29 18:52:20 brad Exp $

# End of file 
