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
CCLD = @CXX@

subpackage = query
include $(top_srcdir)/pkgs/subpackage.am

bin_PROGRAMS = cencalvmcquery

cencalvmcquery_SOURCES = \
	cencalvmcquery.c

cencalvmcquery_LDFLAGS =

INCLUDES +=

cencalvmcquery_LDADD = \
	$(top_builddir)/pkgs/query/liblang/libcencalvmquerylang.la

# ----------------------------------------------------------------------

VELMODEL = USGSBayAreaVM-05.1.0.etree
VELMODELEXT = USGSBayAreaVMExt-05.1.0.etree

testmaxres:
	cencalvmcquery \
		-d $(VELMODEL) \
		-i ../data/sample-05.1.x.in \
		-o test_maxres.out \
		-l test.log

testfixedres:
	cencalvmcquery \
		-d $(VELMODEL) \
		-i ../data/sample-05.1.x.in \
		-o test_fixedres.out \
		-t fixedres \
		-r 200.0 \
		-l test.log

testwaveres:
	cencalvmcquery \
		-d $(VELMODEL) \
		-i ../data/sample-05.1.x.in \
		-o test_waveres.out \
		-t waveres \
		-r 0.2 \
		-l test.log

testextmaxres:
	cencalvmcquery \
		-d $(VELMODEL) \
		-e $(VELMODELEXT) \
		-i ../data/sample_ext-05.1.x.in \
		-o testext_maxres.out \
		-l test.log

testextfixedres:
	cencalvmcquery \
		-d $(VELMODEL) \
		-e $(VELMODELEXT) \
		-i ../data/sample_ext-05.1.x.in \
		-o testext_fixedres.out \
		-t fixedres \
		-r 200.0 \
		-l test.log

testextwaveres:
	cencalvmcquery \
		-d $(VELMODEL) \
		-e $(VELMODELEXT) \
		-i ../data/sample_ext-05.1.x.in \
		-t waveres \
		-r 0.2 \
		-o testext_waveres.out \
		-l test.log

CLEANFILES = \
	test_maxres.out \
	test_fixedres.out \
	test_waveres.out \
	test.log \
	testext_maxres.out \
	testext_fixedres.out \
	testext_waveres.out

# version
# $Id: Makefile.am,v 1.5 2005/10/29 18:52:20 brad Exp $

# End of file 
