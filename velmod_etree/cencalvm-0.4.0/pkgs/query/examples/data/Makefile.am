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

data_DATA = \
	sample-05.0.0.out \
	sample-05.0.1.out \
	sample-05.0.1_fixedres.out \
	sample-05.0.1_waveres.out \
	sample-05.0.x.in \
	sample-05.1.0_maxres.out \
	sample-05.1.0_fixedres.out \
	sample-05.1.0_waveres.out \
	sample-05.1.x.in \
	sample_ext-05.1.x.in \
	sample_ext-05.1.0_maxres.out \
	sample_ext-05.1.0_fixedres.out \
	sample_ext-05.1.0_waveres.out

noinst_HEADERS = \
	$(data_DATA)

# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/pkgs/query/examples/data
export-data: $(data_DATA)
	for f in $(data_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	for f in $(data_DATA); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done

# export
clean-local: clean-data
BUILT_SOURCES = export-data

# version
# $Id: Makefile.am,v 1.3 2005/10/29 18:52:20 brad Exp $

# End of file 
