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

noinst_HEADERS = \
	CHANGES \
	README \
	INSTALL \
	TODO

SUBDIRS = \
	pkgs \
	extensions

if ENABLE_DOCUMENTATION
  SUBDIRS += doc
endif


DISTCHECK_CONFIGURE_FLAGS = \
	--enable-testing \
	--enable-documentation \
	CPPFLAGS="-I${ETREE_INCDIR} -I${PROJ4_INCDIR} -I${CPPUNIT_INCDIR}" \
	LDFLAGS="-L${ETREE_LIBDIR} -L${PROJ4_LIBDIR} -L${CPPUNIT_LIBDIR}"

# version
# $Id: Makefile.am,v 1.3 2005/11/08 18:21:58 brad Exp $

# End of file 
