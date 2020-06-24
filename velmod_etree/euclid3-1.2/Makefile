# -*- Makefile -*-
#
# $Id: Makefile,v 1.13 2005/08/05 00:40:07 tutk Exp $
# Author:	Julio Lopez <jclopez@cs.cmu.edu>
#
# Top level makefile for the etree library.

# Variable assigment:
#
# Override the following variables in user.make or make command line to
# reflect local settings
#

SHELL   = /bin/sh

PACKAGE = euclid3
VERSION = 1.2
TAR	= tar

TMPDIR  = /tmp/


#
# where the distribution files are published
#
PUBDIR  = /afs/cs/project/euclid/www

-include user.make

MODULES = cvm examples output4d tools

#
# Targets and rules
#

default: all


#
# this rule specifies how modules (subdirs) are built
#
libsrc $(MODULES):
	$(MAKE) -C $@

all:	libsrc tools examples

$(MODULES): libsrc


#
# target to clean subdirectories
#
clean cleanall:
	$(foreach dir, libsrc $(MODULES), $(MAKE) -C $(dir) $@;)
	-rm -rf *~


#
# Rules and targets below are for the maintainers to create distribution
# archives
#

usecvs  := true

#
# if no CVSROOT  defined, get a valid CVSROOT from the content of CVS/Root
#
#CVSROOT = bar
ifeq ($(strip $(CVSROOT)),)
    CVSROOT =
    ifneq ($(strip $(wildcard CVS/Root)),)
	CVSROOT := $(shell cat CVS/Root)
    endif

    ifeq ($(strip $(CVSROOT)),)
	CVSROOT := "CVSROOT variable not defined"
	usecvs  := false
    endif
endif


.PHONY:	all cleanall tar cvstag cvstag-latest pub cvs-distdir local-distdir \
	dist libsrc cvm tools examples meshgen output4d


srcdir  = $(CURDIR)
distdir = $(PACKAGE)-$(VERSION)


ifeq ($(usecvs),true)
distdir: cvs-distdir
else
distdir: local-distdir
endif


# replace dot (.) with underscore (_)
ver_str = $(subst .,_,$(VERSION))

# This target creates a tag of the form $(PACKAGE)-$(ver_str),
# for the current version of the source files. E.g., euclid3-1_1
cvstag:
	cvs tag -F $(PACKAGE)-$(ver_str)

cvstag-latest:
	cvs rtag -a -d LATEST-DIST $(PACKAGE)
	cvs tag -F LATEST-DIST


cvsroot:
	echo CVSROOT=$(CVSROOT)

#
# create dist dir using cvs export
#
cvs-distdir: cvstag cvstag-latest
	[ -n "$(TMPDIR)" ] && cd $(TMPDIR); \
	cvs -d $(CVSROOT) export -r LATEST-DIST \
	    -d $(distdir) $(PACKAGE)
	-rm -rf $(TMPDIR)$(distdir)/meshgen
	-rm -f $(TMPDIR)$(distdir)/cvm/Version3.0.tar.gz

#
# create dist dir by copying source dir
#
local-distdir: cleanall
	-rm -rf $(TMPDIR)$(distdir)
	mkdir -p $(TMPDIR)$(distdir)
	chmod 777 $(TMPDIR)$(distdir)

	-cp $(srcdir)/Makefile $(TMPDIR)$(distdir)
	-cp $(srcdir)/README $(TMPDIR)$(distdir)
	for d in $(MODULES) libsrc; do \
	   echo copying $$d; \
	   cp -r $(srcdir)/$$d $(TMPDIR)$(distdir); \
	done
	-rm -f $(TMPDIR)$(distdir)/cvm/Version3.0.tar.gz
	find $(TMPDIR)$(distdir) -name CVS -prune -exec rm \{\}/Root \
	\{\}/Repository \{\}/Entries \; -exec rmdir \{\} \;


#
# create distribution archives
#
dist:	distdir
	chmod -R a+r $(TMPDIR)$(distdir)
	-rm -f $(TMPDIR)$(distdir).tar*
	[ -n "$(TMPDIR)" ] && cd $(TMPDIR); \
	$(TAR) chof $(TMPDIR)$(distdir).tar $(distdir)
	rm -rf $(TMPDIR)$(distdir)
	-gzip -9 -c $(TMPDIR)$(distdir).tar > $(TMPDIR)$(distdir).tar.gz
	-bzip2 -k -c  $(TMPDIR)$(distdir).tar > $(TMPDIR)$(distdir).tar.bz2


# This target untars the dist file and tries a VPATH configuration.  Then
# it guarantees that the distribution is self-contained by making another
# tarfile.
distcheck: dist
	-rm -rf $(TMPDIR)$(distdir)
	[ -n $(TMPDIR) ] && cd $(TMPDIR); \
	$(TAR) zxf $(TMPDIR)$(distdir).tar.gz
	cd $(TMPDIR)$(distdir) \
	  && $(MAKE) \
	  && $(MAKE) TMPDIR= dist
	rm -rf $(TMPDIR)$(distdir)
	@echo "========================"; \
	echo "$$(TMPDIR)(distdir).tar.gz is ready for distribution"; \
	echo "========================"

#
# update the archives in the distribution directory
#
pub: dist
	mv -f $(TMPDIR)$(distdir).tar* $(PUBDIR)
	-rm -f $(PUBDIR)/$(PACKAGE)-latest.tar*
	cd $(PUBDIR); \
	ln -s $(distdir).tar $(PACKAGE)-latest.tar; \
	ln -s $(distdir).tar.gz $(PACKAGE)-latest.tar.gz; \
	ln -s $(distdir).tar.bz2 $(PACKAGE)-latest.tar.bz2



# $Id: Makefile,v 1.13 2005/08/05 00:40:07 tutk Exp $
