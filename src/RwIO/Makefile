#
# Makefile for GENIE/Reweight package
#
# Author: Costas Andreopoulos <c.andreopoulos \at cern.ch>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE_REWEIGHT)/src/make/Make.include

PACKAGE           = RwIO
PACKAGE_ABBREV    = RwIO
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE_REWEIGHT)/src/make/Make.targets

FORCE:
