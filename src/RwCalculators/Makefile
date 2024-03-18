#
# Makefile for GENIE/Reweight package
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE_REWEIGHT)/src/make/Make.include

PACKAGE           = RwCalculators
PACKAGE_ABBREV    = RwClc
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    = $(ROOT_LIBRARIES) -L$(GENIE_REWEIGHT_LIB_PATH) -lGRwFwk \
  -L$(GENIE_LIB_PATH) -lGFwGHEP -lGFwAlg -lGPhXSIg -lGPhHadTransp

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE_REWEIGHT)/src/make/Make.targets

FORCE:
