#
# Makefile for GENIE/Reweight Apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

# Mac OS X -bind_at_load means we must list the low level libraries
# first, and each stage must fully resolved before the next library
LIBRARIES  := $(LIBRARIES) $(CERN_LIBRARIES) $(GENIE_LIBS)

TGT_BASE =  grwght1p   \
            grwghtnp   

TGT = $(addprefix $(GENIE_BIN_PATH)/,$(TGT_BASE))

# this is a nasty-looking construction which runs ldd on an existing library
# to figure out what other libraries it depends on.
# (sort of like g++ -M)
find_libs = ${shell ldd $(GENIE_BIN_PATH)/$(1) 2>&1 | sed -n 's/.*=> \(.*\) (.*/\1/pg' | tr '\n' ' '}

all: $(TGT)

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -MMD -MP -c $(CPP_INCLUDES) $< -o $@

# CLEANING-UP

purge: FORCE
	$(RM) *.o *~ *.d core

clean: FORCE
	$(RM) *.o *~ *.d core
	$(RM) $(TGT)

distclean: FORCE
	$(RM) $(patsubst $(GENIE_BIN_DIR),$(GENIE_BIN_INSTALLATION_PATH),$(TGT))

FORCE:

# DO NOT DELETE

# this should be last because it includes new targets...
-include $(wildcard *.d)

