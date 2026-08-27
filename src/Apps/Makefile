#
# Makefile for GENIE/Reweight Apps
#
# Costas Andreopoulos <c.andreopoulos \at cern.ch>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE_REWEIGHT)/src/make/Make.include

GENIE_LIBS           = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
GENIE_REWEIGHT_LIBS  = $(shell $(GENIE_REWEIGHT)/src/scripts/setup/genie-reweight-config --libs)

# Mac OS X -bind_at_load means we must list the low level libraries
# first, and each stage must fully resolved before the next library
LIBRARIES  := $(LIBRARIES) $(CERN_LIBRARIES) $(GENIE_LIBS) $(GENIE_REWEIGHT_LIBS)

TGT_BASE =  grwght1p   \
            grwghtnp

ifeq ($(strip $(GOPT_ENABLE_PROFESSOR2)),YES)
TGT_BASE += grwproftest
endif

TGT = $(addprefix $(GENIE_REWEIGHT_BIN_PATH)/,$(TGT_BASE))

# this is a nasty-looking construction which runs ldd on an existing library
# to figure out what other libraries it depends on.
# (sort of like g++ -M)
find_libs = ${shell ldd $(GENIE_BIN_PATH)/$(1) 2>&1 | sed -n 's/.*=> \(.*\) (.*/\1/pg' | tr '\n' ' '}

all: $(TGT)


# simple utility for calculating weights in response to tweaking a single systematic param
#
$(GENIE_REWEIGHT_BIN_PATH)/grwght1p: gRwght1Param.o $(call find_libs,grwght1p)
	@echo "** Building grwght1p"
	$(LD) $(LDFLAGS) gRwght1Param.o $(LIBRARIES) -o $(GENIE_REWEIGHT_BIN_PATH)/grwght1p

# simple utility for calculating weights in response to tweaking a multiple (possibly correlated) systematic param
#
$(GENIE_REWEIGHT_BIN_PATH)/grwghtnp: gRwghtNCorrelatedParams.o $(call find_libs,grwghtnp)
	@echo "** Building grwghtnp"
	$(LD) $(LDFLAGS) gRwghtNCorrelatedParams.o $(LIBRARIES) -o $(GENIE_REWEIGHT_BIN_PATH)/grwghtnp

$(GENIE_REWEIGHT_BIN_PATH)/grwproftest: gRwProfTest.o $(call find_libs,grwghtnp)
	@echo "** Building grwghtnp"
	$(LD) $(LDFLAGS) gRwProfTest.o $(LIBRARIES) -o $(GENIE_REWEIGHT_BIN_PATH)/grwproftest


%.o : %.cxx
	$(CXX) $(CXXFLAGS) -MMD -MP -c $(CPP_INCLUDES) $< -o $@

# CLEANING-UP

purge: FORCE
	$(RM) *.o *~ *.d core

clean: FORCE
	$(RM) *.o *~ *.d core
	$(RM) $(TGT)

distclean: FORCE
	$(RM) $(patsubst $(GENIE_REWEIGHT_BIN_PATH)%,$(GENIE_REWEIGHT_BIN_INSTALLATION_PATH)%,$(TGT))

FORCE:

# DO NOT DELETE

# this should be last because it includes new targets...
-include $(wildcard *.d)

