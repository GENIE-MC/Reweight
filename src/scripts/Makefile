# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE_REWEIGHT)/src/make/Make.include

doxygen: FORCE
	$(DOXYGEN) dox/genie-doxygen;

install: FORCE
	cp setup/genie-reweight-config $(GENIE_REWEIGHT_BIN_PATH); \
	chmod 755 $(GENIE_REWEIGHT_BIN_PATH)/genie-reweight-config; \


#
#install_db: FORCE
#	cp setup/gnuvld_dbbootstrap $(GENIE_BIN_PATH); \
#	chmod 755 $(GENIE_BIN_PATH)/gnuvld_dbbootstrap
#

clean: FORCE
	$(RM) $(GENIE_REWEIGHT_BIN_PATH)/genie-reweight-config

#
#clean_db: FORCE
#	$(RM) $(GENIE_BIN_PATH)/gnuvld_dbbootstrap 
#

distclean: FORCE
	$(RM) $(GENIE_REWEIGHT_BIN_INSTALLATION_PATH)/genie-reweight-config

FORCE:

