# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################


# User specifications
# 1. Set PDAFDIR, the directory where the PDAF package resides
#    either in this file, in the 'make' command, or as environment variable
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory PDAFDIR/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package
#PDAFDIR = ../..

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(PDAFDIR)/make.arch/$(PDAF_ARCH).h


# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Source files
SRC_MODEL =  	model_parallel_mod.F90 \
		model_mod.F90 \
		model_init_mod.F90 \
		model_step_mod.F90 \
		model_post_mod.F90 \
		main.F90

#######################################################
# object files
OBJ_MODEL := $(SRC_MODEL:%.F90=%.o)

######################################################

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D model with parallelization";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  model         - just the forward model  (use without -DUSE_PDAF)";
	@echo "";
	@echo "Depending on whether PDAF should be active, one has to define";
	@echo "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo  "List available settings for PDAF_ARCH with 'make listarch'";
	@echo "Path to PDAF (PDAFDIR): " $(PDAFDIR);
	@echo "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo "------------------------------------------------------------------";


model : $(OBJ_MODEL)
	$(LD)  $(OPT_LNK)  -o $@  $(OBJ_MODEL)

######################################################
# use pattern rule to create rules for all object files

%.o: %.F90
	$(FC) $(OPT) $(CPP_DEFS) -I$(PDAFDIR)/include  -c $*.F90

######################################################
# Cleans

cleanall : clean cleandata

clean :
	rm -f *.o *.mod *.a model

cleandata:
	rm -f ens*ana.txt ens*for.txt state_*ana.txt state_*for.txt true*.txt

cleandataq:
	@rm -f ens*ana.txt ens*for.txt state_*ana.txt state_*for.txt

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 $(PDAFDIR)/make.arch | cut -d"." -f1
