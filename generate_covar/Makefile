#####################################################################
# Generic Makefile for to build tool program for PDAF
#
# Set the following enviroment variables
# 1. PDAFDIR=     ... the directory where the PDAF package resides
# 2. PDAF_ARCH=   ... file including the compile definitions
#                    (See directory PDAFDIR/make.arch for files.
#		      PDAF_ARCH is filename without .h)
#
# Example
# make offline_pdaf PDAF_ARCH=linux_gfortran_openmpi PDAFDIR=../PDAF
#####################################################################

# Include machine-specific definitions
# For available include files see PDAF directory make.arch
# To choose a file, set PDAF_ARCH for make or by an
# environment variable.
include $(PDAFDIR)/make.arch/$(PDAF_ARCH).h

######################################################

SRC_COVAR  =	model_pdaf_mod_offline.F90 \
		assimilation_pdaf_mod.F90 \
		parallel_pdaf_mod.F90 \
		initialize_grid_mod.F90 \
		statevector_pdaf_mod.F90 \
		parser_mpi.F90 \
		init_parallel_pdaf_offline.F90 \
		generate_covar.F90


#######################################################
# object files
OBJ_COVAR := $(SRC_COVAR:%.F90=build/%.o)

######################################################

.phony: info
info:
	@echo "Makefile to build program generating a covariance matrix file";
	@echo  "-------------------------------------------------------------------------------";
	@echo  "Use as   make TARGET   where TARGET can be:";
	@echo  "  generate_covar - Generate a file for the covariance matrix from a trajectory";
	@echo  "  clean          - Clean up";
	@echo  "";
	@echo  "Path to PDAF (PDAFDIR): " $(PDAFDIR);
	@echo  "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo  "-------------------------------------------------------------------------------";

generate_covar: libpdaf-d.a $(OBJ_COVAR)
	rm -f $@
	$(LD)  $(OPT) -o $@ $(OBJ_COVAR) $(NC_LIB) \
	-L$(PDAFDIR)/lib -lpdaf-d $(LINK_LIBS)

######################################################
# Build PDAF library

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(PDAFDIR); make;


######################################################
# use pattern rule to create rules for all object files

build/%.o: src/%.F90
	$(FC) $(OPT) $(CPP_DEFS) $(NC_INC) -I$(PDAFDIR)/include -c -o $@ $(MODULEOPT) build/ $<

######################################################
# Cleans

.phony: clean
clean :
	rm -f build/*.o build/*.mod generate_covar
