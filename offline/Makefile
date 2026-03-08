# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################


# User specifications
# 1. Set PDAFDIR, the directory where the PDAF package resides
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

# Moduls used for PDAF
MOD_ASSIM = 	model_pdaf_mod_offline.F90 \
		parallel_pdaf_mod.F90 \
		assimilation_pdaf_mod.F90 \
		statevector_pdaf_mod.F90 \
		initialize_grid_mod.F90 \
		io_pdaf_mod.F90 \
		main_offline.F90 \
		parser_mpi.F90 


# Interface to PDAF - model sided
SRC_PDAF_INT =  init_parallel_pdaf_offline.F90 \
		init_pdaf_offline.F90 \
		assimilate_pdaf_offline.F90 \
		finalize_pdaf.F90 \

# user-supplied call-back routines
SRC_CB  =       init_ens_pdaf_offline.F90 \
		prepoststep_pdaf_offline.F90

# User-supplied call-back routines for localization
SRC_CB_LOCAL =  init_n_domains_pdaf.F90 \
		init_dim_l_pdaf.F90

# Routines for PDAF-OMI observation handling
SRC_CB_OMI = 	synobs_pdaf_mod.F90 \
		obs_A_pdafomi.F90 \
		obs_B_pdafomi.F90 \
		callback_obs_pdafomi.F90

# Additional routines call by other routines
SRC_MISC =      init_pdaf_parse.F90

# Full list of user-supplied routines for online modes
SRC_PDAF_USER = $(MOD_ASSIM) $(SRC_CB_OMI) $(SRC_MISC) $(SRC_PDAF_INT) $(SRC_CB) $(SRC_CB_LOCAL)

#######################################################
# object files
#OBJ_MODEL := $(SRC_MODEL:%.F90=%.o)
OBJ_PDAF_USER := $(SRC_PDAF_USER:%.F90=%.o)

######################################################

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D model with parallelization";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  offline_pdaf   - PDAF offline DA program";
	@echo "  (if compiled without -DUSE_PDAF, both programs behave equally)";
	@echo "";
	@echo "Depending on whether PDAF should be active, one has to define";
	@echo "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo  "List available settings for PDAF_ARCH with 'make listarch'";
	@echo "Path to PDAF (PDAFDIR): " $(PDAFDIR);
	@echo "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo "------------------------------------------------------------------";



offline_pdaf : libpdaf-d.a $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@ \
	$(OBJ_PDAF_USER) \
	-L$(PDAFDIR)/lib -lpdaf-d $(LINK_LIBS) $(NC_LIB)

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(PDAFDIR); make;


######################################################

# use pattern rule to create rules for all object files

%.o: %.F90
	$(FC) $(OPT) $(CPP_DEFS) $(NC_INC) -I$(PDAFDIR)/include  -c $*.F90

######################################################
# Cleans

cleanall : cleanpdaf clean cleandata

clean :
	rm -f *.o *.mod *.a offline_pdaf

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(PDAFDIR); make clean

cleandata:
	rm -f *.nc

cleandataq:
	@rm -f *.nc

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1
