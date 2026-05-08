#####################################################################
# Makefile to build PDAF offline mode assimilation program
#
# Set the following enviroment variables
# 1. PDAFDIR=     ... the directory where the PDAF package resides
# 2. PDAF_ARCH=   ... file including the compile definitions
#                    (See directory PDAFDIR/make.arch for files.
#		      PDAF_ARCH is filename without .h)
#
# Example
# make offline_pdaf PDAF_ARCH=linux_gfortran_openmpi PDAFDIR=../pdaf
#####################################################################

# Include machine-specific definitions
# For available include files see PDAF directory make.arch
# To choose a file, set PDAF_ARCH for make or by an
# environment variable.
include $(PDAFDIR)/make.arch/$(PDAF_ARCH).h


######################################################

# Generic routines and modules used for PDAF
ASSIM_GEN = 	parallel_pdaf_mod.F90 \
		init_parallel_pdaf_offline.F90 \
		assim_pdaf_mod.F90 \
		assimilate_pdaf_offline.F90 \
		finalize_pdaf.F90 \
		main_offline.F90 \
		init_ens_cb_pdaf.F90 \
		prepoststep_cb_pdaf.F90 \
		synobs_pdaf_mod.F90 \
		get_obs_cb_pdaf.F90

# Case-specific modules
MOD_ASSIM = 	model_pdaf_mod_offline.F90 \
		statevector_pdaf_mod.F90 \
		initialize_grid.F90 \
		io_pdaf_mod.F90 \

# Interface to PDAF - model sided
SRC_PDAF_INT =  init_pdaf_offline.F90

# User-supplied call-back routines for localization
SRC_CB_LOCAL =  init_n_domains_cb_pdaf.F90 \
		init_dim_l_cb_pdaf.F90

# Routines for PDAF-OMI observation handling
SRC_CB_OMI = 	obs_A_pdafomi.F90 \
		obs_B_pdafomi.F90 \
		callback_obs_pdafomi.F90

# Additional routines call by other routines
SRC_MISC =      init_pdaf_parse.F90

# Full list of user-supplied routines for online modes
SRC_PDAF_USER = $(MOD_ASSIM) $(SRC_CB_OMI) $(SRC_MISC) $(SRC_PDAF_INT) $(SRC_CB) $(SRC_CB_LOCAL)

#######################################################
# object files
OBJ_PDAF_USER := $(SRC_PDAF_USER:%.F90=build/%.o)

# object files
OBJ_PDAF_GEN := $(ASSIM_GEN:%.F90=build/%.o)

######################################################

.phony: info
info:
	@echo "Makefile to build PDAF tutorial offline implementation";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  offline_pdaf  - PDAF offline DA program";
	@echo "  clean         - clean up compilation";
	@echo "  cleandata     - remove output files";
	@echo "  cleanpdaf     - clean up PDAF library compilation";
	@echo "  listarch      - List available choices for PDAF_ARCH";
	@echo "";
	@echo "Path to PDAF (PDAFDIR): " $(PDAFDIR);
	@echo "Machine-specific definitions (PDAF_ARCH): " $(PDAF_ARCH);
	@echo "------------------------------------------------------------------";


offline_pdaf : libpdaf-d.a $(OBJ_PDAF_GEN) $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@ \
	$(OBJ_PDAF_GEN) $(OBJ_PDAF_USER) \
	-L$(PDAFDIR)/lib -lpdaf-d $(LINK_LIBS) $(NC_LIB)

######################################################
# Build PDAF library

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd $(PDAFDIR); make;


######################################################
# use pattern rule to create rules for all object files

$(OBJ_PDAF_USER): build/%.o: src/%.F90
	$(FC) $(OPT) $(CPP_DEFS) $(NC_INC) -I$(PDAFDIR)/include -c -o $@ $(MODULEOPT) build/ $<

$(OBJ_PDAF_GEN): build/%.o: src_generic/%.F90
	$(FC) $(OPT) $(CPP_DEFS) $(NC_INC) -I$(PDAFDIR)/include -c -o $@ $(MODULEOPT) build/ $<

######################################################
# Cleans

.phony: cleanall cleanpdaf clean cleandata
cleanall : cleanpdaf clean cleandata

clean :
	rm -f build/*.o build/*.mod offline_pdaf

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(PDAFDIR); make clean

cleandata:
	rm -f *.nc

######################################################
# List arch files

.phony: listarch
listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 $(PDAFDIR)/make.arch | cut -d"." -f1


######################################################
# Dependencies
# created using 'mkdepends' in PDAF's directory external/
# Run as: ../pdaf/external/mkdepends/mkdepends ./src ./src_generic build
include Depends
