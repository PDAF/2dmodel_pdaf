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
# make offline_pdaf PDAF_ARCH=linux_gfortran_openmpi PDAFDIR=../PDAF
#####################################################################

# Include machine-specific definitions
# For available include files see PDAF directory make.arch
# To choose a file, set PDAF_ARCH for make or by an
# environment variable.
include $(PDAFDIR)/make.arch/$(PDAF_ARCH).h


######################################################

# Modules used for PDAF
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
OBJ_PDAF_USER := $(SRC_PDAF_USER:%.F90=build/%.o)

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


offline_pdaf : libpdaf-d.a $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@ \
	$(OBJ_PDAF_USER) \
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

######################################################
# Cleans

.phony: cleanall
cleanall : cleanpdaf clean cleandata

.phony: clean
clean :
	rm -f build/*.o build/*.mod offline_pdaf

.phony: cleanpdaf
cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd $(PDAFDIR); make clean

.phony: cleandata
cleandata:
	rm -f *.nc

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1
