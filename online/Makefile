# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

######################################################
# definitions for highlighting outputs
bold := $(shell tput bold)
sgr0 := $(shell tput sgr0)

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
$(info $(bold)use machine-specific definitions in $(PDAFDIR)/make.arch/$(PDAF_ARCH).h$(sgr0))

SRCDIR:=src
OBJDIR:=build
INCDIR:=include
BINDIR:=bin

# End of user specifications
######################################################

# Modules used for the model part
SRC_MODEL =  	model_parallel_mod.F90 \
		model_mod.F90 \
		model_io_mod.F90 \
		model_init_mod.F90 \
		model_step_mod.F90 \
		model_post_mod.F90 \
		main.F90

# Moduls used for PDAF
MOD_ASSIM = 	model_pdaf_mod.F90 \
		parallel_pdaf_mod.F90 \
		assimilation_pdaf_mod.F90 \
		statevector_pdaf_mod.F90 \
		io_pdaf_mod.F90 \
		parser_mpi.F90

# Interface to PDAF - model sided
SRC_PDAF_INT =  init_parallel_pdaf.F90 \
		init_pdaf.F90 \
		assimilate_pdaf.F90 \
		finalize_pdaf.F90

# user-supplied call-back routines
SRC_CB  =       init_ens_pdaf.F90 \
		prepoststep_pdaf.F90 \
		next_observation_pdaf.F90 \
		distribute_state_pdaf.F90 \
		collect_state_pdaf.F90 \
		get_obs_pdaf.F90

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
OBJ_MODEL := $(SRC_MODEL:%.F90=$(OBJDIR)/%.o)
OBJ_PDAF_USER := $(SRC_PDAF_USER:%.F90=$(OBJDIR)/%.o)

######################################################

info:
	@echo "Makefile to build PDAF tutorial online implementation";
	@echo "Example: 2D model with parallelization";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  model         - just the forward model  (use without -DUSE_PDAF)";
	@echo "  model_pdaf    - assimilation model      (use with -DUSE_PDAF)";
	@echo "  (if compiled without -DUSE_PDAF, both programs behave equally)"
	@echo "";
	@echo "Depending on whether PDAF should be active, one has to define";
	@echo "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo  "List available settings for PDAF_ARCH with 'make listarch'";
	@echo "Path to PDAF (PDAFDIR): " $(PDAFDIR);
	@echo "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo "------------------------------------------------------------------";

#######################################################
# compiler instructions

INCLUDE:= $(MPI_INC)\
          $(NC_INC) \
          $(MPI_INC) \
          -I$(PDAFDIR)/include

COMPILE.f90 = $(FC) $(OPT) $(CPP_DEFS) $(INCLUDE) -c -o $@ $(MODULEOPT) $(INCDIR)

#######################################################

# requires compiling without -DUSE_PDAF
.PHONY: model
model: $(BINDIR)/model

# requires compiling with -DUSE_PDAF
.PHONY: model_pdaf
model_pdaf: $(BINDIR)/model_pdaf

######################################################

.PHONY:
libpdaf-d.a:
	$(info $(bold)Generate Filter library$(sgr0))
	make all -C $(PDAFDIR)

######################################################

# use pattern rule to create rules for all object files

$(OBJDIR)/%.o: $(SRCDIR)/%.F90
	$(info $(bold)compile $<$(sgr0))
	$(COMPILE.f90) $<

$(BINDIR)/model: $(OBJ_MODEL) | directories
	$(info $(bold)linking $@$(sgr0))
	$(LD) $(OPT_LNK) -o $@ $(OBJ_MODEL) $(LINK_LIBS) $(NC_LIB) $(LDFLAGS)
	chmod u=rwx $@

PROG_DEPS = $(OBJ_MODEL) $(OBJ_PDAF_USER)

$(BINDIR)/model_pdaf: $(PROG_DEPS) | libpdaf-d.a directories
	$(info $(bold)linking $@$(sgr0))
	$(LD) $(OPT_LNK) -o $@ $(PROG_DEPS) -L$(PDAFDIR)/lib -lpdaf-d $(LINK_LIBS) $(NC_LIB) $(LDFLAGS)
	chmod u=rwx $@

#######################################################
MISSINGDIRS:= $(INCDIR) $(OBJDIR) $(LIBDIR) $(BINDIR)
.PHONY: directories
directories: $(MISSINGDIRS)

$(MISSINGDIRS):
	mkdir -p $@

######################################################
# Cleans

.PHONY: cleanall
cleanall : cleanpdaf clean cleandata

.PHONY: clean
clean :
	$(info $(bold)clean$(sgr0))
	rm -rf $(OBJDIR)/*.o
	rm -rf $(INCDIR)/*.mod
	rm -rf $(BINDIR)/model
	rm -rf $(BINDIR)/model_pdaf

.PHONY: cleanpdaf
cleanpdaf:
	$(info $(bold)Clean up PDAF directory$(sgr0))
	make clean -C $(PDAFDIR)

.PHONY: cleandata
cleandata:
	rm -f *.nc

.PHONY: cleandataq
cleandataq:
	@rm -f ens*.txt state*.txt true*.txt variance*txt

######################################################
# List arch files

.PHONY: listarch
listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1

include Depends
