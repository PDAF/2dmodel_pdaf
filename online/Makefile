#######################################################
# Simple Makefile to build the 2D tutotrial model     #
#######################################################

# Include machine-specific definitions
include compile_settings.txt

# this makefile is intendet to build the model with
# online coupling. When compiling the model without
# PDAF (target `model`) one has to disable -DUSE_PDAF
FFLAGS :=-DUSE_PDAF

#######################################################
# Sourcefiles

SRC_MODEL =  	model_parallel_mod.F90 \
		model_mod.F90 \
		model_io_mod.F90 \
		model_init_mod.F90 \
		model_step_mod.F90 \
		model_post_mod.F90 \
		main.F90

# PDAF: Include source files for PDAF coupling
include srcfiles_pdaf

#######################################################
# Define object files - to be create in build/
OBJ_MODEL := $(SRC_MODEL:%.F90=build/%.o)

# PDAF: Define PDAF-related object files
OBJ_PDAF_USER := $(SRC_PDAF_USER:%.F90=build/%.o)

######################################################

.phony: info
info:
	@echo "Makefile to build advanced 2d tutorial model of PDAF";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  model         - the 2D tutorial model  (use without -DUSE_PDAF)";
	@echo "  model_pdaf    - assimilation model      (use with -DUSE_PDAF)";
	@echo "  clean         - clean up compilation";
	@echo "  cleandata     - remove output files";
	@echo "------------------------------------------------------------------";


model : $(OBJ_MODEL)
	$(LD)  $(OPT_LNK)  -o $@  $(OBJ_MODEL) $(LIBS)

# PDAF: Target to compile assimilative model
model_pdaf : $(OBJ_MODEL) $(OBJ_PDAF_USER)
	$(LD)  $(OPT_LNK)  -o $@ \
	$(OBJ_MODEL) $(OBJ_PDAF_USER) $(LIBS)

######################################################
# use pattern rule to create rules for all object files

$(OBJ_MODEL): build/%.o: src/%.F90
	$(FC) $(FFLAGS) $(OPT) $(INCLUDES) -c -o $@ $(MODULEOPT) build/ $<

# PDAF: Compile rule for user-code in src_pdaf/
$(OBJ_PDAF_USER): build/%.o: src_pdaf/%.F90
	$(FC) $(FFLAGS) $(OPT) $(INCLUDES) -c -o $@ $(MODULEOPT) build/ $<

######################################################
# Cleans

.phony: cleanall
cleanall : clean cleandata

.phony: clean
clean :
	rm -f build/*.o build/*.mod model model_pdaf

.phony: cleandata
cleandata:
	rm -f field*nc

######################################################
# Dependencies
# (created using 'mkdepends' in PDAF's directory external/)
include Depends
