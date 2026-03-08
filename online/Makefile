#######################################################
# Simple Makefile to build the 2D tutotrial model     #
#######################################################

# Include machine-specific definitions
include compile_settings.txt

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

# Adapted for PDAF to separate model and user-code for PDAF
$(OBJ_MODEL): build/%.o: src/%.F90
	$(FC) $(OPT) $(INCLUDES) -c -o $@ $(MODULEOPT) build/ $<

# PDAF: Compile rule for user-code in src_pdaf/
$(OBJ_PDAF_USER): build/%.o: src_pdaf/%.F90
	$(FC) $(OPT) $(INCLUDES) -c -o $@ $(MODULEOPT) build/ $<

######################################################
# Cleans

cleanall : clean cleandata

clean :
	rm -f build/*.o build/*.mod model model_pdaf

cleandata:
	rm -f field*nc


######################################################
# Dependencies
# (created using 'mkdepends' in PDAF's directory external/)

build/main.o: ./src/main.F90 build/model_parallel_mod.o build/model_init_mod.o build/model_step_mod.o build/model_post_mod.o
build/model_init_mod.o: ./src/model_init_mod.F90 build/model_mod.o build/model_parallel_mod.o build/model_io_mod.o
build/model_io_mod.o: ./src/model_io_mod.F90 build/model_parallel_mod.o build/model_mod.o
build/model_mod.o: ./src/model_mod.F90
build/model_parallel_mod.o: ./src/model_parallel_mod.F90
build/model_post_mod.o: ./src/model_post_mod.F90 build/model_mod.o build/model_parallel_mod.o
build/model_step_mod.o: ./src/model_step_mod.F90 build/model_mod.o build/model_parallel_mod.o build/model_io_mod.o
