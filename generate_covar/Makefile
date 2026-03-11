#######################################################
# Simple Makefile to build the generate_covar program #
#######################################################

# Include machine-specific definitions
include ../compile_settings.txt

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
	@echo  "  clean          - clean up compilation";
	@echo  "-------------------------------------------------------------------------------";

generate_covar: libpdaf-d.a $(OBJ_COVAR)
	rm -f $@
	$(LD)  $(FFLAGS) -o $@ $(OBJ_COVAR) $(LIBS)

######################################################
# Build PDAF library

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	make -C $(PDAFDIR)

######################################################
# use pattern rule to create rules for all object files

build/%.o: src/%.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $(MODULEOPT) build/ $<

######################################################
# Cleans

.phony: clean
clean :
	rm -f build/*.o build/*.mod generate_covar
