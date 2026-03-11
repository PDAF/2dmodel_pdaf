#######################################################
# Simple Makefile to build the 2D tutotrial model     #
# offline                                             #
#######################################################

# Include machine-specific definitions
include ../compile_settings.txt

#######################################################
# Sourcefiles

# PDAF: Include source files for PDAF coupling
include srcfiles_pdaf

#######################################################
# object files
OBJ_PDAF_USER := $(SRC_PDAF_USER:%.F90=build/%.o)

######################################################

.phony: info
info:
	@echo "Makefile to build advanced 2d tutorial model of PDAF";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  offline_pdaf    - assimilation model";
	@echo "  clean         - clean up compilation";
	@echo "  cleandata     - remove output files";
	@echo "------------------------------------------------------------------";


offline_pdaf : $(OBJ_PDAF_USER)
	$(LD)  $(FFLAGS)  -o $@ \
	$(OBJ_MODEL) $(OBJ_PDAF_USER) $(LIBS)

######################################################
# use pattern rule to create rules for all object files

$(OBJ_PDAF_USER): build/%.o: src/%.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@  $(MODULEOPT) build/ $<

######################################################
# Cleans

.phony: cleanall
cleanall : cleanpdaf clean cleandata

.phony: clean
clean :
	rm -f build/*.o build/*.mod offline_pdaf

.phony: cleandata
cleandata:
	rm -f *.nc

######################################################
# Dependencies
# (created using 'mkdepends' in PDAF's directory external/)
include Depends
