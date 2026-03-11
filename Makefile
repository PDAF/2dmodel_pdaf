.PHONY: all model online offline generate_covar clean

all: model online offline generate_covar

# definitions for highlighting outputs
bold := $(shell tput bold)
sgr0 := $(shell tput sgr0)

model:
	$(info $(bold)building model$(sgr0))
	make -C 2dmodel model

online:
	$(info $(bold)building online coupled model$(sgr0))
	make -C online model_pdaf

offline:
	$(info $(bold)building offline coupled model$(sgr0))
	make -C offline offline_pdaf

generate_covar:
	$(info $(bold)building generate_covar$(sgr0))
	make -C generate_covar generate_covar

clean:
	make -C 2dmodel/ clean
	make -C online/ clean
	make -C offline/ clean
	make -C generate_covar/ clean
