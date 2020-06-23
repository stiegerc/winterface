LTOOL = ../ltool
RM = rm -f
MATLAB = matlab

.PHONY: default

default: help

help:
	@echo make targets:
	@echo xdir: transport in x direction
	@echo bs: compute local bandstructures
	@echo clean: delete stuff

xdir: clean
	$(LTOOL) -i
	$(MATLAB) -nodisplay -nosplash -nodesktop -r "run('./make_layer_matrix.m');exit;"
	$(LTOOL) -h

bs:
	$(LTOOL) -bbinput

clean:
	$(RM) out/* lattice_dat ph_mat_par *.wbh stump.cmd *.dat H_* *.log bsout/*
