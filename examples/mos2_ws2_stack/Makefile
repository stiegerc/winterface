LTOOL = ../ltool
RM = rm -f

.PHONY: default

default: help

help:
	@echo make targets:
	@echo xdir: transport in x direction
	@echo clean: delete stuff

xdir: clean
	$(LTOOL) -i -h

clean:
	$(RM) out/* lattice_dat ph_mat_par *.wbh stump.cmd *.dat H_* *.log
