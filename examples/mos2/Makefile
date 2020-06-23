LTOOL = ../ltool
RM = rm -f

.PHONY: default

default: help

help:
	@echo make targets:
	@echo xdir: transport in x direction
	@echo ydir: transport in y direction
	@echo clean: delete stuff

all: xdir ydir

xdir: clean
	$(LTOOL) -iwinput_x -hwinput_x

ydir: clean
	$(LTOOL) -iwinput_y -hwinput_y

clean:
	$(RM) out/* lattice_dat ph_mat_par *.wbh stump.cmd *.dat H_* *.log
