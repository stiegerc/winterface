LTOOL = ../ltool
RM = rm -f

.PHONY: default

default: help

help:
	@echo make targets:
	@echo b: compute bandstructure for MoS2
	@echo bl: compute bandstructure in legacy mode
	@echo '    'requires to run 'xdir' in ../mos2 first
	@echo l: manipulate lattices
	@echo p: find primitive cell from super cell
	@echo w: generate wannier bonds hamiltonian from MoS2
	@echo w: generate Wannier bonds Hamiltonian
	@echo clean: delete stuff

all: b bl hx l p w

b:
	$(LTOOL) -binput_b

bl:
	$(LTOOL) -binput_bl

hx:
	$(RM) Layer_Matrix.dat
	$(LTOOL) -iinput_hx1 -hinput_hx1 -binput_hx1b
	$(RM) Layer_Matrix.dat
	$(LTOOL) -iinput_hx2 -hinput_hx2 -binput_hx2b

l:
	$(LTOOL) -finput_lf -linput_l

p:
	$(LTOOL) -finput_pf -linput_pl -pinput_p

w:
	$(LTOOL) -winput_w


clean:
	$(RM) lattice_dat ph_mat_par *.log *.hr stump.cmd *.wad Layer_Matrix.dat
	$(RM) output/b/* output/bl/* output/hx/* output/l/* output/p/* output/w/*
