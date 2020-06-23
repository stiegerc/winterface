LTOOL = ../ltool
RM = rm -f
MATLAB = matlab

.PHONY: default

default: help

help:
	@echo make targets:
	@echo xdir: transport in x direction
	@echo clean: delete stuff

xdir: clean
	./make_wbhs
	$(MATLAB) -nodisplay -nosplash -nodesktop -r "run('./mod_wbhs.m');exit;"
	$(LTOOL) -h

clean:
	$(RM) out/* wbhs/* H_* *.log
