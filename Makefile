include def.mk

.PHONY: all clean veryclean doc

MAKE = make

all: dbg_libmat rls_libmat dbg_liblat rls_liblat dbg_ltool rls_ltool doc

dbg_libmat:
	$(MAKE) -C libmat debug

rls_libmat:
	$(MAKE) -C libmat release

dbg_liblat:
	$(MAKE) -C liblat debug

rls_liblat:
	$(MAKE) -C liblat release

dbg_ltool:
	$(MAKE) -C ltool debug

rls_ltool:
	$(MAKE) -C ltool release

stat_ltool:
	$(MAKE) -C ltool static

doc:
	$(DOXYGEN) doxyfile


clean:
	$(MAKE) -C libmat clean
	$(MAKE) -C liblat clean
	$(MAKE) -C ltool clean

veryclean:
	$(MAKE) -C libmat veryclean
	$(MAKE) -C liblat veryclean
	$(MAKE) -C ltool veryclean
	$(RM) doc/*
