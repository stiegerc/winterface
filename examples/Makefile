.PHONY: all clean

all:
	make -C mos2 xdir
	make -C ws2 xdir
	make -C mos2_ws2_stack xdir
	make -C ws2_mos2_lateral xdir
	make -C ws2_mos2_lateral bs
	make -C mos2_ws2_multiple xdir

clean:
	make -C mos2 clean
	make -C ws2 clean
	make -C mos2_ws2_stack clean
	make -C ws2_mos2_lateral clean
	make -C ws2_mos2_lateral clean
	make -C mos2_ws2_multiple clean
