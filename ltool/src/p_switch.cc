// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "p_switch.h"
#include "ll_io.h"

using namespace lm__;

void p_switch(const p_input& inp, std::ostream& os) {
	
	ll_cell cell(inp.pscin);	// initial cell
	
	// make cell primitive
	set_mtol(inp.tol);
	cell.makePrimitive(inp.r);
	reset_mtol();
	
	// write result POSCAR
	printPOSCAR(inp.pscout,cell,inp.direct,inp.strip);
}
