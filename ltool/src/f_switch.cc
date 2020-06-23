// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "f_switch.h"
#include "aux_parser.h"
#include "ll_io.h"
#include <fstream>
#include <string>

using namespace aux;
using namespace ll__;


void f_switch(const f_input& inp, std::ostream& os) {

	if (!inp.wout.empty()) {
		// read from wout
		const auto B = readB(inp.wout);
		const auto Ap = readAp(inp.wout,false);
		const auto Wp = readWp(inp.wout);

		// write psc
		printPOSCAR(inp.prefix+inp.pscout+".ap",B,inp.direct ? B.leftDivide(Ap.Ap()): Ap.Ap(),
						Ap.id(),1.0,inp.direct);
		printPOSCAR(inp.prefix+inp.pscout+".wp",B,inp.direct ? B.leftDivide(Wp.Wp()): Wp.Wp(),
						idv{Wp.N(),"H"},1.0,inp.direct);
	}

	if (!inp.lattice_dat.empty()) {
		// read from lattice_dat
		auto D = readOlf(inp.lattice_dat);

		// strip ids
		if (inp.strip)
			for (auto& s: D.id) s = ll_cell::stripId(s);

		// write psc
		printPOSCAR(inp.prefix+inp.pscout,ll_cell(D.B,D.B.leftDivide(D.Ap),D.id),inp.direct);
	}
}
