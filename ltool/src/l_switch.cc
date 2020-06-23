// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "l_switch.h"
#include "ll_io.h"
#include "ll_fn.h"


using namespace ll__;

void l_switch(const l_input& inp, std::ostream& os) {
	
	if (det(inp.ROT)!=1.0)
		throw(std::invalid_argument("parameter ROT invalid"));
	if (inp.C!=round(inp.C) || std::abs(det(inp.C))<1.0)
		throw(std::invalid_argument("parameter C invalid"));
	if (inp.R!=round(inp.R))
		throw(std::invalid_argument("parameter R invalid"));
	
	// read raw poscar
	ll__::psc P;
	if (inp.wbh.empty())
		P = readPOSCAR(inp.pscin);
	else {
		auto cell = extractCell(inp.wbh);
		P.N = cell.Ntype();
		P.B = cell.moveB();
		P.Ap = cell.moveAp();
		P.id = cell.moveId();
	}

	// rotate basis
	if (inp.ROT != lm__::eye<fMat>(DIM__,DIM__) ||
	    inp.phi_x || inp.phi_y || inp.phi_z) {
		const double phi_x = (M_PI/180.0) * inp.phi_x;
		const double phi_y = (M_PI/180.0) * inp.phi_y;
		const double phi_z = (M_PI/180.0) * inp.phi_z;

		P.B = 
		fMat({             1.0,		    0.0,             0.0,
		                   0.0, std::cos(phi_x), std::sin(phi_x),
		                   0.0,-std::sin(phi_x), std::cos(phi_x)}
			,DIM__,DIM__).prod(
		fMat({ std::cos(phi_y),             0.0,-std::sin(phi_y),
		                   0.0,             1.0,             0.0,
		       std::sin(phi_y),             0.0, std::cos(phi_y)}
			,DIM__,DIM__).prod(
		fMat({ std::cos(phi_z), std::sin(phi_z),             0.0,
		      -std::sin(phi_z), std::cos(phi_z),             0.0,
				   0.0,             0.0,             1.0}
			,DIM__,DIM__).prod(
		inp.ROT.prod(P.B))));
	}
	
	// strip id
	if (inp.strip)
		for (auto& i: P.id) i = ll_cell::stripId(i);

	// adjust vacuum	
	if (!inp.vac.empty())
	for (size_t i=0; i!=P.B.M(); ++i) {
		if (!inp.r[i]) continue;

		// vacuum in terms of current basis
		const double f = (inp.vac.size()==P.Ap.M() ?
				  inp.vac[i]: inp.vac.front())/norm(P.B.cAt(i));

		// center and check vacuum validity
		auto r = P.Ap.rBegin()+i;
		*r %= 1.0; *r += .5-com(*r).front(); *r %= 1.0;
		if (max(*r)-min(*r)>f)
			throw(std::invalid_argument("insufficient vacuum"));

		// scale, recenter
		*r /= f, *r += .5-.5/f;

		// rescale NB
		P.B.cAt(i) *= f;
	}

	// expand basis
	if (inp.C != lm__::eye<fMat>(DIM__,DIM__)) {
		set_mtol(WTOL__);
		ll_cell cell(std::move(P.B),P.Ap%1.0,P.N,P.id); cell.expand(inp.C);
		reset_mtol();
		P.N = cell.Ntype(); P.B = cell.moveB(); P.Ap = cell.moveAp(); P.id = cell.id();
	}
	
	// add bond centers
	if (inp.bond_factor) {
		ll_cell cell(std::move(P.B),P.Ap%1.0,P.N,P.id);
		cell.merge(cell.getBonds(inp.bond_factor,genNNmat(!inp.r)).getBondCenters());
		P.N = cell.Ntype(); P.B = cell.moveB(); P.Ap = cell.moveAp(); P.id = cell.id();
	}
	
	// generate T vector
	aTv T; T.reserve(std::accumulate(P.N.cbegin(),P.N.cend(),size_t(0)));
	for (size_t n=0; n!=P.N.size(); ++n)
		for (size_t i=0; i!=P.N[n]; ++i) T.push_back(n);

	// get R images
	fMat NAp(P.Ap.M(),0); NAp.reserve(inp.R.N()*P.Ap.N());
	aTv NT; NT.reserve(NAp.ccap());
	
	for (auto i=inp.R.ccBegin(),ie=inp.R.ccEnd(); i!=ie; ++i) {
		auto tmp = P.Ap; cadd(tmp,*i);
		NAp.push_back(tmp);
		NT.insert(NT.end(),T.cbegin(),T.cend());
	}
	
	printPOSCAR(inp.prefix+inp.pscout,P.B,inp.direct ? NAp: P.B.prod(NAp),
			NT,P.id,1.0,inp.direct);
}
