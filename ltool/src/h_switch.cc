// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "h_switch.h"
#include "ll_io.h"
#include "aux_md5.h"

using namespace aux;
using namespace ll__;

void h_switch(const h_input& inp, std::ostream& os) {
	using namespace aux;

	olf OLF;
	try {
		// try reading user defined layer matrix
		OLF = readOlf(inp.layer_matrix);
	} catch(const std::exception& e) {
		
		// complain if no device length is specified
		if (!inp.device_length) throw(e);

		// check device definition
		if (!inp.device_definition.empty() &&
		   ( inp.device_definition.M()!=DIM__ || (inp.device_definition.N()%2 != 0)))
			throw(std::invalid_argument("bad device definition"));
		for (auto j=inp.device_definition.ccBegin(),
			  e=inp.device_definition.ccEnd(); j!=e; j+=2)
			if (any(j->geq(*(j+1))))
				throw(std::invalid_argument("bad device definition"));

		// read lattice file
		const auto LD = readOlf(OLF__);

		// generate device structure
		{
			if (inp.verbosity && PRINTBIT__)
				os << "Generating Layer matrix...\n";

			// number of unit cells
			const size_t NC = std::ceil(inp.device_length / LD.B(0,0));
			if (inp.verbosity && PRINTBIT__)
				os << "number of unit cells: " << BLUE__
				   << NC << RESET__ << "\n";
			if (NC>1000)
				throw(std::invalid_argument(
				"number of cells is very large: "+std::to_string(NC)));
			
			// copy T, bl, nn
			OLF.T = LD.T; OLF.nn = LD.nn; OLF.bl = LD.bl;

			// new basis
			OLF.B = LD.B; OLF.B(0,0) = LD.B(0,0)*NC;
			
			// new positions and ids
			OLF.Ap = fMat(DIM__,0); OLF.Ap.reserve(LD.N()*NC); OLF.id.reserve(LD.N()*NC);
			for (size_t nc=0; nc!=NC; ++nc) {
				const fMat sh = LD.B.cFront() * nc;
				
				auto j = LD.id.cbegin();
				for (auto i=LD.Ap.ccBegin(),e=LD.Ap.ccEnd(); i!=e; ++i,++j)
					if (inp.device_definition.empty()) {
						OLF.Ap.push_back(*i + sh),
						OLF.id.push_back(*j);
					}
					else {
						// filter according to device definition
						const fMat cont = *i + sh;
						bool keep = false;
						for (auto j=inp.device_definition.ccBegin(),
							  e=inp.device_definition.ccEnd();
							  !keep && j!=e; j+=2)
							keep |= all(j->leq(cont) & (j+1)->geq(cont));
						if (keep)
							OLF.Ap.push_back(cont),
							OLF.id.push_back(*j);
					}
			}
		}


		// write layer matrix to file
		const std::string fileName = inp.layer_matrix.empty() ? "Layer_Matrix.dat":
									inp.layer_matrix;
		printOlf(fileName,OLF);
		if (inp.verbosity && PRINTBIT__)
			os << "file '" << fileName << "' written\n";
	}


	// hr format switch
	switch (fnvHash(inp.hr_format.c_str())) {
		case "omen"_h:
		{
			// adapt L according to r
			fMat L = lm__::diag(OLF.B);
			for (size_t i=1; i!=L.size(); ++i)
				if (inp.r[i]) L[i] = 1e8;

			// call OMEN matrix constructor
			omen::hctor(inp,OLF.Ap,OLF.id,L,os);
		}
		break;
		default:
		{
			// call generic matrix constructor
			h_hctor(inp,OLF.B,OLF.Ap,OLF.id,os);
		}
		break;
	}
}



void h_hctor(const h_input& inp, const fMat& B, const fMat& LM,
		const idv& id, std::ostream& os) {
	assert(B.M() == B.N());
	assert(B.M() == LM.M());
	assert(LM.N() == id.size());

	// read wbhs, set tolerance and equalize energy shifts
	ll_hbondss W(inp,os); W.fixDiagonals(1e-4); W.setQueryTolCartesian(inp.sptol);
	
	// convert ids to type indices, adapt W to structure
	const aTv T = W.ind(id);
	adaptWBH(W,inp,LM,T,os);

	// get interaction radius from inp or derived from wbh
	const double IR = std::isnan(inp.IR) ? W.radius(): inp.IR;
	
	// matrix block size
	const size_t Nw = std::accumulate(T.cbegin(),T.cend(),size_t(0),
		[&W](const size_t s, const aT t) -> size_t
			{ return s + W.Norb(t); });
	if (inp.verbosity & PRINTBIT__) {
		os << "\nlayer matrix of " << LM.N() << " positions loaded\n"
		   << "total number of orbitals in layer matrix: " << Nw << "\n"
		   << "using interaction cutoff radius of " << IR << "\n"
		   << "restricted dimensions are: " << inp.r << "\n\n";
	}
	if (inp.verbosity & WRITEBIT__)
		printPOSCAR(inp.prefix+"lm.psc",B,LM,id,1.0,false);

	// get R grid
	const fMat R = !inp.R.empty() ? inp.R:
		        getConnectedGrid(W,B,LM,T,inp.r,IR,!inp.strict_matching);
	if (inp.verbosity & PRINTBIT__)
		os << "using R grid of "
		   << BLUE__ << R.N() << RESET__ << " points\n";
	if (inp.verbosity & VERBOBIT__)
		os << lm__::T(R).print(0,1) << "\n\n";


	// call hctor using appropriate format writer
	switch (fnvHash(inp.hr_format.c_str())) {
		case "wannier90"_h:
		{
			ll_writerW90 writer(inp.hr_out,R,Nw,inp.pprec);
			hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
			if (inp.verbosity & MD5BIT__)
				os << "md5 sum:" << " \'" << writer.fileName() << "\', "
				   << aux::md5(writer.fileName()) << "\n\n";
		}
		break;
		case "hr32r"_h:
		{
			if (Nw<=std::numeric_limits<uint16_t>::max()) {
				ll_writerBIN<uint16_t,float> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			} else {
				ll_writerBIN<uint32_t,float> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			}
		}
		break;
		case "hr32c"_h:
		{
			if (Nw<=std::numeric_limits<uint16_t>::max()) {
				ll_writerBIN<uint16_t,std::complex<float>> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			} else {
				ll_writerBIN<uint32_t,std::complex<float>> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			}
		}
		break;
		case "hr64r"_h:
		{
			if (Nw<=std::numeric_limits<uint16_t>::max()) {
				ll_writerBIN<uint16_t,double> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			} else {
				ll_writerBIN<uint32_t,double> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			}
		}
		break;
		case "hr64c"_h:
		{
			if (Nw<=std::numeric_limits<uint16_t>::max()) {
				ll_writerBIN<uint16_t,std::complex<double>> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			} else {
				ll_writerBIN<uint32_t,std::complex<double>> writer(inp.hr_out,R,Nw);
				hctor(W,B,LM,T,writer,IR,!inp.strict_matching,inp.Nthreads,inp.verbosity,os);
				if (inp.verbosity & MD5BIT__)
					os << "md5 sum:" << " \'" << writer.fileName() << "\', "
					   << aux::md5(writer.fileName()) << "\n\n";
			}
		}
		break;
		default:
			throw(std::invalid_argument("hamiltonian format '"+inp.hr_format+"' not recognized"));
		break;
	}

	if (!inp.force_check) return;

	os << "check forced!\n";
	auto hr = readHrSparse(inp.hr_out);
	for (auto i=hr.ccBegin(),e=hr.ccBegin()+hr.NR()/2+1; i!=e; ++i) {
		const auto itr = std::lower_bound(e-1,hr.ccEnd(),-*i);
		if (itr==hr.ccEnd() || *itr!=-*i)
			os << "warning: inverse of '" << lm__::T(*i).print(0) << "' not found in the hr\n";
		else
			os << "checking block: '" << lm__::T(*i).print(0) << "'"
			   << (hr[size_t(i)] == ll_sparse(hr[size_t(itr)]).T()
				? GREEN__ "   ok!": RED__ "  bad!") << RESET__ << "\n";
	}
}
