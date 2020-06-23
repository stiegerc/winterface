// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_hio.h"

using namespace lm__;
using namespace ll__;


R_H<ll_sparse> ll__::readHrSparse(const std::string& fileName) {

	// open file
	auto file = aux::openFile<std::ifstream>(fileName, std::ios::in | std::ios::binary);

	// read header
	char head[6];
	file.read(head, 5); head[5] = '\0';

	// read dim, Nw, NR
	fMat R; size_t Nw;
	{
		uint32_t dim, Nw_, NR;
		file.read((char*) &dim, sizeof(uint32_t));
		file.read((char*) &Nw_, sizeof(uint32_t));
		file.read((char*) &NR,  sizeof(uint32_t));
		
		R = fMat(dim,NR); Nw = Nw_;
	}

	// read sparse data
	std::vector<ll_sparse> H; H.reserve(R.N());
	for (auto i=R.cBegin(),e=R.cEnd(); i!=e; ++i) {
		// read R
		file.read((char*) i->data(), R.M()*sizeof(double));

		// read nnz
		uint64_t nnz;
		file.read((char*) &nnz, sizeof(uint64_t));
		H.push_back(ll_sparse(Nw,Nw)); H.back().vec_.reserve(nnz);

		switch (aux::fnvHash(head)) {
		case "hr32r"_h:
		{
			typedef ll_sphelBIN<float> sphel;
			sphel cel(0,0,.0);
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &cel, sizeof(sphel));
				H.back().vec_.push_back({size_t(cel.m),size_t(cel.n),cel.h});
			}
		}
		break;
		case "hr32c"_h:
		{
			typedef ll_sphelBIN<std::complex<float>> sphel;
			sphel cel(0,0,.0);
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &cel, sizeof(sphel));
				H.back().vec_.push_back({size_t(cel.m),size_t(cel.n),cel.h});
			}
		}
		break;
		case "hr64r"_h:
		{
			typedef ll_sphelBIN<double> sphel;
			sphel cel(0,0,.0);
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &cel, sizeof(sphel));
				H.back().vec_.push_back({size_t(cel.m),size_t(cel.n),cel.h});
			}
		}
		break;
		case "hr64c"_h:
		{
			typedef ll_sphelBIN<std::complex<double>> sphel;
			sphel cel(0,0,.0);
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &cel, sizeof(sphel));
				H.back().vec_.push_back({size_t(cel.m),size_t(cel.n),cel.h});
			}
		}
		break;
		default:
		throw(std::runtime_error("bad header in file '"+fileName+"'"));
		}

		// sort entries in sparse
		std::sort(H.back().vec_.begin(),H.back().vec_.end());
	}

	return {std::move(R),std::move(H)};
}

void ll__::adaptWBH(ll_hbondss& W, const ll_hbondss_input& inp,
			const fMat& Ap, const aTv& T, std::ostream& os) {

	assert(W.dim()==Ap.M());
	assert(Ap.N()==T.size());
	if (!inp.perimeter_radius) return;

	// get connections	
	const auto B = W.getPerimeterConnections(Ap,T,inp,os);
	const auto S = W.getSubstitutes(Ap,T,inp.perimeter_radius);
	
	// equalize energy shifts
	W.equalizeShifts(B,inp,os);

	// equalize wannier signs
	W.equalizeSigns(B,S,inp,os);

	// dump adapted wbh
	if (inp.wbh_out.size() == W.size())
		for (size_t i=0; i!=W.size(); ++i) {
			(W.cbegin()+i)->writeToFile(inp.prefix+inp.wbh_out[i]);
			if (inp.verbosity & VERBOBIT__)
				os << "file '"+inp.prefix+inp.wbh_out[i]+"' written"
				   << (inp.verbosity & MD5BIT__ ?
					", "+aux::md5(inp.prefix+inp.wbh_out[i]): "")
				   << "\n";
		}
}




void ll__::h_gen_hctor(const h_input& inp, const fMat& B, const fMat& LM,
		const idv& id, std::ostream& os) {
	assert(B.M() == B.N());
	assert(B.M() == LM.M());
	assert(LM.N() == id.size());

	// read wbhs, set tolerance and equalize energy shifts
	ll_hbondss W(inp,os); W.fixDiagonals(1e-4);
	W.setQueryTolCartesian(inp.sptol);
	if (inp.perimeter_radius) {
		// get connections	
		const auto T = W.ind(id);
		const auto B = W.getPerimeterConnections(LM,T,inp,os);
		const auto S = W.getSubstitutes(LM,T,inp.perimeter_radius);
		
		// equalize energy shifts
		W.equalizeShifts(B,inp,os);

		// equalize wannier signs
		W.equalizeSigns(B,S,inp,os);
	
		// dump adapted wbh
		if (inp.wbh_out.size() == W.size())
			for (size_t i=0; i!=W.size(); ++i) {
				(W.cbegin()+i)->writeToFile(inp.prefix+inp.wbh_out[i]);
				if (inp.verbosity & VERBOBIT__)
					os << "file '"+inp.prefix+inp.wbh_out[i]+"' written"
					   << (inp.verbosity & MD5BIT__ ?
						", "+aux::md5(inp.prefix+inp.wbh_out[i]): "")
					   << "\n";
			}
	}
	
	// convert ids to type indices
	const aTv T = W.ind(id);

	// get interaction radius from inp or derived from wbh
	const double IR = std::isnan(inp.IR) ? W.radius(): inp.IR;
	
	// matrix block size
	const size_t Nw = std::accumulate(T.cbegin(),T.cend(),size_t(0),
		[&W](const size_t s, const aT t) -> size_t
			{ return s + W.Norb(t); });
	if (inp.verbosity & PRINTBIT__) {
		os << "\nlayer matrix of " << LM.N() << " positions loaded\n"
		   << "total number of orbitals in layer matrix: " << Nw << "\n"
		   << "using relative tolerance of " << mtol() << "\n"
		   << "using interaction cutoff radius of " << IR << "\n"
		   << "restricted dimensions are: " << inp.r << "\n\n";
	}


	// get R grid
	fMat R(W.dim(),0);
	if (!inp.R.empty()) {
		R = inp.R;
		if (inp.verbosity & PRINTBIT__)
			os << "using user defined R grid of "
			   << BLUE__ << R.N() << RESET__ << " points\n";
		if (inp.verbosity & VERBOBIT__)
			os << lm__::T(R).print(0,1) << "\n\n";
	} else {

		// probe interactions lambda
		const auto probe = [&W,IR,&LM,&T,strict=inp.strict_matching]
				(const fMat& sh) -> bool {
			
			for (auto i1=LM.ccBegin(),e=LM.ccEnd(); i1!=e; ++i1) {
				const aT t1 = T[(size_t)i1];
				
				for (auto i2=LM.ccBegin(); i2!=e; ++i2) {
					const aT t2 = T[(size_t)i2];
					
					// check bond length
					const fMat b = *i2 + sh -*i1;
					if (norm(b)>IR) continue;
					
					// try direct
					if (!W.getInteraction({t1,t2},b).empty())
						return true;
					
					// try approximate
					if (!strict && ! W.getApproximateInteraction({t1,t2},b).H.empty())
						return true;
				}
			}
			return false;
		};


		// recursive scan space lambda
		const size_t N = std::ceil(std::pow((2.0*IR),Nfalse(inp.r)) / std::abs(lm__::det(B)));
		const fMat NN = genNNmat(inp.r);
		R.reserve(N);
		
		std::function<void(const fMat&)> scanner;
		scanner = [&scanner,&probe,&R,&NN,&B]
				(const fMat& cR) -> void {

			// probe this block
			if (!probe(B.prod(cR))) return;

			// insert this block
			R.cInsert(std::lower_bound(R.ccBegin(),R.ccEnd(),cR),cR);
			
			// recursive calls for all neighbor cR
			for (auto i=NN.ccBegin(),e=NN.ccEnd(); i!=e; ++i) {
				const fMat nR = cR + *i;
				const auto itr = std::lower_bound(R.ccBegin(),R.ccEnd(),nR);
				if (itr==R.ccEnd() || *itr!=nR)
					scanner(nR);
			}
		};
		
		// call scanner at origin
		if (inp.verbosity & PRINTBIT__)
			os << "\nprobing grid...";
		scanner(lm__::zeros<fMat>(W.dim(),1));
		if (inp.verbosity & PRINTBIT__)
			os << " done, found " << BLUE__ << R.N() << RESET__ " grid points\n";
		if (inp.verbosity & VERBOBIT__)
			os << lm__::T(R).print(0,1) << "\n\n";
	}
	

	// write Hamiltonian matrices lambda
	const auto wHm = [&LM,&R,&B,&T,IR,&W,&inp,&os](auto& cont) -> void {
		if (inp.verbosity & PRINTBIT__)
			os << "writing Hamiltonian data to file '" << cont.fileName() << "'\n";
		
		// generate Hamiltonians for all entries in the grid
		size_t nnz = 0;
		for (auto r=R.ccBegin(),re=R.ccEnd(); r!=re; ++r) {
			if (inp.verbosity & PRINTBIT__)
				os << "writing block: " << BLUE__ << lm__::T(*r).print(0)
				   << RESET__ << " ..." << std::flush;
			
			writeHamBlock(LM,T,LM,T,W,!inp.strict_matching,cont,inp.Nthreads,IR,
				[sh=B.prod(*r)](const fMat& p1, const fMat& p2)->fMat{return (p2+sh)-p1;});
			cont.flush();
			
			if (inp.verbosity & PRINTBIT__)
				os << " done!, nnz: " << (cont.nnz()-nnz) << "\n" << std::flush;
			nnz = cont.nnz();
		}
	};


	// generate output container and write data
	switch (aux::fnvHash(inp.hr_format.c_str())) {
		case "wannier90"_h:
		{
			ll_w90_hr cont(inp.hr_out,R,Nw,inp.pprec);
			wHm(cont);
		}
		break;
		case "bin32re"_h:
		{
			ll_bin_hr<float> cont(inp.hr_out,R,Nw);
			wHm(cont);
		}
		break;
		case "bin32cpx"_h:
		{
			ll_bin_hr<std::complex<float>> cont(inp.hr_out,R,Nw);
			wHm(cont);
		}
		break;
		case "bin64re"_h:
		{
			ll_bin_hr<double> cont(inp.hr_out,R,Nw);
			wHm(cont);
		}
		break;
		case "bin64cpx"_h:
		{
			ll_bin_hr<std::complex<double>> cont(inp.hr_out,R,Nw);
			wHm(cont);
		}
		break;
		default:
			throw(std::invalid_argument("format '"+inp.hr_format+"' not recognized"));
		break;
	}
}
