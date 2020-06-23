// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_hio.h"
#include "aux_md5.h"

using namespace ll__;
using namespace lm__;


// probe R grid
template <class WT>
fMat ll__::getConnectedGrid(const WT& W, const fMat& B, const fMat& Ap, const aTv& T,
			const rv& r, const double IR, const bool approx) noexcept {
	assert(W.dim()==B.M());
	assert(B.M()==B.N());
	assert(B.M()==Ap.M());
	assert(Ap.N()==T.size());
	assert(r.size()==W.dim());

	// probe interactions lambda
	const auto probe = [&W,IR,&Ap,&T,approx]
			(const lm__::fMat& sh) -> bool {
		
		for (auto i1=Ap.ccBegin(),e=Ap.ccEnd(); i1!=e; ++i1) {
			const aT t1 = T[(size_t)i1];
			
			for (auto i2=Ap.ccBegin(); i2!=e; ++i2) {
				const aT t2 = T[(size_t)i2];
				
				// check bond length
				const fMat b = *i2 + sh -*i1;
				if (norm(b)>IR) continue;
				
				// try direct
				if (!W.getInteraction({t1,t2},b).empty())
					return true;
				
				// try approximate
				if (approx && ! W.getApproximateInteraction({t1,t2},b).H.empty())
					return true;
			}
		}
		return false;
	};
		
	const fMat NN = genNNmat(r);
	fMat res(Ap.M(),0); res.reserve(100);
	
	// recursive scan space lambda
	std::function<void(const fMat&)> scanner;
	scanner = [&scanner,&probe,&res,&NN,&B]
			(const fMat& cR) -> void {

		// probe this block
		if (!probe(B.prod(cR))) return;
		
		// insert this block
		res.cInsert(std::lower_bound(res.ccBegin(),res.ccEnd(),cR),cR);
		
		// recursive calls for all neighbor cR
		for (auto i=NN.ccBegin(),e=NN.ccEnd(); i!=e; ++i) {
			const fMat nR = cR + *i;
			const auto itr = std::lower_bound(res.ccBegin(),res.ccEnd(),nR);
			if (itr==res.ccEnd() || *itr!=nR)
				scanner(nR);
		}
	};

	// call scanner at origin
	scanner(zeros<fMat>(W.dim(),1));
	
	// shrink and return
	res.shrink_to_fit();
	return res;
}
template fMat ll__::getConnectedGrid(const ll_hbonds& W, const fMat& B, const fMat& Ap, const aTv& T,
			const rv& r, const double IR, const bool approx) noexcept;
template fMat ll__::getConnectedGrid(const ll_hbondss& W, const fMat& B, const fMat& Ap, const aTv& T,
			const rv& r, const double IR, const bool approx) noexcept;


// hamiltonian constructor
template <class WT>
void ll__::hctor(const WT& W, const fMat& B, const fMat& Ap, const aTv& T,
		  ll_writer& writer, const double IR, const bool approx,
		  const size_t Nthreads, const size_t verbosity, std::ostream& os) {
	
	assert(W.dim()==Ap.M());
	assert(Ap.N()==T.size());

	// lambda to get all interactions for one atom with all the others
	auto getInteractions = [&Ap,IR,approx,&T,&W]
			(const fMat& p1, const aT t1, const size_t m) -> std::vector<sphel> {
		
		size_t n=0;
		std::vector<sphel> buff; buff.reserve(Ap.N());
		
		auto i2 = Ap.ccBegin();
		for (const auto t2: T) {
			const auto b = *i2 - p1;
			if (norm(b)<=IR) {
			
				// find H, first try direct, then approximate
				auto H = W.getInteraction({t1,t2},b);
				if (H.empty() && approx)
					H = W.getApproximateInteraction({t1,t2},b).H;

				// enter data into buffer
				if (!H.empty()) {
					for (size_t hn=0; hn!=H.N(); ++hn)
					for (size_t hm=0; hm!=H.M(); ++hm)
					if (H(hm,hn)!=0.0)
						buff.push_back({m+hm,n+hn,H(hm,hn)});
				}
			}
			n+=W.Norb(t2); ++i2;
		}
		// sort buffer and return
		std::sort(buff.begin(),buff.end());
		return buff;
	};
	
	// write blocks
	while (!writer.eof()) {

		if (verbosity & PRINTBIT__)
			os << "writing block: \'" << writer.c_descr() << "\'..." << std::flush;

		// initiate block
		writer.newBlock();

		// spawn tolerance stacks for parallel section
		lm__::spawn_mtol_stacks(Nthreads ? Nthreads: omp_get_max_threads());
		
		// scan interactions between atoms and write data to container
		const auto e = Ap.ccEnd();
		#pragma omp parallel for ordered schedule(static,1) num_threads(Nthreads)
		for (auto i1=Ap.ccBegin(); i1<e; ++i1) {
			const auto buff = getInteractions(*i1-B.prod(writer.c_R()),T[size_t(i1)],
				std::accumulate(T.cbegin(),T.cbegin()+size_t(i1),size_t(0),
				[&W](const size_t s, const aT t){return s+W.Norb(t);}));
			#pragma omp ordered
			writer.insert(buff);
		}
		
		// collapse tolerance stacks back to master
		lm__::collapse_mtol_stacks();

		// report results
		if (verbosity & PRINTBIT__)
			os << "   done!\n" << "nnz is " << writer.c_nnz()
			   << ", size in MB: " << (writer.c_bytes()*1e-6) << "\n\n" << std::flush;
		
		// finalize block
		writer.flush();
	}

	if (verbosity & PRINTBIT__)
		os << "all done! nnz check " << (writer.nnzTest() ?
				GREEN__ "passed": RED__   "failed") << RESET__ << "\n";
}
template void ll__::hctor(const ll_hbonds& W, const fMat& B, const fMat& Ap, const aTv& T,
		  ll_writer& writer, const double IR, const bool approx,
		  const size_t Nthreads, const size_t verbosity, std::ostream& os);
template void ll__::hctor(const ll_hbondss& W, const fMat& B, const fMat& Ap, const aTv& T,
		  ll_writer& writer, const double IR, const bool approx,
		  const size_t Nthreads, const size_t verbosity, std::ostream& os);


// read hr sparse files as written by ll_writerBIN
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
		case "hrssr"_h:
		{
			typedef uint16_t IT; typedef float FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
			}
		}
		break;
		case "hrlsr"_h:
		{
			typedef uint32_t IT; typedef float FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
			}
		}
		break;
		case "hrsdr"_h:
		{
			typedef uint16_t IT; typedef double FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
			}
		}
		break;
		case "hrldr"_h:
		{
			typedef uint32_t IT; typedef double FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
			}
		}
		break;
		case "hrssc"_h:
		{
			typedef uint16_t IT; typedef std::complex<float> FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
			}
		}
		break;
		case "hrlsc"_h:
		{
			typedef uint32_t IT; typedef std::complex<float> FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
			}
		}
		break;
		case "hrsdc"_h:
		{
			typedef uint16_t IT; typedef std::complex<double> FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
			}
		}
		break;
		case "hrldc"_h:
		{
			typedef uint32_t IT; typedef std::complex<double> FT;
			IT mn[2]; FT h;
			for (uint32_t i=0; i!=nnz; ++i) {
				file.read((char*) &mn, 2*sizeof(IT));
				file.read((char*) &h, sizeof(FT));
				H.back().vec_.push_back({size_t(mn[0]),size_t(mn[1]),hel(h)});
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


// adapt ll_hbondss to structure
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
