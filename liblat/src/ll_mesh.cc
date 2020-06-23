// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_mesh.h"
#include "lm_fn.h"
#include "ll_lambda.h"

using namespace ll__;


// printing
template<class MT>
void ll_mesh<MT>::writeToFile(const std::string& fileName, const bool noheader) const {
	// open file
	std::ofstream file;
	file.open(fileName,std::ios::binary);
	if (!file.good())
		throw(std::invalid_argument("open file \'"+fileName+"\' failed"));

	// write header
	if (!noheader) {
		const std::string hdr = base_.cpx() ? "cMat": "fMat";
		file.write(hdr.c_str(),5);
	}

	// 64 bit buffer
	uint64_t dat;

	// write dimensions
	file.write((char*) &(dat = meshDim()+1), sizeof(uint64_t));
	file.write((char*) &(dat = M()), sizeof(uint64_t));
	
	for (const auto m: maj())
		file.write((char*) &(dat=D(m)), sizeof(uint64_t));

	// write matrix data
	file.write((char*) base_.data(), base_.L()*sizeof(double));
	file.close();
}
template void ll_mesh<fMat>::writeToFile(const std::string& fileName, const bool noheader) const;
template void ll_mesh<cMat>::writeToFile(const std::string& fileName, const bool noheader) const;


// generic mesh generators
ll_mesh<> ll__::genMesh(const fMat& bounds, std::vector<size_t> maj, std::vector<size_t> D) noexcept {
	if (!D.size()) return ll_mesh<>();
	assert(msize(bounds)==lm_size({D.size(),2}));
	assert(maj.size()==D.size());
	assert(std::all_of(D.cbegin(),D.cend(),[](const size_t i)->bool{return i;}));

	fMat base(D.size(),0); base.reserve(std::accumulate(D.cbegin(),D.cend(),
				size_t(1),std::multiplies<size_t>()));

	fMat steps = bounds.cBack()-bounds.cFront();
	auto j = D.cbegin();
	for (auto i=steps.begin(),e=steps.end(); i!=e; ++i,++j)
		*i = (*j==1 ? 0.0: *i/double(*j-1));

	auto n = zeros<fMat>(D.size(),1);
	std::function<void(const std::vector<size_t>::const_iterator i,
			   const std::vector<size_t>::const_iterator& j)> genBase;
	genBase = [&genBase,&base,&bounds,&steps,&D,&n]
				(const std::vector<size_t>::const_iterator i,
				 const std::vector<size_t>::const_iterator& j)->void {
		assert(j<=i);
		for (n[*i]=.0; n[*i]!=D[*i]; ++n[*i])
			if (i==j) base.push_back(bounds.cFront()+n*steps);
			else	  genBase(i-1,j);
	};

	genBase(maj.cend()-1,maj.cbegin());

	return ll_mesh<>(std::move(base),std::move(maj),std::move(D));
}
ll_mesh<> ll__::genMesh(const double lb, const double ub, const size_t D,
		std::vector<size_t> maj, const rv& r) noexcept {
	if (!D || maj.empty()) return ll_mesh<>();
	assert(maj.size()==r.size());

	fMat bounds(maj.size(),2); bounds.cFront()=lb; bounds.cBack()=ub;
	std::vector<size_t> D_(maj.size(),D);

	for (const auto i: inds(r))
		bounds(i,0)=0.0,
		bounds(i,1)=0.0,
		D_[i]=1;

	return genMesh(bounds,std::move(maj),std::move(D_));
}

// mesh in a cell with given density generator
ll_mesh<> ll__::genMesh_cell(const double rho, const fMat& B, const double lb, const double ub,
		std::vector<size_t> maj, const rv& r) {
	assert(r.size()==maj.size());
	assert(B.square() && B.M()==r.size());
	assert(rank(B)==B.M());
	assert(lb<ub);
	assert(rho>0.0);

	auto I = ninds(r);
	const auto RB = B.get({},I)*(ub-lb);
	
	const double N = rho*std::abs(det(ncat(RB,complement(RB))));
	fMat n = mnorm(RB)*std::pow(N /
			std::accumulate(RB.ccBegin(),RB.ccEnd(),double(1.0),
			[](const double s, const auto& i)->double{return s*norm(i);}),
				1/double(I.size()));

	std::vector<size_t> D(r.size(),1);
	while (!I.empty()) {
		const auto rn = round(n);
		const auto rnd = abs(n-rn);
		const auto pos = std::distance(rnd.cbegin(),
				 std::min_element(rnd.cbegin(),rnd.cend()));

		// set result
		D[I[pos]] = rn[pos];

		// rescale and remove used entries
		n *= std::pow(n[pos]/rn[pos],1.0/double(I.size()-1));
		n.cRm(pos), I.erase(I.cbegin()+pos);
	}

	// get bounds
	fMat bounds(r.size(),2); bounds.cFront()=lb; bounds.cBack()=ub;
	for (const auto i: inds(r))
		bounds(i,0)=0.0,
		bounds(i,1)=0.0;

	// generate standard mesh
	return genMesh(bounds,std::move(maj),std::move(D));
}
	
// mesh on an integer grid
ll_mesh<> ll__::genMesh_int(const fMat& bounds, std::vector<size_t> maj) noexcept {
	assert(bounds==lm__::round(bounds));
	assert(bounds.M()==maj.size());

	std::vector<size_t> D(maj.size());
	for (size_t i=0; i!=D.size(); ++i)
		D[i] = (size_t)std::abs(bounds(i,1)-bounds(i,0))+1;

	return genMesh(bounds,std::move(maj),std::move(D));
}
