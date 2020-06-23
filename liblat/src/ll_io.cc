// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "ll_io.h"
#include "ll_fn.h"
#include "aux_io.h"
#include "aux_sort.h"
#include "aux_io.h"
#include <fstream>
#include <ctime>
#include <iomanip>
#include <cfloat>
#include <omp.h>
#include <limits>

using namespace ll__;


// omp functions
extern int omp_get_thread_num(void);


// dft related
ll__::psc ll__::readPOSCAR(const std::string& fileName, const size_t d) {
	const std::regex r("[\\s]+");

	auto file = aux::openFile<std::ifstream>(fileName);
	
	std::string line;
	std::getline(file,line);			// first line is a comment

	std::getline(file,line);			// 2nd line holds scaling constant
	double a0;
	sscanf(line.c_str(),"%lf",&a0);
	if (a0<0.0) a0=1.0;				// negative a0 means volume, no scaling

	auto B = a0*fMat(file,d,d).T();			// read next d lines into basis B

	std::getline(file,line);			// 6th line holds # atoms per species

	idv id;
	if (std::any_of(line.begin(),line.end(),	// check if ids for atoms are included
	[](const char c){return std::isalpha(c);})) {
		for (std::sregex_token_iterator
		i(line.begin(),line.end(),r,-1), e; i!=e; ++i)
			if (!i->str().empty()) id.push_back(*i);
		std::getline(file,line);				
	}

	std::vector<size_t> N; N.reserve(id.size());
	for (std::sregex_token_iterator
	i(line.begin(),line.end(),r,-1), e; i!=e; ++i)
		if (!i->str().empty()) N.push_back(std::stoul(*i));

	std::getline(file,line);
	if (line[0]=='S' || line[0]=='s')
		std::getline(file,line);
	
	bool cart = line[0]=='C' || line[0]=='c' ||	// check if atomic positions are in cartesian
			line[0]=='K' || line[0]=='k';

	size_t NN = std::accumulate(N.begin(),N.end(),0);
	auto Ap = fMat(file,NN,d).T();			// read next lines into atomic positions Ap
	if (cart) Ap = B.leftDivide(a0*Ap);		// if Ap cartesian transform into direct

	return {std::move(B),std::move(Ap),std::move(N),std::move(id)};
}
void ll__::printPOSCAR(const std::string& fileName, const fMat& B, const fMat& Ap,
		const aTv& T, const idv& id,
		const double a0, const bool direct,
		const size_t prec, const std::string& header) {
	
	using namespace aux;

	// bitching
	if (B.empty())
		throw(std::invalid_argument("B may not be empty"));
	if (!B.square())
		throw(std::invalid_argument("B has dimensions ("+std::to_string(B.M())+","+
				std::to_string(B.N())+"), need B square"));
	if (B.M()!=Ap.M())
		throw(std::invalid_argument("B.M()=="+std::to_string(B.M())+", Ap.M()=="+
				std::to_string(Ap.M())+", need B.M()==Ap.M()"));
	if (T.size() && T.size()!=Ap.N())
		throw(std::invalid_argument("size(T)=="+std::to_string(T.size())+", Ap.N()=="+
				std::to_string(Ap.N())+", need T empty or size(T)==Ap.N()"));
	if (header.find_first_of('\n')!=std::string::npos)
		throw(std::invalid_argument("new line character in header, invalid!"));

	// open file
	auto file = aux::openFile<std::ofstream>(fileName);

	// write header
	file << (header.empty() ? aux::timeStamp(): header) << "\n";

	// write lattice constant
	file << a0 << "\n";

	// write basis and id
	file << lm__::T(B).print(prec) << "\n";

	// write identifiers
	if (!id.empty()) file << id << "\n";

	// write # of atoms according to T, write Ap
	if (T.size()) {
		// sort Ap accoring to T
		auto Ap_ = lm__::T(Ap);
		aux::reorder(Ap_.rBegin(),aux::sorted_order(T.cbegin(),T.cend()));

		// find Ntype vector from type vector
		auto NT = T;
		std::sort(NT.begin(),NT.end());
		NT.erase(std::unique(NT.begin(),NT.end()),NT.end());
		std::for_each(NT.begin(),NT.end(),[&T](size_t& i){i=std::count(T.begin(),T.end(),i);});

		// print
		file << NT << "\n" << (direct ? "Direct\n": "Cartesian\n") << Ap_.print(prec) << "\n";
	} else {
		file << Ap.N() << "\n" << (direct ? "Direct\n": "Cartesian\n");
		file << lm__::T(Ap).print(prec) << "\n";
	}

	file.close();
}
void ll__::printPOSCAR(const std::string& fileName, const fMat& B, const fMat& Ap,
			const idv& id,
			const double a0, const bool direct,
			const size_t prec, const std::string& header) {
	assert(id.size()==Ap.N());
	assert(std::none_of(id.cbegin(),id.cend(),
		[](const auto& s)->bool{return s.empty();}));
	
	// get unique ids
	idv uid = id;
	std::sort(uid.begin(),uid.end());
	uid.resize(std::distance(uid.begin(),std::unique(uid.begin(),uid.end())));

	// get T vector from id and uid
	aTv T; T.reserve(id.size());
	for (const auto& s: id)
		T.push_back(std::distance(uid.begin(),std::lower_bound(uid.begin(),uid.end(),s)));

	printPOSCAR(fileName,B,Ap,T,uid,a0,direct,prec,header);
}


void ll__::printPOSCAR(const std::string& fileName, const ll_cell& inp,
			const bool direct, const bool strip,
			const size_t prec, const std::string& header) {

	// get type vector
	aTv T; T.reserve(inp.N());
	for (size_t i=0; i!=inp.N(); ++i)
		T.push_back(inp.type(i));

	// get id vector
	idv id; id.reserve(inp.id().size());
	for (const auto t: inp.types())
		id.push_back(strip ? inp.stripId(inp.id(t)): inp.id(t));

	// write poscar
	printPOSCAR(fileName,inp.B(),direct ? inp.Ap(): inp.getcAp(),
			T,id,1.0,direct,prec,header);
}
double ll__::readEf(const std::string& fileName) {
	std::string buffer;
	auto file = aux::openFile<std::ifstream>(fileName);
	
	double res;
	do std::getline(file,buffer); while(!file.eof() && buffer.find("E-fermi")==std::string::npos);
	if (file.eof()) throw(std::runtime_error(fileName+": unexpected end of file"));
	sscanf(buffer.c_str(),"%*s %*s %lf",&res);

	return res;
}


// omen related
rv ll__::genr(const size_t d) noexcept {
	assert(d>0 && d<=DIM__);

	rv res(DIM__,false);
	for (size_t i=1; i<=DIM__-d; ++i)
		res[i] = true;
	return res;
}
Ap_T ll__::readLayerMatrix(const std::string& fileName) {
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);

	// count number of lines, stream into sstr and close file
	std::stringstream sstr; size_t cnt=0;
	while (std::getline(file,line)) { ++cnt; sstr << line << '\n'; }
	file.close();

	// create output containers
	fMat Ap(DIM__,cnt);
	aTv T; T.reserve(cnt);

	// read data
	double fbuff;
	for (auto i=Ap.data(),e=Ap.data()+Ap.L(); i!=e; i+=DIM__) {
		std::getline(sstr,line);
		sscanf(line.c_str(),"%lf %lf %lf %lf",i,i+1,i+2,&fbuff);
		T.push_back(size_t(fbuff)-1);
	}

	// return result
	return {std::move(Ap*=10.0),std::move(T)};
}
void ll__::printOmf(const std::string& fileName, const vb_cb& E,
		    const std::vector<size_t>& Norb, const std::vector<double>& mass) {

	// bitching
	if (Norb.empty())
		throw(std::invalid_argument("Norb may not be empty"));
	if (!mass.empty() && mass.size()!=Norb.size())
		throw(std::invalid_argument("Norb.size() is "+std::to_string(Norb.size())+
			", mass.size() is "+std::to_string(mass.size())+
			", need Norb.size()==mass.size()"));

	// open file
	auto file = aux::openFile<std::ofstream>(fileName);
	file.precision(12); file.setf(std::ios::fixed);

	// write #anion #cation
	file << (Norb.size()/2 + Norb.size()%2) << " " << (Norb.size()/2) << "\n";

	// write bandgape, cb edge, vb edge
	file << (E.cb-E.vb) << " " << E.cb << " " << E.vb << "\n";

	// write #orbitals
	for (size_t i=0; i<Norb.size(); ++i)
		file << Norb[i] << ((i+1)%EPC__ ? " ": "\n");
	if (Norb.size()%EPC__) file << "\n";

	// write mass
	if (!mass.empty()) {
		for (size_t i=0; i!=mass.size(); ++i)
			file << std::setprecision(3) << mass[i] << ((i+1)%EPC__ ? " ": "\n");
		if (mass.size()%EPC__) file << "\n";
	}

	// close all
	file.close();
}
void ll__::printOlf(const std::string& fileName, const fMat& B, const fMat& Ap,
		    const idv& id, const size_t nn, const double bl) {

	// bitching
	if (B.empty())
		throw(std::invalid_argument("B may not be empty"));
	if (!B.square())
		throw(std::invalid_argument("B has dimensions ("+std::to_string(B.M())+","+
				std::to_string(B.N())+"), need B square"));
	if (B.M()!=Ap.M())
		throw(std::invalid_argument("B.M()=="+std::to_string(B.M())+", Ap.M()=="+
				std::to_string(Ap.M())+", need B.M()==Ap.M()"));
	if (id.size()!=Ap.N())
		throw(std::invalid_argument("size(id)=="+std::to_string(id.size())+", Ap.N()=="+
				std::to_string(Ap.N())+", need size(id)==Ap.N()"));

	// open file
	auto file = aux::openFile<std::ofstream>(fileName);

	// write header, i.e. #atoms #nn 0 0 0
	file << Ap.N() << " " << nn << " 0 0 0\n\n";

	// print basis
	file << lm__::T(B).print(16) << "\n\n";

	// find width of largest id
	const size_t idw = std::max_element(id.cbegin(),id.cend(),
		[](const auto& i, const auto& j)->bool
		{ return i.size()<j.size(); })->size()+3;

	// print atomic positions and id
	auto j = id.cbegin();
	for (auto i=Ap.ccBegin(),e=Ap.ccEnd(); i!=e; ++i,++j) {
		file << std::setw(idw) << std::left << *j;
		for (const auto f: *i)
			file << std::fixed << std::setw(24) <<
				std::left << std::setprecision(12) << f;
		file << "\n";
	}
	
	// print bond length
	file << "\n" << bl;

	// close file
	file.close();
}
omf ll__::readOmf(const std::string& fileName) {
	// stream and buffer and regex
	std::string line;
	std::regex r("[\\s,;]+");
	
	// open file
	auto file = aux::openFile<std::ifstream>(fileName);

	// create result
	omf res;

	// read number of species end resize result
	{
		size_t N1,N2;
		std::getline(file,line);
		sscanf(line.c_str(),"%lu %lu",&N1,&N2);
		res.Norb.reserve(N1+N2); res.mass.reserve(N1+N2);
	}

	// read vb_cb
	std::getline(file,line);
	sscanf(line.c_str(),"%*f %lf %lf",&res.E.cb,&res.E.vb);

	// read Norb
	std::getline(file,line);
	for (std::sregex_token_iterator i{line.begin(),line.end(),r,-1}, e; i!=e; ++i)
		res.Norb.push_back(std::stoul(*i));

	// read mass if present
	std::getline(file,line);
	if (line.empty() || file.eof())
		return res;

	for (std::sregex_token_iterator i{line.begin(),line.end(),r,-1}, e; i!=e; ++i)
		res.mass.push_back(std::stod(*i));
	res.mass.shrink_to_fit();

	return res;
}
olf ll__::readOlf(const std::string& fileName) {
	// stream and buffers
	std::string line;
	char buff[256];
	
	// open file
	auto file = aux::openFile<std::ifstream>(fileName);

	// read number of atoms, nn
	size_t N,nn;
	std::getline(file,line);
	sscanf(line.c_str(),"%lu %lu",&N,&nn);

	// create result
	olf res = {fMat(DIM__,DIM__),fMat(DIM__,N),{},{},nn,.0};
	res.id.reserve(N);

	// read basis
	std::getline(file,line);
	for (auto i=res.B.cBegin(),e=res.B.cEnd(); i!=e; ++i) {
		std::getline(file,line);
		sscanf(line.c_str(),"%lf %lf %lf",i->data(),i->data()+1,i->data()+2);
	}

	// read positions and id
	std::getline(file,line);
	for (auto i=res.Ap.cBegin(),e=res.Ap.cEnd(); i!=e; ++i) {
		std::getline(file,line);
		sscanf(line.c_str(),"%s %lf %lf %lf",buff,i->data(),i->data()+1,i->data()+2);
		res.id.push_back(buff);
	}

	// read bl
	std::getline(file,line);
	std::getline(file,line);
	sscanf(line.c_str(),"%lf",&res.bl);

	// close file
	file.close();

	// generate T
	idv uid; uid.reserve(res.id.size());
	for (const auto& s: res.id)
		if (std::find(uid.cbegin(),uid.cend(),s)==uid.cend())
			uid.push_back(s);
	res.T.reserve(res.id.size());
	for (const auto& s: res.id)
		res.T.push_back(std::distance(uid.cbegin(),std::find(uid.cbegin(),uid.cend(),s)));

	return res;
}


// wannier90 related
rv ll__::hrDim(const std::string& fileName) {
	// stream and read line buffer
	std::string line;

	auto file = aux::openFile<std::ifstream>(fileName);
	std::getline(file,line);

	// read header
	size_t NR, Nw; file >> Nw >> NR;
	if (!NR || !Nw) throw(std::invalid_argument(fileName+": failed to read header"));

	// skip to data
	for (size_t i=0; i<=std::ceil(NR/15.0); ++i)
		std::getline(file,line);

	// read upto center R vector
	int r0,r1,r2;
	rv res(DIM__,true);
	for (size_t i=0; i!=NR/2; ++i) {
		std::getline(file,line);
		sscanf(line.c_str(),"%i %i %i",&r0,&r1,&r2);
		res[0] = res[0] && !r0;
		res[1] = res[1] && !r1;
		res[2] = res[2] && !r2;
	}

	return res;
}
fMat ll__::readB(const std::string& fileName) {

	// stream and read line buffer
	std::string buffer;

	auto file = aux::openFile<std::ifstream>(fileName);
	
	// read basis vectors   
	fMat B(DIM__);
	do std::getline(file,buffer); while(!file.eof() && buffer.find("a_1")==std::string::npos);
	if (file.eof()) throw(std::runtime_error(fileName+": unexpected end of file"));

	for (size_t i=0; i<B.L(); i+=DIM__) {
		sscanf(buffer.c_str(),"%*s %lf %lf %lf",B.data()+i,B.data()+i+1,B.data()+i+2);
		std::getline(file,buffer);
		if (file.eof()) throw(std::runtime_error(fileName+": unexpected end of file"));
	}

	return B;
}
Ap_id ll__::readAp(const std::string& fileName, const bool direct, const std::vector<size_t>& abl) {

	// stream and read line buffer
	std::stringstream sstr;
	std::string line;

	// open file, skip to relevant part, count number of Ap and stream into sstr
	size_t Na=0;
	{
		auto file = aux::openFile<std::ifstream>(fileName);
		
		do std::getline(file,line);
		while(!file.eof() && line.find("Fractional Coordinate")==std::string::npos);
		std::getline(file,line);

		while (std::getline(file,line) && line.find('*')==std::string::npos) {
			++Na;
			sstr << line << '\n';
		}
		if (file.eof()) throw(std::runtime_error(fileName+": unexpected end of file"));
		file.close();
	}

	// check abl
	if (!std::is_sorted(abl.begin(),abl.end()))
		throw(std::invalid_argument("abl is not sorted"));
	if (std::adjacent_find(abl.begin(),abl.end())!=abl.end())
		throw(std::invalid_argument("abl has duplicate entries"));
	if (std::any_of(abl.begin(),abl.end(),[&Na](const size_t i){return i>=Na;}))
		throw(std::invalid_argument("abl holds illegal entries, need all <"+
					std::to_string(Na)));
	if (Na==abl.size()) return {{},{}};

	
	// read id, Ap
	fMat Ap(DIM__,0); Ap.reserve(Na-abl.size());
	std::vector<std::string> id; id.reserve(Na-abl.size());
	{
		char sbuff[16]; fMat mbuff(DIM__,1);
		const char* format = direct ? "%*s %*s %*i %lf %lf %lf":
					      "%*s %*s %*i %*f %*f %*f %*s %lf %lf %lf";
		auto b=abl.begin();
		for (size_t i=0; i!=Na; ++i) {
			std::getline(sstr,line);
		
			// check blacklist
			if (b!=abl.end() && i==*b) { ++b; continue; }
			
			// read type
			sscanf(line.c_str(),"%*s %s",sbuff);
			id.push_back(sbuff);
		
			// read Ap
			sscanf(line.c_str(),format,mbuff.data(),mbuff.data()+1,mbuff.data()+2);
			Ap.push_back(mbuff);
		}
	}

	// move into container and return
	return {std::move(Ap),std::move(id)};
}
Wp_s ll__::readWp(const std::string& fileName, const std::vector<size_t>& wbl) {

	// stream and read line buffer
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);

	// read number of wannier functions
	size_t Nw;
	{
		int ibuff;
		do std::getline(file,line);
		while(!file.eof() && line.find("Number of Wannier Functions")==std::string::npos);
		if (file.eof()) throw(std::runtime_error(fileName+": unexpected end of file"));
		sscanf(line.c_str(),"%*s %*s %*s %*s %*s %*c %i",&ibuff);
		Nw = size_t(ibuff);
	}
	
	// check wbl
	if (!std::is_sorted(wbl.begin(),wbl.end()))
		throw(std::invalid_argument("wbl is not sorted"));
	if (std::adjacent_find(wbl.begin(),wbl.end())!=wbl.end())
		throw(std::invalid_argument("wbl has duplicate entries"));
	if (std::any_of(wbl.begin(),wbl.end(),[&Nw](const size_t i){return i>=Nw;}))
		throw(std::invalid_argument("abl holds illegal entries, need all <"+
					std::to_string(Nw)));
	if (Nw==wbl.size()) return {{},{}};

	// skip to relevant part
	do std::getline(file,line);
	while(!file.eof() && line.find("Final State")==std::string::npos);
	if (file.eof()) throw(std::runtime_error(fileName+": unexpected end of file"));

	// read wannier centers and spread
	fMat Wp(DIM__,0); Wp.reserve(Nw-wbl.size());
	fMat S(1,0); S.reserve(Nw-wbl.size());
	{
		fMat mbuff(DIM__,1), sbuff(1,1);
		auto b=wbl.begin();
		for (size_t i=0; i<Nw; ++i) {
			std::getline(file,line);
			if (file.eof()) throw(std::runtime_error(fileName+": unexpected end of file"));
	
			// check blacklist
			if (b!=wbl.end() && i==*b) { ++b; continue; }

			// read Wp and spread
			sscanf(line.c_str(),"%*s %*s %*s %*s %*i %*c %lf %*c %lf %*c %lf %*c %lf",
				mbuff.data(),mbuff.data()+1,mbuff.data()+2,sbuff.data());
			Wp.push_back(mbuff);
			S.push_back(sbuff);
		}
	}
	file.close();

	return {std::move(Wp),std::move(S)};
}
R_H<> ll__::readOp(std::istream& strm, const char* frmt,
		const size_t Nw, const size_t NR,
		const std::vector<size_t>& wbl,
		const std::function<bool(const cMat&)>& keep) {
	
	// check wbl
	if (!std::is_sorted(wbl.begin(),wbl.end()))
		throw(std::invalid_argument("wbl is not sorted"));
	if (std::adjacent_find(wbl.begin(),wbl.end())!=wbl.end())
		throw(std::invalid_argument("wbl has duplicate entries"));
	if (std::any_of(wbl.begin(),wbl.end(),[&Nw](const size_t i){return i>=Nw;}))
		throw(std::invalid_argument("abl holds illegal entries, need all <"+
					std::to_string(Nw)));
	if (Nw==wbl.size()) return {{},{}};
	
	// create result containers
	fMat R(DIM__,0); R.reserve(NR);
	std::vector<cMat> H; H.reserve(NR);
	
	// read data
	std::string line;
	fMat rbuff(DIM__,1); cMat mat(Nw); double re,im;
	std::getline(strm,line);
	for (size_t r=0; r!=NR; ++r) {
		
		// read R vector
		sscanf(line.c_str(),"%lf %lf %lf",rbuff.data(),rbuff.data()+1,rbuff.data()+2);

		// operator data
		for (auto& i: mat) {
			// check premature eof
			if (strm.eof()) throw(std::runtime_error("unexpected eof"));
		
			sscanf(line.c_str(),frmt,&re,&im);
			i = {re,im};	
			std::getline(strm,line);
		}

		// apply wbl and check for relevancy
		mat.cRm(wbl); mat.rRm(wbl);
		if (keep(mat)) {
			R.push_back(rbuff);
			H.push_back(mat);
		}
		mat.resize(Nw,Nw);
	}

	return {std::move(R),std::move(H)};
}
R_H<> ll__::readHr(const std::string& fileName,
		const std::vector<size_t>& wbl,
		const std::function<bool(const cMat&)>& keep) {

	// stream and read line buffer
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);
	std::getline(file,line);

	// read Nw
	unsigned int Nw;
	std::getline(file,line);
	sscanf(line.c_str(),"%u",&Nw);
	
	// read NR
	unsigned int NR;
	std::getline(file,line);
	sscanf(line.c_str(),"%u",&NR);

	// skip to data
	for (size_t i=0; i<std::ceil(NR/15.0); ++i)
		std::getline(file,line);

	return readOp(file,"%*i %*i %*i %*i %*i %lf %lf",Nw,NR,wbl,keep);
}
template<class MT>
void ll__::writeHr(const std::string& fileName, const R_H<MT>& hr, const std::string& header) {
	assert(!hr.empty());
	assert(hr.front().square());
	assert(std::all_of(hr.cbegin(),hr.cend(),[&hr](const auto& i)->bool
		{ return size(i)==size(hr.front()); }));
	if (header.find_first_of('\n')!=std::string::npos)
		throw(std::invalid_argument("new line character in header, invalid!"));
	
	// open file
	auto file = aux::openFile<std::ofstream>(fileName);

	// write header
	if (!header.size()) {
		std::locale::global(std::locale("en_US.utf8"));
		std::time_t t = std::time(nullptr);
		char mbstr[100];
		std::strftime(mbstr, sizeof(mbstr), "%c", std::localtime(&t));
		
		file << mbstr << "\n";
	} else
		file << header << "\n";

	// write Nw, Nr
	file << std::setw(12) << hr.Nw() << "\n";
	file << std::setw(12) << hr.N() << "\n";

	// fake degeneracy
	for (size_t i=1; i<=hr.N(); ++i)
		file << std::setw(5) << 1 << (!(i%15) && i!=hr.N() ? "\n": "");

	// write data
	auto i = hr.cbegin();
	for (auto r=hr.ccBegin(),re=hr.ccEnd(); r!=re; ++r,++i) {
		auto j=i->cbegin();
		for (size_t m=1; m<=hr.Nw(); ++m)
		for (size_t n=1; n<=hr.Nw(); ++n,++j) {
			// r vector
			file << "\n";
			for (auto ri=r->cbegin(),rie=r->cend(); ri!=rie; ++ri)
				file << std::setw(5) << int(*ri);

			// wannier indices
			file << std::setw(5) << m;
			file << std::setw(5) << n;

			// wannier data;
			file << std::setw(12) << std::fixed << std::setprecision(6) << std::real(*j);
			file << std::setw(12) << std::fixed << std::setprecision(6) << std::imag(*j);
		}
	}

	// close file
	file.close();
}
template void ll__::writeHr<fMat>(const std::string& fileName, const R_H<fMat>& hr, const std::string& header);
template void ll__::writeHr<cMat>(const std::string& fileName, const R_H<cMat>& hr, const std::string& header);

std::vector<ll__::cMat> ll__::readAMN(const std::string& fileName) {
	
	// stream and read line buffer
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);
	std::getline(file,line);
	
	// read Nb,Nk,Norb
	unsigned int Nb,Nk,Norb;
	std::getline(file,line);
	sscanf(line.c_str(),"%u %u %u",&Nb,&Nk,&Norb);

	// read data
	double im,re;
	std::vector<cMat> res; res.reserve(Nk);
	while (res.size()!=res.capacity()) {
		res.push_back(cMat(Nb,Norb));
		for (auto& i: res.back()) {
			std::getline(file,line);
			sscanf(line.c_str(),"%*i %*i %*i %lf %lf",&re,&im);
			i = {re,im};
		}
	}

	return res;
}

R_H<> ll__::readXYZ(const std::string& fileName, const size_t d,
		const std::vector<size_t>& wbl,
		const std::function<bool(const cMat&)>& keep) {

	// stream and read line buffer
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);
	std::getline(file,line);

	// read Nw
	unsigned int Nw;
	std::getline(file,line);
	sscanf(line.c_str(),"%u",&Nw);
	
	// read NR
	unsigned int NR;
	std::getline(file,line);
	sscanf(line.c_str(),"%u",&NR);

	switch (d) {
		case 0: return readOp(file,"%*i %*i %*i %*i %*i %lf %lf",
				Nw,NR,wbl,keep); break;
		case 1: return readOp(file,"%*i %*i %*i %*i %*i %*lf %*lf %lf %lf",
				Nw,NR,wbl,keep); break;
		case 2: return readOp(file,"%*i %*i %*i %*i %*i %*lf %*lf %*lf %*lf %lf %lf",
				Nw,NR,wbl,keep); break;
		default: return {{},{}}; break;
	}
}
fMat ll__::readEig(const std::string& fileName) {
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);

	// find Nw, Nk
	size_t Nw, Nk;
	{
		std::string lline;
		while (std::getline(file,line))
			lline=std::move(line);
		sscanf(lline.c_str(),"%lu %lu",&Nw,&Nk);
	}
	file.close();


	// reopen file, read data and return result
	file.open(fileName);
	if (!file.is_open()) throw std::invalid_argument("failed to open " + fileName);	

	fMat res(Nw,Nk);
	for (auto i=res.data(),e=res.data()+res.L(); i!=e; ++i) {
		std::getline(file,line);
		sscanf(line.c_str(),"%*i %*i %lf",i);
	}	
	return res;
}
k_U ll__::readWannierTransf(const std::string& fileName) {
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);
	
	// skip header
	std::getline(file,line);

	// read Nk, Nw
	size_t Nk,Nw,Nw2;
	std::getline(file,line);
	sscanf(line.c_str(),"%lu %lu %lu",&Nk,&Nw,&Nw2);
	if (Nw!=Nw2) throw(std::runtime_error("matrices not square in '"+fileName+"'"));

	// result containers
	fMat k(DIM__,Nk);
	std::vector<cMat> U(Nk,cMat(Nw,Nw));

	// read k points and matrices
	auto i = k.cBegin();
	for (auto& u: U) {
		// read empty line
		std::getline(file,line);

		// read k points line
		std::getline(file,line);
		sscanf(line.c_str(),"%lf %lf %lf",i->data(),i->data()+1,i->data()+2);
		++i;

		// read matrix elements
		double re,im;
		for (auto& j: u) {
			std::getline(file,line);
			sscanf(line.c_str(),"%lf %lf",&re,&im);
			j = {re,im};
		}
	}
	file.close();
	
	return {std::move(k),std::move(U)};
}
k_U ll__::readWannierTransf(const std::string& fileNameU, const std::string& fileNameUdis) {
	std::string line;

	// open files
	auto file1 = aux::openFile<std::ifstream>(fileNameU);
	auto file2 = aux::openFile<std::ifstream>(fileNameUdis);

	// skip headers
	std::getline(file1,line);
	std::getline(file2,line);

	// read Nk,Nw,Nb
	size_t Nk,Nw,Nb;
	{
		size_t ck1,ck2;
		
		std::getline(file1,line);
		sscanf(line.c_str(),"%lu %lu %lu",&Nk,&Nw,&ck1);
		if (Nw!=ck1) throw(std::runtime_error("matrices not square in '"+fileNameU+"'"));

		std::getline(file2,line);
		sscanf(line.c_str(),"%lu %lu %lu",&ck1,&ck2,&Nb);
		if (Nk!=ck1) throw(std::runtime_error("Nk in file '"+fileNameU+"' and Nk in file '"+
					fileNameUdis+"' do not match"));
		if (Nw!=ck2) throw(std::runtime_error("Nw in file '"+fileNameU+"' and Nw in file '"+
					fileNameUdis+"' do not match"));
	}

	// result containers
	fMat k(DIM__,Nk);
	std::vector<cMat> U; U.reserve(Nk);

	// read k points and matrices
	fMat ck(DIM__,1);
	for (auto i=k.cBegin(),e=k.cEnd(); i!=e; ++i) {
		// read empty line
		std::getline(file1,line);
		std::getline(file2,line);

		// read k points line
		std::getline(file1,line);
		sscanf(line.c_str(),"%lf %lf %lf",i->data(),i->data()+1,i->data()+2);
		std::getline(file2,line);
		sscanf(line.c_str(),"%lf %lf %lf",ck.data(),ck.data()+1,ck.data()+2);
		if (*i!=ck) throw(std::invalid_argument("k points mismatch in file '"+
				fileNameU+"' and file '"+fileNameUdis+"'"));

		// read matrix elements
		double re,im;
		cMat cU(Nw,Nw);
		for (auto& j: cU) {
			std::getline(file1,line);
			sscanf(line.c_str(),"%lf %lf",&re,&im);
			j = {re,im};
		}
		cMat Udis(Nb,Nw);
		for (auto& j: Udis) {
			std::getline(file2,line);
			sscanf(line.c_str(),"%lf %lf",&re,&im);
			j = {re,im};
		}

		// extend result
		U.push_back(Udis.prod(cU));
	}

	return {std::move(k),std::move(U)};
}
k_U ll__::readChk(const std::string& fileName) {
	std::string line;
	auto file = aux::openFile<std::ifstream>(fileName);
	std::getline(file,line);
	
	// read num_bands
	size_t Nb;
	std::getline(file,line);
	sscanf(line.c_str(),"%lu",&Nb);

	// read num_exclude_bands
	size_t Neb;
	std::getline(file,line);
	sscanf(line.c_str(),"%lu",&Neb);
	for (size_t i=0; i!=Neb; ++i)
		std::getline(file,line);
	std::getline(file,line);
	std::getline(file,line);

	// read num_kpt
	size_t Nk;
	std::getline(file,line);
	sscanf(line.c_str(),"%lu",&Nk);
	std::getline(file,line);

	// result containers
	fMat k(DIM__,Nk);
	std::vector<cMat> U; U.reserve(Nk);

	// read kpoints
	for (auto i=k.cBegin(),e=k.cEnd(); i!=e; ++i) {
		std::getline(file,line);
		sscanf(line.c_str(),"%lf %lf %lf",i->data(),i->data()+1,i->data()+2);
	}
	std::getline(file,line);

	// read num_wann
	size_t Nw;
	std::getline(file,line);
	sscanf(line.c_str(),"%lu",&Nw);
	std::getline(file,line);

	// read have_disentangled
	int dis;
	std::getline(file,line);
	sscanf(line.c_str(),"%i",&dis);

	// read matrices
	double re,im;
	if (dis) {
		// skip stuff
		for (size_t i=0, e=(Nb+1)*Nk+1; i!=e; ++i)
			std::getline(file,line);

		for (size_t ik=0; ik!=Nk; ++ik) {
			U.push_back(cMat(Nb,Nw));
			for (auto& i: U.back()) {
				std::getline(file,line);
				sscanf(line.c_str(),"%lf %lf",&re,&im);
				i = {re,im};
			}
		}
		cMat tmp(Nw,Nw);
		for (size_t ik=0; ik!=Nk; ++ik) {
			for (auto& i: tmp) {
				std::getline(file,line);
				sscanf(line.c_str(),"%lf %lf",&re,&im);
				i = {re,im};
			}
			U[ik] = U[ik].prod(tmp);
		}
	} else {
		cMat tmp(Nw,Nw);
		for (size_t ik=0; ik!=Nk; ++ik) {
			U.push_back(cMat(Nb,Nw));
			for (auto& i: U.back()) {
				std::getline(file,line);
				sscanf(line.c_str(),"%lf %lf",&re,&im);
				i = {re,im};
			}
		}
	}

	// return result
	return {std::move(k),std::move(U)};
}
void ll__::printWiToFile(const std::string& fileName, const wi& I) {
	using namespace aux;

	// open file
	auto file = aux::openFile<std::ofstream>(fileName);

	// write wi
	file << I;	
}

// gnu plot related
void ll__::writeGnuPlot(const fMat& xdat, const fMat& ydat,
		const std::string& seed, const std::string& prefix, std::vector<double> window,
		const std::vector<double>& xticks, const idv& xlabels,
		const std::vector<double>& yticks, const idv& ylabels) {
	
	assert(xdat.size()==ydat.size());
	assert((window.size()==4 && window[0]<=window[1] && window[2]<=window[3]) || !window.size());
	assert(xticks.size()==xlabels.size() || !xlabels.size());
	assert(yticks.size()==ylabels.size() || !ylabels.size());

	// find window if not specified
	if (!window.size()) {
		window.reserve(4);
		window.push_back(min(xdat));
		window.push_back(max(xdat));
		window.push_back(min(ydat));
		window.push_back(max(ydat));
	}

	// open gnu file
	auto file = aux::openFile<std::ofstream>(prefix+seed+".gnu");

	// write gnu file
	file << "set style data dots\n";
	file << "set nokey\n";
	file << "set xrange [" << std::setprecision(6) << window[0] << " : "
			       << std::setprecision(6) << window[1] << "]\n";
	file << "set yrange [" << std::setprecision(6) << window[2] << " : "
			       << std::setprecision(6) << window[3] << "]\n";
	if (xticks.size()) {
		file << "set xtics (";
		if (xlabels.size()) {
			for (size_t i=0; i!=xticks.size(); ++i)
				file << "\"" << xlabels[i] << "\" " << std::setprecision(3) << xticks[i]
				<< (i!=xticks.size()-1 ? ",": ")\n");
		} else {
			for (size_t i=0; i!=xticks.size(); ++i)
				file << std::setprecision(3) << xticks[i]
				     << (i!=xticks.size()-1 ? ",": ")\n");
		}
	}
	if (yticks.size()) {
		file << "set ytics (";
		if (ylabels.size()) {
			for (size_t i=0; i!=yticks.size(); ++i)
				file << "\"" << ylabels[i] << "\" " << std::setprecision(3) << yticks[i]
				<< (i!=yticks.size()-1 ? ",": ")\n");
		} else {
			for (size_t i=0; i!=yticks.size(); ++i)
				file << std::setprecision(3) << yticks[i]
				     << (i!=yticks.size()-1 ? ",": ")\n");
		}
	}
	file << "plot \"" << seed << ".dat" << "\"\n";

	// close gnu file, open dat file
	file.close();
	file.open(prefix+seed+".dat");
	if (!file.good()) throw(std::invalid_argument("failed to open '"+prefix+seed+".dat'"));

	// dump xdat,ydat into dat file
	for (auto i=xdat.cbegin(),j=ydat.cbegin(),e=xdat.cend(); i!=e; ++i,++j)
		file << std::setw(24) << std::setprecision(16) << std::fixed << *i
		     << std::setw(24) << std::setprecision(16) << std::fixed << *j << "\n";

	// close dat file
	file.close();
}


// print structure report
template<class WT>
std::ostream& ll__::printStructureReport(const struct_report& rep, const size_t mask,
			const WT& W, std::ostream& os) noexcept {

	if (mask & size_t(1)) {
		os << "exact matches ["
		   << std::accumulate(rep.exact.cbegin(),rep.exact.cend(),size_t(0),
			[](const size_t s, const auto& i) -> size_t { return s+i.N; })
		   << "]";
		for (const auto& m: rep.exact)
			os << "\n "
			   << std::setw(16) << std::setfill(' ') << std::right << W.id(m.i1)
			   << "   ->   "
			   << std::setw(16) << std::setfill(' ') << std::right << W.id(m.i2)
			   << std::setw(12) << std::setfill(' ') << std::right << m.N << " times";
		os << "\n\n";
	}
	if (mask & size_t(2)) {
		os << "approximate matches ["
		   << std::accumulate(rep.approx.cbegin(),rep.approx.cend(),size_t(0),
			[](const size_t s, const auto& i) -> size_t { return s+i.N; })
		   << "]";
		for (const auto& m: rep.approx)
			os << "\n "
			   << std::setw(16) << std::setfill(' ') << std::right << W.id(m.i1)
			   << "   ->   "
			   << std::setw(16) << std::setfill(' ') << std::right << W.id(m.i2)
			   << std::setw(18) << std::setfill(' ') << std::left  << (" ("+W.id(m.j)+")")
			   << std::setw(8)  << std::setfill(' ') << std::right << m.N << " times";
		os << "\n\n";
	}
	if (mask & size_t(4))
		os << "unmatched [" << rep.unmatched.size() << "]\n"
		   << " " << rep.unmatched << "\n";

	return os;
}
template std::ostream& ll__::printStructureReport(
	const struct_report& rep, const size_t mask,
	const ll_hbonds& W, std::ostream& os) noexcept;
template std::ostream& ll__::printStructureReport(
	const struct_report& rep, const size_t mask,
	const ll_hbondss& W, std::ostream& os) noexcept;


// extract from wbh files
ll_cell ll__::extractCell(const std::string& fileName) {
	
	// open file
	auto file = aux::openFile<std::ifstream>(fileName,std::ios::binary);
	
	// read header
	{
		char head[5];
		file.read((char*) &head, 5*sizeof(char));
		if (strncmp(head,"wad90",5)) throw(std::runtime_error("bad header in file '"+fileName+"'"));
	}

	// read cell
	uint32_t D, N;
	file.read((char*) &D, sizeof(uint32_t));
	file.read((char*) &N, sizeof(uint32_t));
	
	fMat B(D,D), Ap(D,N);
	file.read((char*)  B.data(), D*D*sizeof(double));
	file.read((char*) Ap.data(), D*N*sizeof(double));
		
	idv id; id.reserve(N);
	uint8_t l; char buff[(int)std::numeric_limits<uint8_t>::max()+1];
	while (id.size()<id.capacity()) {
		file.read((char*) &l, sizeof(uint8_t));
		file.read((char*) buff, l*sizeof(char));
		id.push_back(buff);
	}

	set_mtol(WTOL__);
	ll_cell res(std::move(B),std::move(Ap),std::move(id));
	reset_mtol();

	return res;
}

