// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include "c_switch.h"
#include "aux_parser.h"
#include "ll_io.h"
#include "ll_hbonds.h"
#include <fstream>
#include <string>

using namespace aux;
using namespace ll__;

void c_switch(const c_input& inp, std::ostream& os) {

	// read old 'bh.wad' file
	std::vector<size_t> Norb;
	fMat B, bnds;
	std::vector<ll__::i_i_R_H> dat;
	{
	
		// open file
		std::ifstream file;
		file.open(inp.bh, std::ios::binary);
		if (!file.good()) throw(std::invalid_argument("failed to open '"+inp.bh+"'"));
		
		// read header
		{
			char head[5];
			file.read((char*) &head, 5*sizeof(char));
			if (strncmp(head,"wad90",5)) throw(std::runtime_error(
						"bad header in file '"+inp.bh+"'"));
		}

		// read dim and #bonds
		uint32_t D, Nb;
		file.read((char*) &D, sizeof(uint32_t));
		file.read((char*) &Nb, sizeof(uint32_t));

		// read Norb
		{
			uint32_t N;
			file.read((char*) &N, sizeof(uint32_t));

			Norb.reserve(N);
			while (Norb.size()<Norb.capacity()) {
				file.read((char*) &N, sizeof(uint32_t));
				Norb.push_back(N);
			}
		}

		// read B and bonds
		B = fMat(D,D); bnds = fMat(D,Nb);
		file.read((char*) B.data(), D*D*sizeof(double));
		file.read((char*) bnds.data(), D*Nb*sizeof(double));

		// read pairs, R and H
		dat.reserve(Nb);
		while (dat.size()<dat.capacity()) {
			uint32_t tmp[3]; // i1,i2,size
			file.read((char*) tmp, 3*sizeof(uint32_t));

			fMat R(D,tmp[2]);
			file.read((char*) R.data(), R.size()*sizeof(double));

			std::vector<cMat> H; H.reserve(tmp[2]);
			while (H.size()<H.capacity()) {
				cMat cH(Norb[tmp[0]],Norb[tmp[1]]);
				file.read((char*) cH.data(), cH.size()*sizeof(std::complex<double>));
				H.push_back(cH);
			}

			dat.push_back(i_i_R_H(tmp[0],tmp[1],std::move(R),std::move(H)));
		}

		// close file
		file.close();
	}

	// read cell.psc file
	auto P = readPOSCAR(inp.psc);
	const auto J = aux::sorted_order(P.id.cbegin(),P.id.cend());
	aux::reorder(P.Ap.cBegin(),J); aux::reorder(P.id.begin(),J);	
	std::transform(P.id.cbegin(),P.id.cend(),P.id.begin(),
		[](const std::string& id) -> std::string {
			return id.substr(0,id.find_first_of("0123456789"));
		}
	);

	const ll_cell cell(std::move(B),std::move(P.Ap),
		std::vector<size_t>(P.Ap.N(),1),std::move(P.id));

	// brutally shove everything into poor ll_hbonds
	ll_hbonds(std::move(cell),std::move(Norb),std::move(dat)).
		writeToFile(inp.wbh);
}
