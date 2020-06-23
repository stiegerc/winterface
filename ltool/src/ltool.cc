// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include "b_switch.h"
#include "f_switch.h"
#include "l_switch.h"
#include "p_switch.h"
#include "w_switch.h"
#include "c_switch.h"
#include "ll_omen.h"
#include "h_switch.h"
#include "ll_io.h"
#include "ll_fn.h"
#include "aux_io.h"


using namespace aux;
using namespace lm__;
using namespace ll__;

#define OPTS__ ":b::l::f::i::h::p::q::w::z::c::"
#define ERRMSG__ "valid options are:\n" \
		" -i: to generate OMEN input\n" \
		" -h: to write matrices for OMEN\n" \
		" -b: to calculate bandstructures\n" \
		" -l: to manipulate lattices\n" \
		" -f: to convert files to poscars\n" \
		" -p: to find primitive unit cells\n" \
		" -q: to query the contents of a wannier bonds hamiltonian\n" \
		" -w: to generate a wannier bonds hamiltonian from wannier90 output\n" \
		" -c: to convert deprecated 'bh.wad' into 'wannier90.wbh'\n\n" \
		"the following options will overwrite their input script counterparts:\n" \
		" --direct: output in direct coordinates\n" \
		" --cartesian: output in cartesian coordinates\n" \
		" --pscin: input POSCAR file\n" \
		" --pscout: output POSCAR file\n" \
		" --r: restriction vector, used to specify vacuum dimensions\n" \
		" --vac: new vacuum level\n" \
		" --bond_factor: factor used when generating bonds"
		


// complain and quit macros
#define QUIT_1__(msg) {os << RED__ "ERROR: " << RESET__ << msg << "\n"; return 0;}
#define QUIT_2__(inp,msg) {os << inp << "\n\n"; QUIT_1__(msg);}
#define QUIT_SELECT__(_1,_2,NAME,...) NAME
#define QUIT__(...) QUIT_SELECT__(__VA_ARGS__, \
		QUIT_2__, QUIT_1__)(__VA_ARGS__)


int main(const int argc, char** argv) {

#ifndef NLOGFILE_
//	auto fs = aux::openFile<std::ofstream>(aux::timeStamp()+".log");
	auto fs = aux::openFile<std::ofstream>("ltool.log",std::ios::app);
	fs << "\n\n" << aux::timeStamp() << "\n\n";
	aux_tee os(std::cout,fs);
#else
	std::ostream& os = std::cout;
#endif

	// print errmsg if there are no arguments
	if (argc==1) os << ERRMSG__ << "\n";

	// commands and arguments
	struct co_ {
		int cmd;
		std::string opt;
	};
	std::vector<co_> cmds;


	// flags and direct options
	int direct_flag = -1;
	int strip_flag = -1;
	std::string prefix="";
	std::string pscin="", pscout="";
	std::string wout="", hrdat="";
	std::string wbh="";
	size_t verbosity=NPOS__;
	rv r; r.reserve(DIM__);
	double vac=-1.0;
	double bond_factor=0.0;
	int all_range_flag = 0;


	// parse commands
	int cmd;
	while (true) {
		
		// read options
		int options_index = -1;
		static struct option long_options[] = {
			{"direct", no_argument, &direct_flag, 1},
			{"cartesian", no_argument, &direct_flag, 0},
			{"strip", no_argument, &strip_flag, 1},
			{"nostrip", no_argument, &strip_flag, 0},
			{"prefix", required_argument, nullptr, 0},
			{"pscin", required_argument, nullptr, 0},
			{"pscout", required_argument, nullptr, 0},
			{"wout", required_argument, nullptr, 0},
			{"hrdat", required_argument, nullptr, 0},
			{"wbh", required_argument, nullptr, 0},
			{"verbosity", required_argument, nullptr, 0},
			{"r", required_argument, nullptr, 0},
			{"vac", required_argument, nullptr, 0},
			{"bond_factor", required_argument, nullptr, 0},
			{"all_range", no_argument, &all_range_flag, 1}
		};
		cmd = getopt_long(argc,argv,OPTS__,
				long_options,&options_index);

		// parse long options
		if (options_index!=-1) {
			switch (fnvHash(long_options[options_index].name)) {
				case "prefix"_h: prefix = optarg; break;
				case "pscin"_h: pscin = optarg; break;
				case "pscout"_h: pscout = optarg; break;
				case "wout"_h: wout = optarg; break;
				case "hrdat"_h: hrdat = optarg; break;
				case "wbh"_h: wbh = optarg; break;
				case "verbosity"_h:
					try { verbosity = std::stoul(optarg); }
					catch (...) { QUIT__("--verbosity: bad argument"); }
				break;
				case "r"_h:
					for (auto c=optarg; *c; ++c) {
						switch (*c) {
						case '0': r.push_back(false); break;
						case '1': r.push_back(true); break;
						default: QUIT__("--r: bad argument");
						}
					}
				break;
				case "vac"_h:
					try { vac = std::stod(optarg); }
					catch (...) { QUIT__("--vac: bad argument"); }
				break;
				case "bond_factor"_h:
					try { bond_factor = std::stod(optarg); }
					catch (...) { QUIT__("--bond_factor: bad argument"); }
				break;
			}

			continue;
		}
	
		// store commands
		if (cmd==-1) break;
		if (cmd=='?') QUIT__(std::string("-")+std::string(1,optopt)+
					std::string(": unknown option\n\n")+std::string(ERRMSG__));
		cmds.push_back(co_{cmd,optarg!=nullptr?optarg:""});
	}


	// run commands
	for (const auto& c: cmds)
	switch (c.cmd) {
	
	case 'b':
	try {
		b_input inp = c.opt.empty() ? b_input(): parseFile<b_input>(c.opt,os);
		if (!prefix.empty()) inp.prefix = prefix;
		if (!wout.empty()) inp.wout = wout;
		if (!hrdat.empty()) inp.hrdat = hrdat;
		if (!wbh.empty()) inp.wbh = {wbh};
		if (verbosity!=NPOS__) inp.verbosity = verbosity;
		b_switch(inp,os);
	} catch(const std::exception& e) {
		QUIT__(b_input(),e.what());
	}
	break;
	case 'l':
	try {
		// read input script and overwrite entries if needed
		l_input inp = c.opt.empty() ? l_input(): parseFile<l_input>(c.opt,os);
		if (direct_flag!=-1) inp.direct = direct_flag;
		if (strip_flag!=-1) inp.strip = strip_flag;
		if (!pscin.empty()) inp.pscin = pscin;
		if (!pscout.empty()) inp.pscout = pscout;
		if (!wbh.empty()) inp.wbh = wbh;
		if (!r.empty()) inp.r = r;
		if (vac!=-1.0) inp.vac = {vac};
		if (bond_factor!=0.0) inp.bond_factor = bond_factor;
		if (!prefix.empty()) inp.prefix = prefix;
		if (verbosity!=NPOS__) inp.verbosity = verbosity;

		l_switch(inp,os);
	} catch(const std::exception& e) {
		QUIT__(l_input(),e.what());
	}
	break;
	case 'f':
	try {
		// detect file lambda
		const auto df = [](const std::string& fileName) -> const char* {
			
			// open file
			std::ifstream file;
			file.open(fileName);
			if (!file.good())
				throw(std::invalid_argument("open file \'"+
						std::string(fileName)+"\' failed"));

			std::string line;
			std::getline(file,line);
			
			// if first line starts with number assume lattice_dat type
			if (!line.empty() && std::isdigit(line[0]))
				return "omen";

			std::getline(file,line);

			// if second line holds wannier90 style output assume wout type
			if (!line.empty() && line.find("+---")!=std::string::npos)
				return "wannier";

			// else assume input script
			return "script";
		};

		f_input inp;
		switch (fnvHash(df(c.opt))) {
		case "wannier"_h: inp.wout = c.opt; break;
		case "omen"_h: inp.lattice_dat = c.opt; break;
		case "script"_h: inp = parseFile<f_input>(c.opt,os); break;
		}	
	
		if (direct_flag!=-1) inp.direct = direct_flag;
		if (strip_flag!=-1) inp.strip = strip_flag;
		if (!pscout.empty()) inp.pscout = pscout;
		if (!prefix.empty()) inp.prefix = prefix;
		if (verbosity!=NPOS__) inp.verbosity = verbosity;
		f_switch(inp,os);
		
	} catch(const std::exception& e) {
		QUIT__(f_input(),e.what());
	}
	break;
	case 'p':
	try {
		p_input inp = c.opt.empty() ? p_input(): parseFile<p_input>(c.opt,os);
		if (direct_flag!=-1) inp.direct = direct_flag;
		if (strip_flag!=-1) inp.strip = strip_flag;
		if (!pscin.empty()) inp.pscin = pscin;
		if (!pscout.empty()) inp.pscout = pscout;
		if (!r.empty()) inp.r = r;
		if (!prefix.empty()) inp.prefix = prefix;
		if (verbosity!=NPOS__) inp.verbosity = verbosity;
		p_switch(inp,os);

	} catch(const std::exception& e) {
		QUIT__(p_input(),e.what());
	}
	break;
	case 'i':
	try {
		omen::writeInput(parseFile<h_input>(c.opt.empty() ? "winput": c.opt,os),os);
	} catch(const std::exception& e) {
		QUIT__(h_input(),e.what());
	}
	break;
	case 'h':
	try {
		h_switch(parseFile<h_input>(c.opt.empty() ? "winput": c.opt,os),os);
	} catch(const std::exception& e) {
		QUIT__(h_input(),e.what());
	}
	break;
	case 'w':
	try {
		w_input inp = c.opt.empty() ? w_input():
			aux::parseFile<w_input>(c.opt);
		if (!prefix.empty()) inp.prefix = prefix;
		if (verbosity!=NPOS__) inp.verbosity = verbosity;
		if (!wbh.empty()) inp.wbh = wbh;

		// get fermi energy
		const double Ef = inp.outcar.empty() ? inp.Ef: readEf(inp.outcar);

		// generate the wbh
		ll_hbonds W(inp,os); W.writeToFile(inp.wbh,Ef);

		// run BS tests
		if (!std::isnan(Ef)) {
			const auto E = findBandEdges(Ef,readEig(inp.weig));
			meshBStest(W,inp.C,W.r(),E,inp,os);
			traceBStest(W,inp.C,W.r(),inp,os);
		}
	} catch(const std::exception& e) {
		os << e.what() << "\n";
	}
	break;
	case 'q':
	try {
		if (!all_range_flag)
			os << extractCell(c.opt.empty() ? WBH__: c.opt).
				print(direct_flag) << "\n";
		else {
			const ll_hbonds W(c.opt.empty() ? WBH__: c.opt);

			// find all R vecs for each type
			std::vector<fMat> RR(W.cell().N(),fMat(W.dim(),0));
			for (const auto& i: W) {
				const fMat& cR1 = RR[i.i2()], cR2 = i.R();
				
				fMat nR(W.dim(),cR1.N()+cR2.N());
				nR.resize(distance(nR.cBegin(),
					std::set_union(cR1.ccBegin(),cR1.ccEnd(),
						       cR2.ccBegin(),cR2.ccEnd(),
						        nR.cBegin())));
				nR.shrink_to_fit();
				RR[i.i2()] = std::move(nR);
			}

			// print as POSCAR
			os << "all range contents\n1.0\n"
			   << T(W.cell().B()) << "\n"
			   << W.cell().id() << "\n";
			for (const auto& i: RR)
				os << i.N() << " ";
			
			fMat Ap(W.dim(),0); Ap.reserve(
				std::accumulate(RR.cbegin(),RR.cend(),size_t(0),
					[](const size_t s, const auto& i){ return s+i.size(); }));
			for (size_t i=0; i!=RR.size(); ++i) {
				const fMat p = W.cell().cAt(i);
				for (auto j=RR[i].ccBegin(),e=RR[i].ccEnd(); j!=e; ++j)
					Ap.push_back(p + *j);
			}
			os << "\n" << (direct_flag ? "Direct": "Cartesian")
			   << "\n" << (direct_flag ? T(Ap): T(W.cell().B().prod(Ap))) << "\n";
		}
	} catch(const std::exception& e) {
		os << e.what() << "\n";
	}
	break;
	case 'c':
	try {
		const c_input inp = aux::parseFile<c_input>(c.opt);
		c_switch(inp);
	} catch(const std::exception& e) {
		os << e.what() << "\n";
	}
	break;
	}
}
