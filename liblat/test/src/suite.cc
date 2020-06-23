// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include <cppunit/ui/text/TestRunner.h>
#include "test_aux_for_each.h"
#include "test_aux_sort.h"
#include "test_compound_all.h"
#include "test_cell_ctor.h"
#include "test_cell_assign.h"
#include "test_cell_iterators.h"
#include "test_cell_data_access.h"
#include "test_cell_memory_management.h"
#include "test_cell_information.h"
#include "test_cell_type_information.h"
#include "test_cell_modification.h"
#include "test_cell_conversion.h"
#include "test_cell_comparison.h"
#include "test_io.h"
#include "test_hio.h"
#include "test_fn_basis_finding.h"
#include "test_fn_metrics_and_clustering.h"
#include "test_fn_linear_path_generation.h"
#include "test_fn_generate_hamiltonian.h"
#include "test_fn_wannier_matching.h"
#include "test_fn_wannier_tools.h"
#include "test_fn_bandstructure_calc.h"
#include "test_bonds_all.h"
#include "test_hbonds_all.h"
#include "test_hbondss_all.h"
#include "test_omen_prepper_all.h"
#include "test_BStest_all.h"
#include "test_mesh_all.h"
#include "test_parser_all.h"


int main(int argc, char** argv) {
	CppUnit::TextUi::TestRunner runner;

	// aux
	runner.addTest(test_aux_sort::suite());

	// compound
	runner.addTest(test_compound_all::suite());	

	// cell
	runner.addTest(test_cell_ctor::suite());
	runner.addTest(test_cell_assign::suite());
	runner.addTest(test_cell_iterators::suite());
	runner.addTest(test_cell_data_access::suite());
	runner.addTest(test_cell_memory_management::suite());
	runner.addTest(test_cell_information::suite());
	runner.addTest(test_cell_type_information::suite());
	runner.addTest(test_cell_modification::suite());
	runner.addTest(test_cell_conversion::suite());
	runner.addTest(test_cell_comparison::suite());

	// io
	runner.addTest(test_io::suite());
	runner.addTest(test_hio::suite());

	// fn
	runner.addTest(test_fn_basis_finding::suite());
	runner.addTest(test_fn_metrics_and_clustering::suite());
	runner.addTest(test_fn_linear_path_generation::suite());
	runner.addTest(test_fn_generate_hamiltonian::suite());
	runner.addTest(test_fn_wannier_matching::suite());
	runner.addTest(test_fn_wannier_tools::suite());
	runner.addTest(test_fn_bandstructure_calc::suite());

	// bonds
	runner.addTest(test_bonds_all::suite());
	runner.addTest(test_hbonds_all::suite());
	runner.addTest(test_hbondss_all::suite());

	// omen prepper and BStest
	runner.addTest(test_omen_prepper_all::suite());
	runner.addTest(test_BStest_all::suite());

	// mesh	
	runner.addTest(test_mesh_all::suite());

	// parser
	runner.addTest(test_parser_all::suite());

	if (argc==1)
//		for (size_t i=0; i!=100000; ++i)
			runner.run("",false,true,false);
	else
		for (int i=1; i<argc; ++i)
			runner.run(argv[i],false,true,false);

	return 0;
}
