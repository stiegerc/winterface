// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

#include <cppunit/ui/text/TestRunner.h>
#include "lm_defs.h"
#include "test_mtol.h"
#include "test_tItr_all.h"
#include "test_c_tItr_all.h"
#include "test_cpxItr_all.h"
#include "test_c_tVecItr_all.h"
#include "test_tVecItr_all.h"
#include "test_cr_tVecItr_all.h"
#include "test_r_tVecItr_all.h"
#include "test_tCol_all.h"
#include "test_tRow_all.h"
#include "test_tMat_ctor.h"
#include "test_tMat_assign.h"
#include "test_tMat_data_access_itrs.h"
#include "test_tMat_memory_management.h"
#include "test_tMat_basic_properties_info.h"
#include "test_tMat_basic_modification.h"
#include "test_tMat_modification.h"
#include "test_tMat_conversion.h"
#include "test_tMat_logical.h"
#include "test_tMat_comparison.h"
#include "test_tMat_el_arithmetic.h"
#include "test_tMat_matrix_arithmetic.h"
#include "test_tFn_math_functions.h"
#include "test_fn_eigenvalue_computation.h"
#include "test_tFn_orth.h"
#include "test_tFn_mat_gen.h"
#include "test_tFn_functionals.h"
#include "test_fn_dot_products.h"
#include "test_fn_vector_sums.h"
#include "test_fn_R3_only.h"
#include "test_tFn_conversion.h"
#include "test_tFn_information.h"
#include "test_tFn_comparison.h"


int main(int argc, char** argv) {
	CppUnit::TextUi::TestRunner runner;

	// mtol
	runner.addTest(test_mtol::suite());

	// iterators
	runner.addTest(test_c_tItr_all<RE__>::suite());
	runner.addTest(test_c_tItr_all<CPX__>::suite());
	runner.addTest(test_tItr_all<RE__>::suite());
	runner.addTest(test_tItr_all<CPX__>::suite());
	runner.addTest(test_cpxItr_all<0>::suite());
	runner.addTest(test_cpxItr_all<1>::suite());
	runner.addTest(test_c_tVecItr_all<RE__,RE__,CPX__,lm_tRow<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_c_tVecItr_all<CPX__,RE__,CPX__,lm_tRow<CPX__,RE__,CPX__>>::suite());
	runner.addTest(test_c_tVecItr_all<RE__,RE__,CPX__,lm_tCol<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_c_tVecItr_all<CPX__,RE__,CPX__,lm_tCol<CPX__,RE__,CPX__>>::suite());
	runner.addTest(test_tVecItr_all<RE__,RE__,CPX__,lm_tRow<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_tVecItr_all<CPX__,RE__,CPX__,lm_tRow<CPX__,RE__,CPX__>>::suite());
	runner.addTest(test_tVecItr_all<RE__,RE__,CPX__,lm_tCol<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_tVecItr_all<CPX__,RE__,CPX__,lm_tCol<CPX__,RE__,CPX__>>::suite());
	runner.addTest(test_cr_tVecItr_all<RE__,RE__,CPX__,lm_tRow<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_cr_tVecItr_all<CPX__,RE__,CPX__,lm_tRow<CPX__,RE__,CPX__>>::suite());
	runner.addTest(test_cr_tVecItr_all<RE__,RE__,CPX__,lm_tCol<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_cr_tVecItr_all<CPX__,RE__,CPX__,lm_tCol<CPX__,RE__,CPX__>>::suite());
	runner.addTest(test_r_tVecItr_all<RE__,RE__,CPX__,lm_tRow<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tRow<CPX__,RE__,CPX__>>::suite());
	runner.addTest(test_r_tVecItr_all<RE__,RE__,CPX__,lm_tCol<RE__,RE__,CPX__>>::suite());
	runner.addTest(test_r_tVecItr_all<CPX__,RE__,CPX__,lm_tCol<CPX__,RE__,CPX__>>::suite());

	// tVec
	runner.addTest(test_tCol_all<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tCol_all<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tRow_all<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tRow_all<CPX__,RE__,CPX__>::suite());

	// tMat
	runner.addTest(test_tMat_ctor<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_ctor<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_assign<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_assign<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_data_access_itrs<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_data_access_itrs<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_memory_management<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_memory_management<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_basic_properties_info<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_basic_properties_info<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_basic_modification<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_basic_modification<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_modification<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_modification<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_conversion<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_conversion<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_logical<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_logical<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_comparison<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_comparison<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_el_arithmetic<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_el_arithmetic<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_matrix_arithmetic<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tMat_matrix_arithmetic<CPX__,RE__,CPX__>::suite());

	// fn
	runner.addTest(test_tFn_math_functions<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_math_functions<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_fn_eigenvalue_computation::suite());
	runner.addTest(test_tFn_orth<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_orth<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_mat_gen<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_mat_gen<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_functionals<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_functionals<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_fn_dot_products::suite());
	runner.addTest(test_fn_vector_sums::suite());
	runner.addTest(test_tFn_conversion<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_conversion<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_information<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_information<CPX__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_comparison<RE__,RE__,CPX__>::suite());
	runner.addTest(test_tFn_comparison<CPX__,RE__,CPX__>::suite());
	
	// R^3 only
	runner.addTest(test_fn_R3_only::suite());

	if (argc==1)
//		for (size_t i=0; i!=1000; ++i)
		runner.run("",false,true,false);
	else
		for (int i=1; i<argc; ++i)
			runner.run(argv[i],false,true,false);

	return 0;
}
