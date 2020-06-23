include def.mk

.PHONY: all run clean

SUITE		= suite
SRCS		= $(SUITE).cc \
		  test_mtol.cc \
		  test_c_fItr_all.cc test_c_cItr_all.cc \
		  test_fItr_all.cc test_cItr_all.cc \
		  test_reItr_all.cc test_imItr_all.cc \
		  test_c_fRowItr_all.cc test_c_cRowItr_all.cc \
		  test_c_fColItr_all.cc test_c_cColItr_all.cc \
		  test_fRowItr_all.cc test_cRowItr_all.cc \
		  test_fColItr_all.cc test_cColItr_all.cc \
		  test_cr_fRowItr_all.cc test_cr_cRowItr_all.cc \
		  test_cr_fColItr_all.cc test_cr_cColItr_all.cc \
		  test_r_fRowItr_all.cc test_r_cRowItr_all.cc \
		  test_r_fColItr_all.cc test_r_cColItr_all.cc \
		  test_fCol_all.cc test_cCol_all.cc \
		  test_fRow_all.cc test_cRow_all.cc \
		  test_fMat_ctor.cc test_cMat_ctor.cc \
		  test_fMat_assign.cc test_cMat_assign.cc \
		  test_fMat_data_access_itrs.cc test_cMat_data_access_itrs.cc \
		  test_fMat_memory_management.cc test_cMat_memory_management.cc \
		  test_fMat_basic_properties_info.cc test_cMat_basic_properties_info.cc \
		  test_fMat_basic_modification.cc test_cMat_basic_modification.cc \
		  test_fMat_modification.cc test_cMat_modification.cc \
		  test_fMat_conversion.cc test_cMat_conversion.cc \
		  test_fMat_logical.cc test_cMat_logical.cc \
		  test_fMat_comparison.cc test_cMat_comparison.cc \
		  test_fMat_el_arithmetic.cc test_cMat_el_arithmetic.cc \
		  test_fMat_matrix_arithmetic.cc test_cMat_matrix_arithmetic.cc \
		  test_fFn_math_functions.cc test_cFn_math_functions.cc \
		  test_fn_eigenvalue_computation.cc \
		  test_fFn_orth.cc test_cFn_orth.cc \
		  test_fFn_mat_gen.cc test_cFn_mat_gen.cc \
		  test_fFn_functionals.cc test_cFn_functionals.cc \
		  test_fn_dot_products.cc \
		  test_fn_vector_sums.cc \
		  test_fn_R3_only.cc \
		  test_fFn_conversion.cc test_cFn_conversion.cc \
		  test_fFn_information.cc test_cFn_information.cc \
		  test_fFn_comparison.cc test_cFn_comparison.cc

LIBDIR		= $(shell pwd)/../$(DBGDIR)
LCPPUNIT	= -L$(CPPUNITHOME)/lib -lcppunit
LMAT		= -L$(LIBDIR) -Wl,-rpath=$(LIBDIR) -lmat
LSHARED		= $(LDFLAGS) $(LCPPUNIT) $(LMAT)
CXXFLAGS	+= -I$(INCLDIR) -I$(SRCDIR) -I../$(INCLDIR)

DEPDIR		:= .d
DEPFLAGS	= -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE		= $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS)
POSTCOMPILE	= mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d
LINK		= $(CXX) $(LSHARED) 


all: $(SUITE)

run: all
	./$(SUITE)

clean:
	$(RM) $(BINDIR)/*.o
	$(RM) $(DEPDIR)/*.d
	$(RM) $(SUITE)

veryclean: clean

$(SUITE): $(addprefix $(BINDIR)/, $(SRCS:.cc=.o))
	$(LINK) $^ -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.cc $(DEPDIR)/%.d
	$(COMPILE) -c $< -o $@
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d


-include $(addprefix $(DEPDIR)/, $(SRCS:.cc=.d))
