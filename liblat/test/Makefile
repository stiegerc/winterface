include def.mk

.PHONY: all run clean

SUITE		= suite
SRCS		= $(SUITE).cc \
		  test_aux_sort.cc \
		  test_cell_ctor.cc test_cell_assign.cc \
		  test_cell_iterators.cc \
		  test_cell_data_access.cc test_cell_memory_management.cc \
		  test_cell_information.cc test_cell_type_information.cc \
		  test_cell_modification.cc \
		  test_cell_conversion.cc \
		  test_cell_comparison.cc \
		  test_io.cc test_hio.cc \
		  test_fn_basis_finding.cc test_fn_metrics_and_clustering.cc \
		  test_fn_linear_path_generation.cc test_fn_generate_hamiltonian.cc \
		  test_fn_wannier_matching.cc test_fn_wannier_tools.cc \
		  test_fn_bandstructure_calc.cc \
		  test_compound_all.cc \
		  test_bonds_all.cc test_hbonds_all.cc test_hbondss_all.cc \
		  test_omen_prepper_all.cc test_BStest_all.cc \
		  test_mesh_all.cc \
		  test_parser_all.cc

MATLIBDIR	= $(shell pwd)/libmat/$(DBGDIR)
LATLIBDIR	= $(shell pwd)/../$(DBGDIR)
LCPPUNIT	= -L$(CPPUNITHOME)/lib -lcppunit
LLAT		= -L$(LATLIBDIR) -Wl,-rpath=$(LATLIBDIR) -llat
LMAT		= -L$(MATLIBDIR) -Wl,-rpath=$(MATLIBDIR) -lmat
LSHARED		= $(LDFLAGS) $(LCPPUNIT) $(LMAT) $(LLAT) $(LMKL) $(LINTEL)
CXXFLAGS	+= -I$(INCLDIR) -I$(SRCDIR) -I../$(INCLDIR) -Ilibmat/$(DBGDIR)/$(INCLDIR)

DEPDIR		:= .d
DEPFLAGS	= -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE		= $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS)
POSTCOMPILE	= mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d
LINK		= $(CXX) $(LSHARED) 
#LINK		= $(CXX) $(LDFLAGS) $(LCPPUNIT) $(LMAT) $(LMKL) $(LINTEL)


all: $(SUITE)

run: all
	./$(SUITE)

veryclean: clean

clean:
	$(RM) $(BINDIR)/*.o
	$(RM) $(DEPDIR)/*.d
	$(RM) $(SUITE)
	$(RM) outp/*


$(SUITE): $(addprefix $(BINDIR)/, $(SRCS:.cc=.o))
	$(LINK) $^ -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.cc $(DEPDIR)/%.d
	$(COMPILE) -c $< -o $@
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d


-include $(addprefix $(DEPDIR)/, $(SRCS:.cc=.d))
