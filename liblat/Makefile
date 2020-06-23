include def.mk

.PHONY: default all clean veryclean
.PHONY: release_fpic release_lib release_headers release
.PHONY: debug_fpic debug_lib debug_headers debug

LIBMAT		= libmat
LIBNAME		= lat

SRCS		= ll_types.cc \
		  ll_cell.cc \
		  ll_io.cc ll_fn.cc \
		  ll_hio.cc \
		  ll_bonds.cc ll_hbonds.cc ll_hbondss.cc \
		  aux_parser.cc \
		  ll_mesh.cc \
		  ll_omen.cc ll_BStest.cc
HDRS		= ll_defs.h \
		  ll_types.h ll_compound.h \
		  ll_cell.h \
		  ll_io.h ll_fn.h \
		  ll_hio.h \
		  ll_bonds.h ll_hbonds.h ll_hbondss.h \
		  aux_parser.h \
		  ll_sparse.h ll_mesh.h \
		  ll_omen.h ll_BStest.h



#SRCS		= ll_types.cc \
		  ll_cell.cc \
		  ll_io.cc ll_fn.cc \
		  ll_Hcont.cc \
		  ll_input.cc \
		  ll_bonds.cc \
		  aux_parser.cc \
		  omen_input.cc omen_io.cc omen_binding.cc \
		  ll_mesh.cc
#HDRS		= ll_defs.h \
		  ll_types.h ll_compound.h \
		  ll_cell.h \
		  ll_io.h ll_fn.h \
		  ll_Hcont.h \
		  ll_input.h \
		  ll_bonds.h \
		  ll_sparse.h ll_mesh.h \
		  aux_parser.h \
		  omen_input.h omen_io.h omen_binding.h

OBJ		= $(SRCS:.cc=.o)
OBJ_DBG		= $(addprefix dbg_, $(OBJ))
OBJ_DBG_FPIC	= $(addprefix dbg_fpic_, $(OBJ))
OBJ_RLS		= $(addprefix rls_, $(OBJ))
OBJ_RLS_FPIC	= $(addprefix rls_fpic_, $(OBJ))

CXXFLAGS	+= -I$(INCLDIR) -I$(SRCDIR)
RLSFLAGS	= $(CXXFLAGS) $(CPPFLAGS) -I$(LIBMAT)/$(RLSDIR)/$(INCLDIR) -DNDEBUG
DBGFLAGS	= $(CXXFLAGS) $(CPPFLAGS) -I$(LIBMAT)/$(DBGDIR)/$(INCLDIR)

DEPDIR		:= .d
DEPFLAGS	= -MT $@ -MMD -MP -MF

default: debug_lib debug_fpic debug_headers

all: release debug

clean:
	$(RM) $(BINDIR)/*.o
	$(RM) $(DEPDIR)/*.d
	$(RM) $(DBGDIR)/lib$(LIBNAME).a
	$(RM) $(DBGDIR)/lib$(LIBNAME).so
	$(RM) $(DBGDIR)/$(INCLDIR)/*.h

veryclean: clean
	$(RM) $(RLSDIR)/lib$(LIBNAME).a
	$(RM) $(RLSDIR)/lib$(LIBNAME).so
	$(RM) $(RLSDIR)/$(INCLDIR)/*.h
	$(MAKE) -C test clean


release: release_lib release_fpic release_headers

debug: debug_lib debug_fpic debug_headers


release_lib: $(RLSDIR)/lib$(LIBNAME).a

release_fpic: $(RLSDIR)/lib$(LIBNAME).so

release_headers:
	$(CP) $(addprefix $(INCLDIR)/, $(HDRS)) $(RLSDIR)/$(INCLDIR)


debug_lib: $(DBGDIR)/lib$(LIBNAME).a

debug_fpic: $(DBGDIR)/lib$(LIBNAME).so

debug_headers:
	$(CP) $(addprefix $(INCLDIR)/, $(HDRS)) $(DBGDIR)/$(INCLDIR)


$(RLSDIR)/lib$(LIBNAME).a: $(addprefix $(BINDIR)/, $(OBJ_RLS))
	$(AR) $@  $^

$(RLSDIR)/lib$(LIBNAME).so: $(addprefix $(BINDIR)/, $(OBJ_RLS_FPIC))
	$(CXX) $^ -shared -o $@

$(DBGDIR)/lib$(LIBNAME).a: $(addprefix $(BINDIR)/, $(OBJ_DBG))
	$(AR) $@  $^

$(DBGDIR)/lib$(LIBNAME).so: $(addprefix $(BINDIR)/, $(OBJ_DBG_FPIC))
	$(CXX) $^ -shared -o $@


$(BINDIR)/rls_%.o: $(SRCDIR)/%.cc $(DEPDIR)/rls_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/rls_$*.Td $(RLSFLAGS) -c $< -o $@
	mv -f $(DEPDIR)/rls_$*.Td $(DEPDIR)/rls_$*.d

$(BINDIR)/rls_fpic_%.o: $(SRCDIR)/%.cc $(DEPDIR)/rls_fpic_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/rls_fpic_$*.Td $(RLSFLAGS) -fPIC -c $< -o $@
	mv -f $(DEPDIR)/rls_fpic_$*.Td $(DEPDIR)/rls_fpic_$*.d

$(BINDIR)/dbg_%.o: $(SRCDIR)/%.cc $(DEPDIR)/dbg_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/dbg_$*.Td $(DBGFLAGS) -c $< -o $@
	mv -f $(DEPDIR)/dbg_$*.Td $(DEPDIR)/dbg_$*.d

$(BINDIR)/dbg_fpic_%.o: $(SRCDIR)/%.cc $(DEPDIR)/dbg_fpic_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/dbg_fpic_$*.Td $(DBGFLAGS) -fPIC -c $< -o $@
	mv -f $(DEPDIR)/dbg_fpic_$*.Td $(DEPDIR)/dbg_fpic_$*.d


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d


-include $(addprefix $(DEPDIR)/, $(OBJ_DBG:.o=.d))
-include $(addprefix $(DEPDIR)/, $(OBJ_DBG_FPIC:.o=.d))
-include $(addprefix $(DEPDIR)/, $(OBJ_RLS:.o=.d))
-include $(addprefix $(DEPDIR)/, $(OBJ_RLS_FPIC:.o=.d))
