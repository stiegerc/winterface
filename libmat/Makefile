include def.mk

.PHONY: default all clean veryclean
.PHONY: release_fpic release_lib release_headers release
.PHONY: debug_fpic debug_lib debug_headers debug

LIBNAME		= mat

SRCS		= lm_fArray.cc lm_cArray.cc \
		  lm_fRow.cc lm_cRow.cc lm_fCol.cc lm_cCol.cc \
		  lm_fMat.cc lm_cMat.cc lm_fn.cc
HDRS		= lm_tItr.h lm_cpxItr.h lm_tVecItr.h \
		  lm_tArray.h lm_ref_tArray.h \
		  lm_tCol.h lm_tRow.h \
		  lm_tMat.h lm_fn.h \
		  lm_lambda.h \
		  lm_defs.h lm_types.h lm_ops.h lib$(LIBNAME).h

OBJ		= $(SRCS:.cc=.o)
OBJ_DBG		= $(addprefix dbg_, $(OBJ))
OBJ_DBG_FPIC	= $(addprefix dbg_fpic_, $(OBJ))
OBJ_RLS		= $(addprefix rls_, $(OBJ))
OBJ_RLS_FPIC	= $(addprefix rls_fpic_, $(OBJ))

CXXFLAGS	+= -I$(INCLDIR)

DEPDIR		:= .d
DEPFLAGS	= -MT $@ -MMD -MP -MF

default: debug

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
	$(CXX) $(DEPFLAGS) $(DEPDIR)/rls_$*.Td $(CXXFLAGS) $(CPPFLAGS) -DNDEBUG -c $< -o $@
	mv -f $(DEPDIR)/rls_$*.Td $(DEPDIR)/rls_$*.d

$(BINDIR)/rls_fpic_%.o: $(SRCDIR)/%.cc $(DEPDIR)/rls_fpic_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/rls_fpic_$*.Td $(CXXFLAGS) $(CPPFLAGS) -DNDEBUG -fPIC -c $< -o $@
	mv -f $(DEPDIR)/rls_fpic_$*.Td $(DEPDIR)/rls_fpic_$*.d

$(BINDIR)/dbg_%.o: $(SRCDIR)/%.cc $(DEPDIR)/dbg_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/dbg_$*.Td $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
	mv -f $(DEPDIR)/dbg_$*.Td $(DEPDIR)/dbg_$*.d

$(BINDIR)/dbg_fpic_%.o: $(SRCDIR)/%.cc $(DEPDIR)/dbg_fpic_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/dbg_fpic_$*.Td $(CXXFLAGS) $(CPPFLAGS) -fPIC -c $< -o $@
	mv -f $(DEPDIR)/dbg_fpic_$*.Td $(DEPDIR)/dbg_fpic_$*.d


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d


-include $(addprefix $(DEPDIR)/, $(OBJ_DBG:.o=.d))
-include $(addprefix $(DEPDIR)/, $(OBJ_DBG_FPIC:.o=.d))
-include $(addprefix $(DEPDIR)/, $(OBJ_RLS:.o=.d))
-include $(addprefix $(DEPDIR)/, $(OBJ_RLS_FPIC:.o=.d))
