include def.mk

.PHONY: default all clean veryclean
.PHONY: release_fpic release_lib release_headers release
.PHONY: debug_fpic debug_lib debug_headers debug

TOOL		= ltool
TOOL_DBG	= dbg_$(TOOL)
TOOL_RLS	= rls_$(TOOL)
TOOL_STATIC	= stat_$(TOOL)

LIBMAT		= libmat
LIBLAT		= liblat
SRCS		= b_switch.cc l_switch.cc f_switch.cc\
		  p_switch.cc c_switch.cc h_switch.cc\
		  ltool.cc

HDRS		= b_switch.h l_switch.h f_switch.h\
		  p_switch.h h_switch.h c_switch.h


LIBMATDIR_DBG	= $(shell pwd)/$(LIBMAT)/$(DBGDIR)
LIBLATDIR_DBG	= $(shell pwd)/$(LIBLAT)/$(DBGDIR)
LIBMATDIR_RLS	= $(shell pwd)/$(LIBMAT)/$(RLSDIR)
LIBLATDIR_RLS	= $(shell pwd)/$(LIBLAT)/$(RLSDIR)

LSHARED		= $(LDFLAGS)
LLAT_DBG	= -L$(LIBLATDIR_DBG) -Wl,-rpath=$(LIBLATDIR_DBG) -llat
LMAT_DBG	= -L$(LIBMATDIR_DBG) -Wl,-rpath=$(LIBMATDIR_DBG) -lmat
LLAT_RLS	= -L$(LIBLATDIR_RLS) -Wl,-rpath=$(LIBLATDIR_RLS) -llat
LMAT_RLS	= -L$(LIBMATDIR_RLS) -Wl,-rpath=$(LIBMATDIR_RLS) -lmat
LINK_DBG	= $(LD) $(LSHARED) $(LMAT_DBG) $(LLAT_DBG)
LINK_RLS	= $(LD) $(LSHARED) $(LMAT_RLS) $(LLAT_RLS)
LINK_STATIC	= $(LD_STATIC) $(LIBMATDIR_DBG)/$(LIBMAT).a $(LIBLATDIR_DBG)/$(LIBLAT).a

OBJ		= $(SRCS:.cc=.o)
OBJ_DBG		= $(addprefix dbg_, $(OBJ))
OBJ_RLS		= $(addprefix rls_, $(OBJ))

CXXFLAGS	+= -I$(INCLDIR) -I$(LIBLAT)/$(INCLDIR)
RLSFLAGS	= $(CXXFLAGS) $(CPPFLAGS) -I$(LIBMAT)/$(RLSDIR)/$(INCLDIR) -DNDEBUG
DBGFLAGS	= $(CXXFLAGS) $(CPPFLAGS) -I$(LIBMAT)/$(DBGDIR)/$(INCLDIR)

DEPDIR		:= .d
DEPFLAGS	= -MT $@ -MMD -MP -MF

default: debug

all: release debug

clean:
	$(RM) $(BINDIR)/*.o
	$(RM) $(DEPDIR)/*.d
	$(RM) $(TOOL_DBG)

veryclean: clean
	$(RM) $(TOOL_RLS)
	$(RM) $(TOOL_STATIC)


debug: $(TOOL_DBG)

static: $(TOOL_STATIC)

release: $(TOOL_RLS)
	

$(TOOL_DBG): $(addprefix $(BINDIR)/, $(OBJ_DBG)) \
	$(LIBMATDIR_DBG)/$(LIBMAT).so $(LIBLATDIR_DBG)/$(LIBLAT).so
	$(LINK_DBG) $^ -o $@

$(TOOL_RLS): $(addprefix $(BINDIR)/, $(OBJ_RLS)) \
	$(LIBMATDIR_RLS)/$(LIBMAT).so $(LIBLATDIR_RLS)/$(LIBLAT).so
	$(LINK_RLS) $^ -o $@

$(TOOL_STATIC): $(addprefix $(BINDIR)/, $(OBJ_DBG)) \
	$(LIBMATDIR_DBG)/$(LIBMAT).a $(LIBLATDIR_DBG)/$(LIBLAT).a
	$(LINK_STATIC) $^ $(LDFLAGS_STATIC) -o $@


$(BINDIR)/rls_%.o: $(SRCDIR)/%.cc $(DEPDIR)/rls_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/rls_$*.Td $(RLSFLAGS) -c $< -o $@
	mv -f $(DEPDIR)/rls_$*.Td $(DEPDIR)/rls_$*.d

$(BINDIR)/dbg_%.o: $(SRCDIR)/%.cc $(DEPDIR)/dbg_%.d
	$(CXX) $(DEPFLAGS) $(DEPDIR)/dbg_$*.Td $(DBGFLAGS) -c $< -o $@
	mv -f $(DEPDIR)/dbg_$*.Td $(DEPDIR)/dbg_$*.d


$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d


-include $(addprefix $(DEPDIR)/, $(OBJ_DBG:.o=.d))
-include $(addprefix $(DEPDIR)/, $(OBJ_RLS:.o=.d))
