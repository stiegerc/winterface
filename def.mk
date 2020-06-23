# directory structure
INCLDIR 	= include
SRCDIR 		= src
FPICDIR		= fpic
DBGDIR		= debug
RLSDIR		= release
BINDIR		= bin
TESTDIR		= test
DOCDIR		= doc

# misc
RM		= rm -rf
CP		= cp -f
AR		= ar rvs
DOXYGEN		= doxygen


# environment, uncomment one
#ENV = iis-ethz
ENV = home


ifeq ($(ENV),home)
	# blas lapack
	LBLAS		= -lblas
	LLAPACK		= -llapack

	# compiler
	CXX		= g++

	# linker
	LD		= $(CXX)
	
	# flags
	CXXFLAGS	+= -O2 -std=c++17 -fopenmp -DDOUBLE__
	CXXFLAGS	+= -Wall -Wextra -Werror -Wno-unused-parameter
	CXXFLAGS	+= -Wno-virtual-move-assign
	#CXXFLAGS	+= -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS 	+= -lpthread -ldl -lgomp $(LBLAS) $(LLAPACK) -lcrypto
	#LDFLAGS 	+= $(CXXFLAGS)
endif

ifeq ($(ENV),iis-ethz)
	# cppunit
	CPPUNITHOME	= /home/stiegerc/.system/

	# openssl
	SSLHOME		= /usr/include/

	# intel
	INTELROOT	= /usr/pack/intel_compiler-2015.0.090-af/Linux-x86_64
	INTELHOME	= $(INTELROOT)/compiler/lib/intel64
	INTELLIBS	= imf intlc iomp5
	LINTEL		= -L$(INTELHOME) -Wl,-rpath=$(INTELHOME) -limf -lintlc -liomp5
	LINTELSTATIC	= $(addprefix $(INTELHOME)/lib, imf.a iomp5.a irc.a)
	MKLHOME		= $(INTELROOT)/mkl/lib/intel64
	MKLLIBS		= mkl_intel_lp64 mkl_core mkl_sequential
	LMKL		= -L$(MKLHOME) -Wl,-rpath=$(MKLHOME) $(addprefix -l, $(MKLLIBS))
	LMKLSTATIC	= -Wl,--start-group $(addprefix $(MKLHOME)/lib, $(MKLLIBS:=.a)) -Wl,--end-group

	# blas lapack
	LBLAS		= $(LMKL) $(LINTEL)
	#LBLAS		= -L/home/stiegerc/.system/lib64 -lblas
	#LLAPACK		= -llapack

	# open MP
	LGOMP_STATIC	= /home/stiegerc/.system/lib64/libgomp.a

	# compiler
	CXX		= g++-8.2.0-af

	# linker
	LD		= $(CXX)
	LD_STATIC	= icpc-2015.0.090-af -static-intel -mkl=sequential

	# flags
	LDFLAGS 	+= -lpthread -ldl -lgomp $(LBLAS) $(LLAPACK) -L$(CPPUNITHOME)/lib/ -lcrypto
	LDFLAGS_STATIC  += -static-libgstdc++ -static-libgcc -lpthread -lm -lrt $(LGOMP_STATIC) -lcrypto
	CXXFLAGS	+= -O2 -std=c++14 -Wall -fopenmp -I$(CPPUNITHOME)/include -I$(SSLHOME)
	#CXXFLAGS	+= -g -O0 -Q -std=c++14 -Wall -fopenmp
endif
