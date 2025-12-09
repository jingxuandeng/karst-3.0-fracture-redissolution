

#MUMPS linking
MDIR    = /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/MUMPS_5.6.0

BREW=/usr/local/Cellar

CC	= clang       
CXX 	= clang++                                            
FL	= gfortran   
AR      = ar vr 

# LSCOTCH   =  -L /Users/jingxuandeng/phd/KRG/research/network_dissol_project/karst_code/Libraries/scotch/lib -lesmumps -lscotch -lscotcherr
LSCOTCHDIR = $(BREW)/scotch/7.0.7/lib
LSCOTCH   = -L$(LSCOTCHDIR) -lesmumps -lscotch -lscotcherr
LMETIS    =  -L $(BREW)/metis/5.1.0/lib -lmetis
LPORD = -L $(MDIR)/lib -lpord

# Corresponding variables reused later
LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)

# LLAPACKDIR = $(BREW)/lapack/3.12.0/lib
# LAPACK =  -L $(LLAPACKDIR) -llapack
LAPACK = -L /usr/local/opt/lapack/lib  -llapack

LIBSEQ  = $(LAPACK) -L$(MDIR)/libseq -lmpiseq

# LBLASDIR = $(BREW)/openblas/0.3.29/lib
# LIBBLAS = -L $(LBLASDIR) -lblas
LIBBLAS = -L /usr/local/opt/openblas/lib -lblas

LIBOTHERS = -lpthread
LFORTRAN = -L /usr/local/opt/gcc/lib/gcc/15  -lgfortran
LIBMUMPS = -L $(MDIR)/lib/ -ldmumps -L $(MDIR)/lib/ -lmumps_common



MINCLUDES = -I $(MDIR)/include -I $(MDIR)/libseq 
ALGDIR  = ../algorithms
ARMINCDIR = -I ../include -I $(ALGDIR)

COPTS	= -Wall -pedantic -O2 
CFLAGS	= $(COPTS) $(MINCLUDES) $(ARMINCDIR)

CXXOPTS = -std=c++17 -Weffc++ -Wall -Wextra -Wfatal-errors -Wunused -pedantic -w -Winline -fno-rtti -ggdb -D_GLIBCXX_DEBUG -lstdc++ $(LSFML) $(LFORTRAN)
CXXFLAGS = $(CXXOPTS) 


MLIBS	= $(LIBMUMPS) $(LORDERINGS) $(LIBSEQ) $(LIBBLAS) $(LIBOTHERS) -lstdc++
MFLAGS	= $(COPTS) $(MLIBS) 

TRICK	= 2>&1 >/dev/null | grep -v -e '^/var/folders/*' -e '^[[:space:]]*\.section' -e '^[[:space:]]*\^[[:space:]]*~*'


OBJ     = network.o main.o printing.o printing_net.o node.o pore.o grain.o merging.o initial_modifications.o\
evolution.o network_maintenance.o setup.o dissolution.o control_and_verification.o import_export.o network_generation.o\
tests.o save_in_VTK.o analysis.o pattern_analysis.o
INC     = network.h node.h pore.h grain.h printing.h constants.h algorithms_cc.h


ALGLIBS      = -lm -L $(ALGDIR) -lalgorithms
LFLAGS    = $(ALGLIBS) $(MFLAGS) 


network_generation.o: network_generation.cc
	$(CXX) $(CXXFLAGS) $(ARMINCDIR) -c network_generation.cc
.c.o: 
	$(CC)  $(CFLAGS) -c $< 

.cc.o: 
#	$(CXX) $(CXXFLAGS) $(CINCLUDES) -c $< 
	$(CXX)  $(CFLAGS) -c $< 


karst: $(OBJ) 	
	$(CXX) -v -o ../bin/karst $^ $(LFLAGS) $(CXXFLAGS)

algorithms:
	$(MAKE) all -C $(ALGDIR)

clean:
	rm -f *.o


all: 
	mkdir -p ../bin
	make algorithms
	make karst

