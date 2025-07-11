
#MUMPS linking
MDIR    = /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/MUMPS_5.6.0

CC	= clang       
CXX 	= clang++                                            
FL	= gfortran   
AR      = ar vr

BREW=/usr/local/Cellar

LSCOTCH   =  -L /Users/jingxuandeng/phd/KRG/research/network_dissol_project/karst_code/Libraries/scotch/lib -lesmumps -lscotch -lscotcherr
LMETIS    =  -L $(BREW)/metis/5.1.0/lib -lmetis
LPORD = -L $(MDIR)/lib -lpord

# Corresponding variables reused later
LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)

LLAPACKDIR = $(BREW)/lapack/3.12.0/lib
LAPACK =  -L $(LLAPACKDIR) -llapack

LIBSEQ  = $(LAPACK) -L$(MDIR)/libseq -lmpiseq

LBLASDIR = $(BREW)/openblas/0.3.29/lib
LIBBLAS = -L $(LBLASDIR) -lblas

LIBOTHERS = -lpthread

LFORTRAN = -L /usr/local/opt/gcc/lib/gcc/15  -lgfortran
LIBMUMPS = -L $(MDIR)/lib/ -ldmumps -L $(MDIR)/lib/ -lmumps_common

MINCLUDES = -I $(MDIR)/include -I $(MDIR)/libseq 

COPTS	= -Wall -pedantic -O2 
CFLAGS	= $(COPTS) $(MINCLUDES) 

CXXOPTS = -std=c++11 -Weffc++ -Wall -Wextra -Wfatal-errors -Wunused -pedantic -w -Winline -fno-rtti -ggdb -D_GLIBCXX_DEBUG $(LSFML) $(LFORTRAN)
CXXFLAGS = $(CXXOPTS) $(MINCLUDE)


MLIBS	= $(LIBMUMPS) $(LORDERINGS) $(LIBSEQ) $(LIBBLAS) $(LIBOTHERS) -lstdc++
MFLAGS	= $(COPTS) $(MLIBS) 


OBJ	= matrix_solvers.o triangularization.o main.o
INC	= algorithms.h 

.cc.o: 
	$(CXX) $(CXXFLAGS) -c $< 

.c.o: 
	$(CC) $(CFLAGS) -c $< 


algorithms: $(OBJ)	
	#$(FL) -o ../bin/algorithms $^ $(MFLAGS)
	$(AR) libalgorithms.a matrix_solvers.o triangularization.o 
	$(CXX) -v -o ../bin/algorithms $^ $(MFLAGS) $(CXXFLAGS)


clean:
	rm -f *.o *.a


all: 
	mkdir -p ../bin
	make algorithms
