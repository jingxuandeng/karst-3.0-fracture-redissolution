UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC		= clang        
	CXX 	= clang++                                            
	FL		= gfortran 
	MDIR	=  ~/.nix-profile/opt/mumps
	FLIBS	= -L ~/.nix-profile/lib -lgfortran -lquadmath -lblas
else 
	CC      = /opt/homebrew/bin/gcc-12
	CXX 	= /opt/homebrew/bin/g++-12
	FL      = /opt/homebrew/bin/gfortran
	MDIR	= ~/Desktop/KARST/MUMPS_5.6.0
	FLIBS	= -L/opt/homebrew/lib -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/12 -lblas -llapack -lesmumps -lscotch -lscotcherr -lmumps_common -lmetis -lgomp -lgfortran
endif	

#MUMPS linking
MLIBS	= -ldmumps -lmpiseq -lpthread -lpord  -lstdc++ $(FLIBS)
ALGDIR  = ../algorithms

COPTS	= -Wall -pedantic -O2  
CINCLUDES	= -I ../include -I ../algorithms 
MINCLUDES	= -I ../include -I $(MDIR)/include -I $(MDIR)/libseq

CFLAGS	 = $(COPTS) $(CINCLUDES) $(MINCLUDES) 
CXXFLAGS = -O2 -std=c++11 -Weffc++ -Wall -Wextra -Wfatal-errors -Wunused -pedantic -w -Winline -O2
TMP      = -fno-rtti -ggdb -D_GLIBCXX_DEBUG -lsfml-graphics -lsfml-window -lsfml-system


OBJ	= network.o main.o printing.o printing_net.o node.o pore.o grain.o merging.o initial_modifications.o\
evolution.o network_maintenance.o setup.o dissolution.o control_and_verification.o import_export.o network_generation.o\
tests.o analysis.o pattern_analysis.o save_in_VTK.o
INC	= network.h node.h pore.h grain.h printing.h constants.h algorithms_cc.h


LINKOPTS  = $(COPTS) -lstdc++ $(FLIBS)  -lgfortran
LIBS      = -lm -L ../algorithms -lalgorithms
MFLAGS    = $(MLIBS) -L $(MDIR)/lib -L $(MDIR)/libseq
LFLAGS    = $(LINKOPTS) $(LIBS) $(MFLAGS) 


## Default target
#.DEFAULT:
#    @echo "No target specified. To build a program, use: make <program_name>"

.c.o:
	$(CC) $(CFLAGS) -c $< 

.cc.o: $(INC)
	$(CXX) $(CXXFLAGS) $(CINCLUDES) -c $<
   	 

karst: $(OBJ)
	$(CXX)  -o ../bin/karst $^ $(LFLAGS) $(CXXFLAGS)

algorithms:
	$(MAKE) all -C $(ALGDIR)

clean:
	rm -f *.o

clean_alg:
	rm -f ../algorithms/*.o ../algorithms/*.a

all:
	mkdir -p ../bin
	make algorithms
	make karst

default: all

