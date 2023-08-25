#include "algorithms.h"
 

//solving set of linear equations
int solve_matrix(int W, int R_N, int* w_r, int* w_c, double* B, double* y){

	DMUMPS_STRUC_C	id;		// struktura przez ktora komunikujemy sie z mumpsem 
	
	// inicjalizacja mumpsa 
	id.job=-1; id.par=1; id.sym=0; id.comm_fortran=-987654;	
	// id.sym = 2 dla symetrycznej, 1 dla symetrycznej dodatnio okreslonej
	dmumps_c(&id);
	
	//dodawanie jedynki, bo to fortran
	for(int i=0;i<R_N;i++)	{w_r[i]++; w_c[i]++;}

	// wprowadzenie danych 
	id.n = W;	    // rzad macierzy
	id.nz = R_N;	// ilosc niezerowych elementow
	id.irn = w_r;	// tablica indeksow wierszy
	id.jcn = w_c;	// tablica indeksow kolumn
	id.a = B;	    // tablica wartosci
	id.rhs = y;	    // prawa strona
	
	// analiza i rozwiazanie 
	id.icntl[1-1]=-1;	// zabijamy wszystkie wyjscia
	id.icntl[2-1]=-1;	// ...
	id.icntl[3-1]=-1;	// ...
	id.icntl[4-1]=0;	// poziom diagnostyki (max 4)
	id.job=6;		    // 6 - robie wsystko, jak w pralce
	dmumps_c(&id);
	
	// koniec zadania 
	id.job=-2;		    // -2 - konczy zadanie
	dmumps_c(&id);

	
	return id.info[0];
}

