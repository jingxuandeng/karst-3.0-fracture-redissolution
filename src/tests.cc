#include "network.h"
#include "printing.h"

int test_triangulation(){

	//int N_test = 10;
	while(true){

		//cerr <<endl<< "Testing triangulation: " << i<<endl<<endl;
		Network *S = new Network("config.txt");	//reading initial parameters form the file
		S->evolution(1);
		delete S;
	}
	return 0;
}
