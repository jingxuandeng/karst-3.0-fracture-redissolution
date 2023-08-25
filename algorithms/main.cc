#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>


#include "algorithms_cc.h"
#include "../include/triangulation.h"

//#define VISUALIZATION //only for debugging

#ifdef VISUALIZATION
#include <SFML/Graphics.hpp>
#endif

using namespace std;



/**
* This function executes the main test for this project.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void test_matrix_solving(){

	int* w_k  = new int    [100];
   	int* w_i  = new int    [100];
   	double* B = new double [100];
	double* y = new double [100];

	for(int i=0;i<10;i++) for(int j=0;j<10;j++) {w_k[i*10+j]=i; w_i[i*10+j]=j; B[i*10+j]=12;  y[i*10+j]=666;}


	solve_matrix(10, 100,  w_k,  w_i, B,  y);

	//for(int i=0;i<100;i++) cerr<<"y[i] = "<<y[i]<<endl;
	for(int i=0;i<100;i++) cerr<<"y[i] = "<<y[i]<<" B[i] = "<<B[i]<<endl; 

	return;
}


/**
* This function executes some additional test for this project.
* It generates the triangulation.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void test_triangulation(){
	int NN,N_p,N_g;
	std::vector<Vector2 <double> >  points;
	std::vector<Vector2 <double> >  points_tmp;
	std::vector<Triangle<double> >  triangles;
	std::vector<Edge    <double> >  edges;

	triangulation(5, 5, points, points_tmp, edges, triangles, false, true);
}

/**
* In main function by now only tests are executed.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
int main(){
	//test_matrix_solving();
	test_triangulation();
}
