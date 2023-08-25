#include "network.h"
#include "printing.h"

/**
* This function creates initial cut in the network
* - the region connected to the inlet of the system with different permeability that the rest of the system.
*
* @param cut_w width of the cut
* @param cut_l length of the cut
* @param factor factor for recycling d in the cut (by default eq 1)
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::create_an_inlet_cut(int cut_w, int cut_l, double factor){

	cerr<<"Adding an inlet cut..."<<endl;

	if(cut_l >=N_y) {cerr<<"WARNING : Cut length is to large. No cut is made."<<endl; return;}
	if(type_of_topology == "hexagonal"){
		for(int i=0; i<cut_l*N_x;i++){
			if(i%N_x >= (N_x/2. - cut_w/2.)  && i%N_x < (N_x/2. + cut_w/2.))
				for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0) n[i]->p[j]->d = d0*factor;}
	}

	else{
		for(int i=0;i<NN;i++){
			if(n[i]->xy.y <cut_l && n[i]->xy.x >= (N_x/2. - cut_w/2.) && n[i]->xy.x< (N_x/2. + cut_w/2.))
				for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0) n[i]->p[j]->d = d0*factor;}
	}
}


/**
* This function creates an initial pattern.
* It is useful when debugging.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::create_an_initial_pattern(){

//	cerr<<"Adding an initial pattern..."<<endl;
//
//	for (int i=0; i<NN; i++) if( i==7 || i == 13) p[2*i+0]->d=3*d0;
//	p[2*8+1]->d=3*d0;
//	p[2*7+1]->d=0.1*d0;

}
