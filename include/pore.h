/**
 * @file
 * @author  Agnieszka Budek
 * @date 25/09/2019
 *
 * @class Pore
 *
 * @section DESCRIPTION
 *
 * Class Pore represents the single pore in a porous material.
 * Is has attributes such as diameter and length (the pore is represented here as a tube).
 *
 */

#ifndef PORE_H
#define PORE_H 

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>

#include "constants.h"

#include "node.h"
#include "grain.h"
#include "network.h"

class Node;
class Network;
class Grain;
class Point;

using namespace std;


/// represents one pore
class Pore{

	public:

		double d;		///< pore diameter
		double l;		///< pore length
		double q;		///< flow through the pore
		double c_in;	///< concentration at the pore inlet
		int    a;       ///< pore number (name)

		int bG;		    ///< number of grains in vicinity


		Node  *n[2];	///< node list
		Grain **g;	    ///< grain list

		signed char x;  ///< additional info about pore (if connected to the dissolution pattern)
		double tmp;		///< temporal informations (info to be printed)


	public:

		
		Pore (double dd = 0.02, double ll = 1, float name=0, int bb=2);
		~Pore ()						{if(bG>0) delete[] g; g=NULL; }

		double perm(double mu_0){return M_PI*pow(d,4)/(128*mu_0*l);}	///< permeability of a particular pore
		//void   diss (double Va, double Ve);							///< precipitation and dissolution of the material: calculate change of d and l
		double calculate_inlet_cb();									///< calculate inlet concentration of the species B
		double calculate_outlet_cb();									///< calculate outlet concentration of the species B
		double calculate_inlet_cc();									///< calculate inlet concentration of the species C
		double calculate_outlet_cc();									///< calculate outlet concentration of the species C

		void   calculate_maximal_length(Network *S = NULL, double l_max=10, double l_0=1);	///< calculate maximal pore length
		void   calculate_actual_length (Network *S = NULL, double l_max=10, double l_0=1);	///< calculate initial pore length
		double local_G        (Network* S);			///< dissolution parameters
		double local_Da_eff   (Network* S); 		///< dissolution parameter
		double local_G_2      (Network* S);      	///< precipitation parameter
		double local_Da_eff_2 (Network* S);      	///< precipitation parameter
		bool   is_Va_left();						///< return false if there is no Va material left
		double default_dd_plus(Network*S);		///< change in diameter as a result of dissolution
		double default_dd_minus(Network*S);    ///< default change in diameter as a result of precipitation (no space condition is checked)


		void   remove_info_from_attached_nodes(); 			///< remove this pore from the list of connected nodes
		void   remove_info_from_attached_grains();			///< remove this pore from the list of connected grains
		void   remove_grain(Grain *g);						///< remove g from the list of this pore
		void   add_grain   (Grain *g_tmp);					///< add g_tmp to the list of this pore
		void   change_pore_neighbours(Node * n2, Node *n1); ///< change pore n1 to n2 for this pore
		bool   is_contain_Node  (Node *n_tmp);				///< check if node n_tmp belongs to the neighbors of this pore
		bool   is_contain_Grain (Grain *g_tmp);				///< check if node n_tmp belongs to the neighbors of this pore
		Node*  find_common_node (Pore *p1);					///< check if node n_tmp belongs to the neighbors of this pore
		Node*  find_closest_node(Node *n1);					///< check if node n_tmp belongs to the neighbors of this pore


};

ofstream_txt & operator << (ofstream_txt & stream, Pore &p); ///< print info about pore p in one line
ostream      & operator << (ostream      & stream, Pore &p); ///< print more elaborate info about pore p

/**
* This function swaps two pores in a main list of pores
*/
inline void swap_pores(Pore **p1, Pore **p2){


	Pore *p_tmp = *p1;
	*p1=*p2;
	*p2=p_tmp;
	double tmptmp = (*p2)->tmp;
	(*p2)->tmp = (*p1)->tmp;
	(*p1)->tmp = tmptmp;

}

#endif


