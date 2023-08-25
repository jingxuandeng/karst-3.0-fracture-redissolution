/**
 * @file
 * @author  Agnieszka Budek
 * @date 25/09/2019
 *
 * @class Node
 *
 * @section DESCRIPTION
 *
 * Class Node represents the single node in a network of connected pores.
 * Is has attributes such as pressure and concentration.
 *
 */

#ifndef NODE_H
#define NODE_H 

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>

#include "pore.h"
#include "grain.h"
#include "printing.h"
#include "constants.h"

using namespace std;

class Pore;
class Grain;
class Point;
class Network;

/// represents one node
class Node{

	public:

		double u;		///< pressure in this node
		double cb;		///< concentration of acid (dissolving species) in this node
		double cc;		///< concentration of precipitating species in this node
		
		int b;			///< number of neighbors
		int bG;			///< number of attached grains

		Node  **n;		///< list of neighboring nodes
		Pore  **p;		///< list of neighboring pores
		Grain **g;		///< list of neighboring grains
		Point  xy;		///< node position (for visualization purpose mostly)
	
		int a;			///< pore number (name)
		double tmp;		///< temporal informations about the node;
		signed char t;	///< type of node: 1 -inlet node, 0 - normal node, -1 - outlet node
		int         x;  ///< additional info about pore (if connected to the dissolution pattern)

	public:

		Node  (int bb, float tt = 0);
		Node  (int name, int b_tmp, int t_tmp, Point point);
		~Node (){delete[] n; n=NULL; delete[] p; p=NULL; if(bG>0) delete[] g; g=NULL;}



		void change_node_neighbours(Node *n2,Node *n1, Pore *p);///< change neighbor n2 to n1
		void remove_form_IO_list (Network *S);					///< removing info about node from list of inlet and outlet nodes
		bool is_connected_to_Node(Node *n_tmp);					///< if true node n_tmp is a neighbor
		bool is_connected_to_Pore(Pore *p_tmp);					///< if true pore p_tmp is connected to the node
		int  neighbor_number     (Node *n_tmp);                 ///< returns nr of neighbour in the node list, if not connected, returnes -1
		void add_new_Grain   (Grain *g_tmp);					///< add g_tmp to the list of grains (without stretching the list)
		void add_Grain   	 (Grain *g_tmp);					///< add g_tmp to the list of grains
		void remove_Grain    (Grain *g_tmp);					///< remove g_tmp to the list of grains
		void add_neighbor    (Node *n_tmp, Pore *p_tmp);		///< add n_tmp and p_tmp to the list of neighbors
		void remove_neighbor (Node *n_tmp);						///< remove n_tmp and proper pore form the list of neighbors
		void remove_neighbor (Pore *p_tmp);						///< remove p_tmp and proper node from the list of neighbors

		double total_abs_flow();								///< returns total absolute flow through neighbor pores
		double total_pores_d ();								///< returns total diameter of neighbor pores

		//pattern analysis
		void check_diss_pattern(double treshold);   			///< set x=1 for pores that are connected to the dissolution pattern
		void find_forks();                                      ///< decides if the nod is a for and what is the fork generation
		int cluster_size (int l);                               ///<returns cluster size with the root in the node and length given by l ()
		void check_cluster(int l, int &sum, bool &if_l);        ///<cheching the cluster
};

ofstream_txt & operator << (ofstream_txt & stream, Node &p);
ostream      & operator << (ostream      & stream, Node &p);



inline void swap_nodes(Node **n1, Node **n2){
	Node *n_tmp;
	n_tmp=*n1;
	*n1=*n2;
	*n2=n_tmp;
	double tmptmp = (*n2)->tmp;
	(*n2)->tmp = (*n1)->tmp;
	(*n1)->tmp = tmptmp;
}

#endif


