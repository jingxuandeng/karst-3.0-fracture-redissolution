/**
 * @file
 * @author  Agnieszka Budek
 * @date 25/09/2019
 *
 * @section DESCRIPTION
 *
 * In this file you will find functions connected to the merging.
 *
 */


#include "network.h"
#include "printing.h"



/**
* This is a main function doing merging.
* It checks what type of merging to execute.
* If type_of_merging is "none", then no merging is done.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::do_merging(){

	//merging
	if      (type_of_merging == "none") 					cerr<<"No merging."<<endl;
	else if (type_of_merging == "merge_empty_grains") 		merge_empty_grains();
	else                               						cerr<<"WARNING: type_of_merging: "<<type_of_merging<<" not implemented yet."<<endl;

	//checking network connections
	if(type_of_merging != "none")     check_network_connections();


}

/**
* This function merge all empty grains.
* By empty grain we mean those with Va+Ve<=0. This condition can be change  by editing Network::to_be_merge().
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_empty_grains(){

	for(int i=0;i<NN;i++) n[i]->tmp = i;
	for(int i=0;i<NP;i++) p[i]->tmp = i;
	for(int i=0;i<NG;i++) g[i]->tmp = i;

	cerr<<"Merging empty grains..."<<endl;

	if(if_debugging_printing && if_verbose) debugging_for_merging ("Before merging one grain: s = " + to_string(tot_steps));

	for(int i=0;i<NG;i++) if(g[i]->to_be_merge()) { //condition for merging, can be adapted later
		cerr<<"I am merging this guy: "<<*g[i]<<endl;

		if(g[i]->bP<=1) merge_one_pore_grain(g[i]);
		else            merge_one_grain(g[i]);

		//additional checking of network structures
		if(if_debugging_printing && if_verbose) debugging_for_merging ("After merging one grain: s = " + to_string(tot_steps));
	}

	//Clear pathological pores: (probably unnecessary)
	for(int i=0; i<NP;i++) if(p[i]->n[0]==p[i]->n[1]) {
		Pore *p_tmp = p[i];
		cerr<<"Clearing pathological pores: "<<*p_tmp<<endl;
		p_tmp->remove_info_from_attached_nodes();
		p_tmp->remove_info_from_attached_grains();
		move_pore_to_the_end(p_tmp->tmp,NP);
	}

	if(if_debugging_printing && if_verbose) debugging_for_merging ("After clearing pathological pores: s = " + to_string(tot_steps));

	//Clear pathological nodes: (probably unnecessary)
	for(int i=0; i<NN;i++) if(n[i]->b==0) {
		Node *n_tmp = n[i];
		cerr<<"Clearing pathological nodes: "<<*n_tmp<<endl;
		n_tmp->remove_form_IO_list(this);
		move_node_to_the_end(i,NN);
	}

	if(if_debugging_printing && if_verbose) debugging_for_merging ("After clearing pathological nodes: s = " + to_string(tot_steps));


	//Merging pathological grains
	for(int i=0;i<NG;i++) if(g[i]->is_pathological()){
		if(g[i]->bP<=1) merge_one_pore_grain(g[i]);
		else            merge_one_grain(g[i]);
		if(if_debugging_printing && if_verbose) debugging_for_merging ("After merging pathological grain: s = " + to_string(tot_steps));
	}

}

/**
* This function merge all pathological grains connected only to aone pore.
* Those grains can result from merging normal grain.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_one_pore_grain (Grain *gg){


	if(if_verbose) cerr<<endl<<endl<<"Merging one-pore-grain: "<<*gg<<endl<<"g_tmp = "<<gg->tmp<<endl;
	for(int i=0;i<gg->bN;i++)	gg->n[i]->remove_Grain(gg);
	for(int i=0;i<gg->bP;i++)   gg->p[i]->remove_grain(gg);

	//updating the total volume of A species
	if (gg->bP>0){
		if(gg->p[0]->bG>0)  for(int i=0; i<p[0]->bG;i++) p[0]->g[i]->Va+=gg->Va/p[0]->bG;
		else 				Va_tot -= gg->Va;}
	else          			Va_tot -= gg->Va;
	move_grain_to_the_end(gg->tmp,NG);

}

/**
* This function merge one grain.
* In a result this grain will vanish and all connected pores will projected to one of them, the same with nodes.
*
* @param gg grain to be merged
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_one_grain(Grain *gg){


	if(if_verbose)cerr<<endl<<endl<<"Merging grain: "<<*gg<<endl;
	//choose main pore
	Pore *p1=NULL;  //master pore (the one that will remain)
	double q_max=0;
	for(int b=0; b<gg->bP;b++) if(fabs(gg->p[b]->q)>=q_max){
		p1 = gg->p[b];
		q_max = fabs(gg->p[b]->q);
	}
	if(if_verbose)cerr<<"Master pore is "<<*p1<<endl;
	gg->set_effective_d_and_l(p1,this);  //calculate d and l basing on S and R of the entire pore


	//temporary list of grains
	int b_tmp=0;
	for(int b=0;b<gg->bP;b++){b_tmp+=gg->p[b]->bG;}
	Grain **g_tmp = new Grain *[b_tmp];
	b_tmp=0;
	for(int b=0;b<p1->bG;b++) if(p1->g[b]!=gg){ g_tmp[b_tmp++]=p1->g[b]; } //copy list of pores and nodes

	//clearing unnecessary pores
	Grain gg_copy = *gg;
	for(int b=0;b<gg_copy.bP;b++) if(gg_copy.p[b]!=p1) {
		Pore * p_tmp=gg_copy.p[b];
		//checking if new grain will be attached to p1 (if p_tmp nodes are not projected to the same p1 node)
		Node *n1_tmp  = p1->find_closest_node(p_tmp->n[0]);
		Node *n2_tmp  = p1->find_closest_node(p_tmp->n[1]);
		//adding new grains to the list of p1 grains
		if(n1_tmp!=n2_tmp) for(int bb=0;bb<p_tmp->bG;bb++) {
			bool if_new = true;
			if(p_tmp->g[bb]==gg)								if_new = false;
			for(int j=0;j<b_tmp;j++) if(p_tmp->g[bb]==g_tmp[j]) if_new = false;
			if(if_new) { g_tmp[b_tmp++] = p_tmp->g[bb]; p_tmp->g[bb]->change_pores(p_tmp,p1);}
		}
		//adding info about flow to master pore (for better printing only)
		p1->q = p1->q + _sign(p1->q)*abs(p_tmp->q)/3;

		//clearing unused pore
		if(if_verbose) cerr<<"Clearing pore: "<<*p_tmp<<endl;
		p_tmp->remove_info_from_attached_nodes ();
		p_tmp->remove_info_from_attached_grains();
		move_pore_to_the_end(p_tmp->tmp,NP);
	}

	p1->bG = b_tmp;
	delete [] p1->g;  p1->g = g_tmp;
	//cerr<<"Glupi element po usuwaniu porow: "<<*glupi_element<<endl;

	//merging all nodes but two (p1->n[0] and p1->n[1])
	gg_copy = *gg;
	for(int b=0;b<gg_copy.bN;b++) if(gg_copy.n[b]!=p1->n[0] && gg_copy.n[b]!=p1->n[1]) {
		Node * n_master = p1->find_closest_node(gg_copy.n[b]);    //closer p1 node
		merge_nodes(n_master,gg_copy.n[b]);
	}
	//cerr<<"Glupi element po laczeniu nodow: "<<*glupi_element<<endl<<flush;
	//changing info about grains in master nodes
	p1->n[0]->remove_Grain(gg);
	p1->n[1]->remove_Grain(gg);
	//cerr<<"Glupie element niedaleko konca: "<<*glupi_element<<endl<<flush;


	if (p1->bG>0) for(int i=0;i<p1->bG;i++) p1->g[i]->Va+=gg->Va/p1->bG;
	else          Va_tot -= gg->Va;
	move_grain_to_the_end(gg->tmp,NG);

}


/**
* This function do special debugging for merging
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::debugging_for_merging(string text){

	check_network_connections();

	for(int i=0;i<NN;i++) n[i]->tmp = n[i]->a;
	for(int i=0;i<NP;i++) p[i]->tmp = p[i]->a;
	for(int i=0;i<NG;i++) g[i]->tmp = g[i]->Va;
	description_note = text;
	net_ps<<*this;
	for(int i=0;i<NN;i++) n[i]->tmp = i;
	for(int i=0;i<NP;i++) p[i]->tmp = i;
	for(int i=0;i<NG;i++) g[i]->tmp = i;

}

/**
* This function merge two nodes connected to the grain that will be deleted..
*
* @param n1 master node
* @param n2 node to be deleted
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_nodes(Node *n1, Node *n2){

	if(if_verbose) cerr<<"Merging nodes: "<<endl<<*n1<<endl<<*n2<<endl;
	if(n1->a<0 || n2->a<0) cerr<<"ERROR in merging nodes!!!!"<<endl;


	//Updating list of neighbors
	Node  ** n_tmp = new Node  *[n1->b +n2->b];
	Pore  ** p_tmp = new Pore  *[n1->b +n2->b];
	Grain ** g_tmp = new Grain *[n1->bG+n2->bG];
	int b_tmp=0; int bg_tmp=0;
    for(int b=0;b<n1->b;b++)  if(n1->n[b]!=n2){ n_tmp[b_tmp]=n1->n[b];  p_tmp[b_tmp++]=n1->p[b];} //copy list of pores and nodes
	Node n2_copy = *n2;
    for(int b=0;b<n2_copy.b;b++)
		if(n2_copy.n[b]!=n1){
			n_tmp[b_tmp]   = n2_copy.n[b];
			p_tmp[b_tmp++] = n2_copy.p[b];
			n2_copy.n[b] -> change_node_neighbours(n2,n1,n2_copy.p[b]);
			n2_copy.p[b] -> change_pore_neighbours(n2,n1);
			}
		else{ // To not leave unattached pores
			Pore *p_del = n2_copy.p[b];
			if(if_verbose) cerr<<"Clearing pore: "<<*p_del<<endl;
			//p_del->remove_info_from_attached_nodes();  //
			p_del->remove_info_from_attached_grains();
			move_pore_to_the_end(p_del->tmp,NP);
		}

	for(int b=0;b<n1->bG;b++) g_tmp[bg_tmp++]=n1->g[b];   //copy list of grains
	for(int b=0;b<n2->bG;b++) {
		bool if_new=true;
		for(int j=0;j<bg_tmp;j++) if(g_tmp[j]==n2->g[b]) if_new = false;
		if(if_new) g_tmp[bg_tmp++]=n2->g[b];
		n2->g[b]->change_nodes(n2,n1);
	}

	//update neighbors for n1
	n1->b = b_tmp; n1->bG = bg_tmp;
	delete [] n1->n; n1->n = n_tmp;
	delete [] n1->p; n1->p = p_tmp;
	delete [] n1->g; n1->g = g_tmp;

	//change position of node n1 (if the distance is not to large - see periodic boundary conditions)
	if(printing_mode!="grains") if(n1->xy - n2->xy<N_x/2)  {
		double q1 = n1->total_pores_d();
		double q2 = n2->total_pores_d();
		if(q1+q2 == 0) {q1 = 1; q2 = 1;}
		n1->xy = ((2*q1/(q1+q2)) * n1->xy) * ((2*q2/(q1+q2))*n2->xy);
	}

	//change type of master node if necessary
	if(n2->t!=0){
		if(n1->t==0){
			n1->t = n2->t;
			if(n2->t==1)  for(int i=0; i<N_wi; i++) if(wi[i]==n2) {wi[i] = n1; break;}
			if(n2->t==-1) for(int i=0; i<N_wo; i++) if(wo[i]==n2) {wo[i] = n1; break;}
		}
		else  n2->remove_form_IO_list(this); //remove n2 from wi/wo list
	}


	//clear node n2
	move_node_to_the_end(n2->tmp,NN);

	if(if_verbose) cerr<<"After merging nodes: "<<endl<<*n1<<endl<<*n2<<endl<<flush;
}


