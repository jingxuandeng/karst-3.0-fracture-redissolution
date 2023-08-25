#include "pore.h"
#include "network.h"
#include "printing.h"
#include "grain.h"

Pore::Pore (double dd, double ll, float name, int bb){
	d = dd; l = ll; a=name; tmp=name; q=0; x=1; bG=bb; c_in=0;

	n[0]=NULL; n[1]=NULL;	
	if(bG>0){
		g = new Grain*[bG];
		for (int i=0;i<bG;i++) g[i]=NULL;
	}
}



/**
* This function returns the inlet concentration of species B for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_inlet_cb(){
	if (q>0)  return n[0]->cb;
	else      return n[1]->cb;
}

/**
* This function returns the outlet concentration of species B for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_outlet_cb(){
	if (q>0) return n[1]->cb;
	else     return n[0]->cb;
}

/**
* This function returns the inlet concentration of species C for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_inlet_cc(){
	if (q>0)  return n[0]->cc;
	else      return n[1]->cc;
}

/**
* This function returns the outlet concentration of species C for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_outlet_cc(){
	if (q>0) return n[1]->cc;
	else     return n[0]->cc;
}


/**
* This function returns true if the species A is still present in neighboring grains.
* This information is important when deciding if the dissolution can still occur in the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Pore::is_Va_left(){
	double Va_tot=0;
	for (int b=0;b<bG;b++) Va_tot+=g[b]->Va;
	if(Va_tot>0)   return true;
	else		   return false;
}


/**
* This function returns the maximal length the pore can obtain taking into account geometric constrains.
*
* WARTNING: Periodic boundary conditions are traced correctly for "hexagonal" and "triangulation" type of topology.
* Other geometries need their own tracing of periodicity.
*
* ATTENTION: Maximal length can be exceed by the pore when merging in on
* -- the geometric constrains are less important to us than the local reaction area and resistivity of the system
* @param S pointer to the network
* @param l_max maximal length for the network, important for periodic boundary conditions
* @param l_0 characteristic pore's length, if there is a problem with calculation max length l_0 is returned
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::calculate_maximal_length(Network *S, double l_max, double l_0){

	cerr<<"WARNING: function calculate_maximal_lengthid is deprecated, use Network::point_distance(Piont*, Piont *) insted."<<endl;
	return;

	l = S->point_distance(n[0]->xy, n[1]->xy);

}

/**
* This function returns the actual length the pore if the if_dynamical_length flag is on.
*
* In this case the length of a pore changes in time depending on dissolution and precipitation processes.
* @param S pointer to the network
* @param l_max maximal length for the network, important for periodic boundary conditions
* @param l_0 characteristic pore's length, if there is a problem with calculation max length l_0 is returned
* @author Agnieszka Budek
* @date 26/09/2019
*/
void Pore::calculate_actual_length(Network *S, double l_max, double l_0){

	l = S->point_distance(n[0]->xy, n[1]->xy);

	if (d == 0) return;

	if(!S->if_dynamical_length || !S->if_track_grains) return;

	double factor=0; double V_max=0; double V_act=0;
	for (int s=0; s<bG; s++){
		V_max = g[s]->calculate_maximal_volume(S);
		V_act = g[s]->Va + g[s]->Ve + g[s]->Vx;
		if(V_max>0 && V_act>0) factor += pow(V_act/V_max,1./3);  //One may ask about exponent for semi 2D network, maybe it should be two?
		}

	if(factor>0) 	l = l*factor/bG;
	else			l = S->l_min;

	if(l<=S->l_min) l = S->l_min;

	if(l<=S->l_min) {
		if(S->if_verbose) cerr<<"l = l_min for Pore:"<<*this<<endl;
		for (int s=0; s<bG; s++){
			if(S->if_verbose) cerr<<"The following gains is going to be empty."<<endl\
					<<"Grain:"<<*g[s]<<endl;
			S->Va_tot -= g[s]->Va; g[s]->Va =0;
            S->Ve_tot -= g[s]->Ve; g[s]->Ve =0;


		}
	}
}

/**
* This function returns the local value of G parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_G(Network* S){

	if     (S->G1==0)	return	0;						         // reaction limited case, G = 0
	else if(S->G1>0 )	return  S->G1*d/S->d0;	             // mixed case: k1 ~ DD1
	else			    return  -1;							     // diffusion limited case, convention: G<0 => G = Inf

}

/**
* This function returns the local value of G2 (used in precipitation) parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_G_2(Network* S){

	if     (S->G2==0)	return	0;						         // reaction limited case, G = 0
	else if(S->G2>0 )	return  S->G2*d/S->d0;	             // mixed case: k1 ~ DD1
	else			    return  -1;							     // diffusion limited case, convention: G<0 => G = Inf

}

/**
* This function returns the local value of Da_eff parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_Da_eff(Network* S){

	if (q==0) return -1;
	double G = this->local_G(S);

	if      (G>0)    return S->Da*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G1)/(1+G));
	else if (G==0)   return S->Da*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	else             return S->Da*(l/S->l0)*(S->q_in_0/fabs(q));
}

/**
* This function returns the local value of Da_eff_2 (used in precipitation) parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_Da_eff_2(Network* S){

	if (q==0) return -1;
	double G = this->local_G_2(S);

	if      (G>0)    return S->Da2*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G2)/(1+G));
	else if (G==0)   return S->Da2*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	else             return S->Da2*(l/S->l0)*(S->q_in_0/fabs(q));
}



/**
* This function returns the change in diameter due to dissolution in one time step.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 14/03/2020
*/
double Pore::default_dd_plus(Network*S){

	if(S->if_track_grains && !is_Va_left())  return 0;   //no reaction if there is no A species available
	if(d==0 || q ==0)  return 0;   //pore with no flow
	if(l<=S->l_min)    return 0;   //no reaction in tiny grain

	//dissolution parameters
	double f1      = local_Da_eff(S);
	double g       = local_G(S);
	double c0;
	if(S->if_streamtube_mixing) c0 = c_in;
	else                        c0 = calculate_inlet_cb();

	double dd_plus = 0; 		//diameter change


	//finding dissolution contribution
	if      (f1==0)      dd_plus = 0;
	else if (S->G1 >=0)  dd_plus = S->dt*c0*(1-exp(-f1))/(1+g)/f1;
	else        	     dd_plus = S->dt*c0*(1-exp(-f1))/f1/d;


	return dd_plus;
}


/**
* This function returns the change in diameter due to precipitation in one time step (no condition for left space is checked here).
* @param S pointer to the network
* @author Agnieszka Budek
* @date 14/03/2020
*/
double Pore::default_dd_minus(Network*S){

	if(d==0 || q ==0)  return 0;   //pore with no flow
	if(l==S->l_min)    return 0;   //no reaction in tiny grain
	if(d<=S->d_min && (!is_Va_left())) return 0;

	//dissolution parameters
	double f1      = local_Da_eff(S);
	double g       = local_G(S);
	double c0;
	if(S->if_streamtube_mixing) c0 = c_in;
	else                        c0 = calculate_inlet_cb();


	//precipitation parameters
	double f2       = local_Da_eff_2(S);
	double c0_c     = calculate_inlet_cc();
	double dd_minus = 0; 		//diameter change


	//finding precipitation contribution
	if      (f2==0)         dd_minus = 0;
	else if (f1==f2 && is_Va_left())        dd_minus = S->gamma*S->dt/(1+g)/f1*((c0_c + c0)*(1-exp(-f1)) -c0*exp(-f1)*f1);
	else if (!is_Va_left()) 				dd_minus = S->gamma*S->dt/(1+g)/f1*  c0_c*      (1-exp(-f2));
	else                    				dd_minus = S->gamma*S->dt/(1+g)/f1*(\
								       	   	   	   c0  * (f1*(1-exp(-f2)) - f2*(1-exp(-f1)))/(f1-f2)+\
												   c0_c*  (1-exp(-f2)) );

	return dd_minus;
}

/**
* This function removes the pore form the list of neighboring nodes.
* This function is used when deleting the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::remove_info_from_attached_nodes(){

	for(int j=0;j<2;j++) n[j]->remove_neighbor(this);
}

/**
* This function removes the pore form the list of neighboring grains.
* This function is used when deleting the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::remove_info_from_attached_grains(){

	for(int j=0;j<bG;j++) g[j]->remove_Pore(this);
}

/**
* This function adds the g_tmp grain to the pore's list of neighboring grains.
* Important for merging and dynamic topology.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::add_grain (Grain *g_tmp){
	//if grain already in the list do nothing
	for(int i=0; i<bG; i++) if(g[i]==g_tmp) return;

	Grain ** g_new = new Grain* [bG+1];
	for(int i=0; i<bG; i++) g_new[i]=g[i];
	g_new[bG] = g_tmp;
	delete [] g; g = g_new;
	bG++;
}

/**
* This function removes the grain from the pore's list of neighboring grains.
* Important for merging and dynamic topology.
* @param gg grain to be removed
* @author Agnieszka Budek
* @date 25/09/2019
*/
void  Pore::remove_grain(Grain *gg){
	Grain ** g_new = new Grain * [bG];
	int b_new = 0;
	for(int bb=0;bb<bG;bb++) if(g[bb]!=gg) 	g_new[b_new++] = g[bb];

	delete [] g; g = g_new;
	bG = b_new;
}

/**
* This function returns true if a node belongs to the list of neighboring nodes.
* @param n_tmp node to be checked
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Pore::is_contain_Node  (Node *n_tmp){
	if(n[0]==n_tmp || n[1]==n_tmp) 	return true;
	else                			return false;
}

/**
* This function returns true if a grain belongs to the list of neighboring grains.
* @param g_tmp grain to be checked
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Pore::is_contain_Grain (Grain *g_tmp){
	for (int bb=0;bb<bG;bb++) if(g[bb]==g_tmp) return true;
	return false;
}

/**
* This function returns the common node of this pore and pore p1.
* @param p1 second pore
* @author Agnieszka Budek
* @date 25/09/2019
*/
Node*  Pore::find_common_node(Pore *p1){
	for(int i=0;i<2;i++) for(int j=0;j<2;j++) if(n[i]==p1->n[j]) return n[i];
	cerr<<"WARNING: Pores"<<*p1<<" and "<<*this<<" don't have common node."<<endl;
	return NULL;
}

/**
* This function returns the closest of two neighboring nodes.
* @param n_tmp second node
* @author Agnieszka Budek
* @date 25/09/2019
*/
Node*  Pore::find_closest_node(Node *n_tmp){

	double u_0 = fabs(n[0]->u - n_tmp->u);    //distance to node 0
	double u_1 = fabs(n[1]->u - n_tmp->u);    //distance to node 1
	if(u_1>u_0) return  n[0];
	else        return  n[1];

}

/**
* This function switches n_old to n_new in the list of neighboring nodes.
* @param n_old old node
* @param n_new new node
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::change_pore_neighbours(Node * n_old, Node *n_new){
	//changing info about nodes in grains
	for (int i=0;i<bG;i++) g[i]->change_nodes(n_old,n_new);
	//changing info about nodes in pores
	for (int j=0; j<2;j++) if(n[j]==n_old) {
		n[j] = n_new;
		if(n[(j+1)%2]==n_new) {
			cerr<<"WARNING: In changing pore neighbors generating pore with two identical nodes."<<endl<<*this<<endl;
			//d=0;
		}
		return;}
	cerr<<"WARNING: Problem with changing neighbors in pore "<<*this<<endl;
}




ofstream_txt & operator << (ofstream_txt & stream, Pore &p){

	stream <<setw(12)<<p.a<<setw(12)<<p.d<<setw(12)<<p.l<<setw(8)<<p.bG<<setw(12)<<p.q;
	return stream;
}



ostream & operator << (ostream & stream, Pore &p){

	stream<<"Pore("<<p.a<<"):  tmp = "<<p.tmp<<"  n = ("<<p.n[0]->a<<","<<p.n[1]->a<<")  g = (";
	for(int b=0;b<p.bG;b++){
		stream<<p.g[b]->a;
		if(b<p.bG-1) stream<<",";
	}
	stream<<")";
	return stream;

}
