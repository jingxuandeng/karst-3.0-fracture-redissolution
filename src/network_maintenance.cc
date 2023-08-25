#include "grain.h"
#include "network.h"
#include "printing.h"


/**
* This function returns the local Da_eff in a given pore.
* Currently not used.
*
* @param p a given pore
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::Da_eff_in_pore(Pore* p){   //currently not used

	if (p->q==0) return -1;
	double k_eff = k_eff_in_pore(p);
	return M_PI*p->d*k_eff*p->l/fabs(p->q);
}


/**
* This function calculates reaction rate for dissolution for given D_eff and G.
* Currently not used.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::recalculate_k1(){
	if (G1 < 0) k1  = -1; 		// diffusion limited case
	else        k1  = (Da*q_in_0 + Da*G1*q_in_0)/(d0*l0*M_PI);
	cerr<<"Reaction rate has been calculated: k1 = "<<k1<<endl;
	if(G1>=0) {dt_unit = 2*k1*gamma_1/d0; cerr<<"Unit of time is set to: [2*k1 * gamma_1/d0] = "<<dt_unit<<endl;}
}



/**
* This function calculates diffusion coefficient for dissolution for given D_eff and G.
* Currently not used.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::recalculate_DD1(){
	if (G1 == 0)     DD1 = -1;            								// diffusion limited case
	else if (G1<0)   DD1 = (Da*q_in_0)/(M_PI*Sh*l0);      				// reaction limited case; convention: G<0 => G = Inf
	else             DD1 = (Da*q_in_0 + Da*G1*q_in_0)/(G1*l0*M_PI*Sh);	// mixed case
	cerr<<"Transversal diffusion has been calculated: DD1 = "<<DD1<<endl;
	if(G1<0) {dt_unit = 2*DD1*Sh*gamma_1/d0; cerr<<"Unit of time is set to: [2*DD1*Sh*gamma_1/d0^2] = "<<dt_unit<<endl;}

}


/**
* This function returns the common pore for two nodes.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
Pore* Network::findPore(Node* n1, Node* n2){  //memory optimization, for runtime optimization should be done differently
	for(int i=0;i<n1->b;i++) for(int j=0;j<n2->b;j++) if(n1->p[i]==n2->p[j] && n2->p[j]!= NULL) if(n2->p[j]->n[0]!= NULL)
		return n1->p[i];
	return NULL;
}


/**
* This function returns the common pore for two nodes.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
Pore* findPore(Node* n1, Node* n2){                //memory optimization, for runtime optimization should be done differently
	for(int i=0;i<n1->b;i++) for(int j=0;j<n2->b;j++) if(n1->p[i]==n2->p[j] && n2->p[j]!= NULL) if(n2->p[j]->n[0]!= NULL)
		return n1->p[i];
	return NULL;
}




/**
* This function returns the common grain for tree nodes.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
Grain* Network::findGrain_T(Node* n1, Node* n2, Node *n3){     //memory optimization, for runtime optimization should be done differently
	for(int i=0;i<NG;i++)
		for(int b=0;b<3;b++)
			if(g[i]->n[b]==n1)
				for(int bb=0;bb<3;bb++)
					if(g[i]->n[bb]==n2)
						for(int bbb=0;bbb<3;bbb++)
							if(g[i]->n[bbb]==n3)
								if (g[i]!=NULL) if(g[i]->n[0]!=NULL)
								return g[i];
	return NULL;
}

/**
* This function returns the common grain for tree nodes.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
Grain* Network::findGrain(Node* n1, Node* n2, Node *n3){     //memory optimization, for runtime optimization should be done differently
	for(int b=0;b<n1->b;b++) for(int bb=0;bb<n1->p[b]->bG;bb++){
		 Grain * gg = n1->p[b]->g[bb];
		 for(int b=0;b<3;b++)
		 			if(gg->n[b]==n1)
		 				for(int bb=0;bb<3;bb++)
		 					if(gg->n[bb]==n2)
		 						for(int bbb=0;bbb<3;bbb++)
		 							if(gg->n[bbb]==n3)
		 								return gg;
	}
	return NULL;
}

/**
* This function returns current efficient reaction rate in a pore.
* Currently not used.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::k_eff_in_pore(Pore* p){ //currently not used

	if(G1==0)		return	k1;				                 // reaction limited case, DD1<<k1
	else if(G1>0)	return  (k1*Sh*DD1)/(Sh*DD1+k1*p->d);	 // mixed case: k1 ~ DD1
	else			return  DD1*Sh/p->d;				     // diffusion limited case, convention: G<0 => G = Inf

}

/**
* This function returns current G parameter in a pore.
* Currently not used.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::G_in_pore(Pore* p){     //currently not used

	if(G1==0)		return	0;						         // reaction limited case, G = 0
	else if(G1>0)	return  (k1*p->d)/(Sh*DD1);	             // mixed case: k1 ~ DD1
	else			return  -1;							     // diffusion limited case, convention: G<0 => G = Inf

}


/**
* This function checks if the time-step length should be changed.
* @param dd diameter change in a single pore in one time step
* @param dV volume change in a single pore in one time step
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::set_adaptive_dt(double dd,double dV){

	if (set_new_dt == 1) return;
	if (d_d_max > 0 && dd > d_d_max)   {set_new_dt = 1; return;}
	if (d_V_max > 0 && dV > d_V_max)   {set_new_dt = 1; return;}
	if (set_new_dt == -1 && d_d_min > 0 && dd > d_d_min) set_new_dt = 0;
	if (set_new_dt == -1 && d_V_min > 0 && dV > d_V_min) set_new_dt = 0;
}

/**
* This function updates the time-step length, dt.
* The dt is elongated if the dissolution occurs to slow or shorten if dissolution is to slow.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::adapt_dt(){
	if (set_new_dt == 0) {cerr<<"Time step set to, dt = "<<dt<<"."<<endl; return;}
	if (set_new_dt == 1)   dt = 0.75*dt;
	if (set_new_dt == -1)  dt = 1.2*dt;
	set_new_dt = -1;
	cerr<<"Time step has been adapted, dt = "<<dt<<"."<<endl;
}

/**
* This function checks if the system is fully dissolved.
* This happens when the dissolution pattern has obtained the end of the system.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::check_if_dissolved(){

	//condition for dissolution
	check_diss_pattern(d_d_dis);
	for(int i=0;i<N_wo;i++){
		Node* nn = wo[i];
		for (int j=0; j<nn->b;j++) if (nn->p[j]->d > d0*d_d_dis && nn->p[j]->x==1){
			cerr<<"\nSystem is fully dissolved\nDissolution finished."<<endl;
			if_system_dissolved = true;
			return;}
	}

	//condition for clogging
	int nr_of_red_grains=0;
	if(if_precipitation) for(int i=0;i<N_wo;i++){
		Node* nn = wo[i];
		for (int j=0; j<nn->bG;j++) if(nn->g[j]->Va == 0 && nn->g[j]->bP>2 && nn->g[j]->Ve>l0*l0/10) nr_of_red_grains++;

	}
	 if(nr_of_red_grains>double(N_wo)/20.){
		cerr<<"\nSystem is dissolved\nDissolution and precipitation is finished."<<endl;
		if_system_dissolved = true;
	}

	//condition for pressure drop
	 if(Q_tot>0 && u_min>0) if (wi[0]->u/N_y * Q_tot/(2*N_x)<u_min){
		 cerr<<"\nSystem is dissolved: condition for pressure has been fulfilled."<<endl;
		 if_system_dissolved = true;
	 }



}

//setting variable x=1 in pores that are connected to the dissolution pattern
/**
* This function checks position of a dissolution pattern defined by the threshold.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::check_diss_pattern (double threshold){
	for (int i =0;i<NN;i++)   n[i]->x=0;
	for (int i =0;i<NP;i++)   p[i]->x=0;
	for (int i =0;i<N_wi;i++) wi[i]->check_diss_pattern(threshold*d0);
	for (int i =0;i<N_wi;i++) wi[i]->x=2;
}

/**
* This function checks if the dissolution pattern has reached the n-th row
*
* @param row row to be checked
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Network::check_diss_front(double threshold, int row){
	check_diss_pattern(threshold);
	for(int i=0;i<NP;i++){
		int x = (p[i]->n[0]->xy.y + p[i]->n[1]->xy.y)/2;
		if(x==row && p[i]->x==1) {cerr<<"Printing row: "<<row<<endl;  return true;}
	}
	return false;
}



/**
* This function deletes pores that are not connected to the inlet of the system.
* It may be useful when reading network from the NMR scan.
* In the general case there may be pores not connected to the rest of the system.
* Unconnected pores and nodes may cause a problem while calculating a pressure filed.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::clear_unconeccted_pores(){

	check_diss_pattern(0);
	for(int i=0;i<NP;i++) if(p[i]->x == 0) p[i]->d=0;
	clear_unused_pores();
	clear_unused_nodes();
	clear_unused_grains();

}



/**
* This function deletes unused pores.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::clear_unused_pores(){
	int NP_tmp = NP;
	for(int i=0;i<NP_tmp;i++) if(p[i]->d==0 && p[i]->a>=0){
		cerr<<p[i]->a<<" ";
		p[i]->remove_info_from_attached_nodes();
		p[i]->remove_info_from_attached_grains();
		move_pore_to_the_end(i,NP_tmp);
        i=-1;
	}
	cerr<<endl;
}


/**
* This function moves unused pores to the end of the pores list.
* Those pores are treated as deleted.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::move_pore_to_the_end(int i, int NP_tmp){

	p[i]->a = -p[i]->a-1;
	int w = NP_tmp-1;
	while(p[w]->a<0 && w>0) w--;
	if(w>i) swap_pores(&p[i], &p[w]);
	NP--;

}

/**
* This function deletes unused nodes.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::clear_unused_nodes(){
	int NN_tmp = NN;
	for(int i=0;i<NN_tmp;i++) if(n[i]->b==0 && n[i]->a>=0){
		Node *nn=n[i];
		cerr<<n[i]->a<<" ";
		nn->remove_form_IO_list(this);
		for(int b=0;b<nn->bG;b++) nn->g[b]->remove_Node(nn);
		move_node_to_the_end(i,NN_tmp);
		i=-1;
	}
	for(int i=0;i<NN_tmp;i++) {
		if(n[i]->b==0) {
			n[i]->t=2;}
	}
		cerr<<endl;
}

/**
* This function moves unused nodes to the end of the nodes list.
* Those nodes are treated as deleted.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::move_node_to_the_end(int i, int NN_tmp){
	n[i]->a = -n[i]->a-1;
	int w = NN_tmp-1;
	//while(n[w]->b==0&& w>0) w--;
	while(n[w]->a<0&& w>0) w--;
	if(w>i) swap_nodes(&n[i], &n[w]);
	NN--;
}


/**
* This function deletes unused grains.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::clear_unused_grains(){

	int NG_tmp = NG;
	for(int i=0;i<NG_tmp;i++) if(g[i]->bP==0 && g[i]->a>=0){
		cerr<<g[i]->a<<" ";
		for(int b=0;b<g[i]->bN;b++) g[i]->n[b]->remove_Grain(g[i]);
		move_grain_to_the_end(i, NG_tmp);
		i=-1;
	}
	cerr<<endl;

}

/**
* This function moves unused grains to the end of the grains list.
* Those grains are treated as deleted.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::move_grain_to_the_end(int i, int NG_tmp){
	g[i]->a = -g[i]->a-1;
	int w = NG_tmp-1;
	while(g[w]->a<0 && w>0) w--;
	if(w>i) swap_grains(&g[i], &g[w]);
	NG--;

}


//return distance between nodes including boundary conditions
/**
* This function returns the distance between two nodes.
* Periodic boundary conditions are taking into account.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::node_distance(Node *n1, Node  *n2){
	return point_distance(n1->xy, n2->xy);
}


//return distance between points including boundary conditions
/**
* This function returns the distance between two points.
* Periodic boundary conditions are taking into account (but only in x directions and N_x must be set correctly).
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::point_distance(Point p1, Point  p2){

	double N_x_tmp = N_x;
	double N_y_tmp = N_y;
	if(type_of_topology=="hexagonal") N_y_tmp = N_y*sqrt(3)/2;

	double xx1= pow(p1.x - p2.x,2);
	double xx2;
	if (p1.x < p2.x ) 	xx2 = pow(p1.x - (p2.x - N_x_tmp),2);
	else  				xx2 = pow(p2.x - (p1.x - N_x_tmp),2);

	double yy1 = pow(p1.y - p2.y,2);
	double yy2;
	if (p1.y < p2.y ) 	yy2 = pow(p1.y - (p2.y - N_y_tmp),2);
	else  				yy2 = pow(p2.y - (p1.y - N_y_tmp),2);

	double xx;
	double yy;
	double zz = pow(p1.z - p2.z,2);

	if (xx1<xx2 || !if_periodic_bc) 	xx = xx1;
	else								xx = xx2;


	if (yy1<yy2 || !if_periodic_bc) 	yy = yy1;
	else								yy = yy2;

	return sqrt(xx+yy+zz);
}


/**
* This function saves information about initial nodes position.
* It is used only for printing in grain style mode.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::save_info_about_initial_node_positions (){ //for printing in grain style only

	cerr<< "Saving info about initial node positions..."<<endl;

	if(NG==0) return;

	if (initial_xy==NULL) {
		initial_xy = new Point* [NG];
		for(int i=0; i<NG; i++)initial_xy[i]=NULL;
	}

	for(int i=0; i<NG; i++){
		Grain* gg = g[i];
		if (initial_xy[i]==NULL)    initial_xy[i] = new Point[gg->bN];
		for(int j=0; j<gg->bN; j++) initial_xy[i][j] = gg->n[j]->xy;
	}

}
