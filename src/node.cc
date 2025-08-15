#include "node.h"
#include "network.h"

Node::Node (int bb, float name)	{
	b=bb; a=name; tmp=name; t=0; x=1;
    is_fracture=false; is_LHS=false;
	cb=0; cc=0; u=0;
	bG = 0; g=NULL;
	n = new Node*[b];
	p = new Pore*[b];
	for(int i=0;i<b;i++){ n[i]=NULL; p[i]=NULL;}
	//if(b!=6) cerr<<"Problem with hexagonal network!"<<endl;
	//	cerr<<"Creation of node nr: "<<int(tt)<<" b = "<<(int) b<<" tmp= "<<int(tmp)<<endl;
}


Node::Node  (int name, int b_tmp, int t_tmp, Point point){
	b=b_tmp; a=name; tmp=name; t=t_tmp; x=1;
    is_fracture=false; is_LHS=false;
	cb=0; cc=0; u=0;
	bG = 0; g=NULL;
	xy = point;
	n = new Node*[b];
	p = new Pore*[b];
	for(int i=0;i<b;i++){ n[i]=NULL; p[i]=NULL;}
}

/**
* This function returns the sum of absolute flow through all neighboring pores.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Node::total_abs_flow(){
	double Q_tmp = 0;
	for(int i=0; i<b; i++) Q_tmp+= abs(p[i]->q);
	return Q_tmp;
}


/**
* This function returns the diameters sum of all neighboring pores.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Node::total_pores_d(){
	double d_tmp = 0;
	for(int i=0; i<b; i++) d_tmp+= p[i]->d;
	return d_tmp;
}

double Node::total_pores_q(){
    double q_tmp = 0;
    for(int i=0; i<b; i++) q_tmp+= fabs(p[i]->q);
    return q_tmp;
}

/**
* This function removes node form input/output list of nodes.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::remove_form_IO_list(Network *S){

	if(t==1) for(int i=0; i<S->N_wi; i++) if(S->wi[i]==this) {
		S->wi[i] = S->wi[S->N_wi-1];
		S->wi[S->N_wi-1] = NULL;
		S->N_wi--;
		break;
	}

	if(t==-1) for(int i=0; i<S->N_wo; i++) if(S->wo[i]==this) {
			S->wo[i] = S->wo[S->N_wo-1];
			S->wo[S->N_wo-1] = NULL;
			S->N_wo--;
			break;
		}
}

/**
* This function checks if the node belongs to the dissolution pattern.
* The pattern is defined as a compact collection of all pores for which d>d_0 * threshold which is connected to the inlet of the system.
* @param threshold threshold of diameter defining the dissolution pattern.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::check_diss_pattern(double threshold){
	x = 1;  			//the node has been checked
	for(int i=0;i<b;i++) if(p[i]->d>threshold || p[i]->d<=-threshold){
		p[i]->x = 1;  	//pore has been set as connected to the dissolution pattern
		x = 2;  		//node is connected to the pattern
		if (n[i]->x==0) n[i]->check_diss_pattern(threshold);
	}
}
void Node::check_preci_pattern(double factor){
    x = 1;  			//the node has been checked
    double max_distance = 10;
    for(int i=0;i<bG;i++)
        if(g[i]->bN==3)             //TODO : Make it work for general bG, not only 3
            if(fabs(g[i]->n[0]->xy.y-g[i]->n[1]->xy.y<max_distance) &&
               fabs(g[i]->n[1]->xy.y-g[i]->n[2]->xy.y<max_distance) &&
               fabs(g[i]->n[2]->xy.y-g[i]->n[0]->xy.y<max_distance))
                if(g[i]->Ve > g[i]->Va*factor) { // && g[i]->Ve>0.000001) {
                    bool if_fracture=false;
                    for(int j=0;j<g[i]->bP;j++) if(g[i]->p[j]->is_fracture) if_fracture=true;
                    if(if_fracture) continue;
                    g[i]->x = 1;  	//grain has been set as connected to the precipitation pattern
					for (int s=0;s<g[i]->bP; s++) g[i]->p[s]->x = 1;			//pore has been a part of a precipitation pattern
                    x = 2;  		//node is connected to the pattern
                    for (int b=0; b<g[i]->bN;b++)
                        if (g[i]->n[b]->x==0)
                            g[i]->n[b]->check_preci_pattern(factor);
                }
}
//Old version numbering fork form the root


void Node::find_forks(){

	int nr_of_children = 0;
	double old_tmp     = 0;

	//
	for(int s=0;s<b;s++){
		if(p[s]->x !=2) continue;
		if(n[s]->x == x+1)                     nr_of_children++;
		else if(n[s]->x == x-1 && n[s]->x>0)   old_tmp = fabs(n[s]->tmp);}

	//update tmp
	if(nr_of_children>1) tmp = old_tmp+1;
	else                 tmp = -old_tmp;

	for(int s=0;s<b;s++) if(n[s]->x == x+1 &&  p[s]->x==2)  n[s]->find_forks();
}

/**
* This function return the cluster size that origins in the node.
* @param l length of the cluster
* if cluster is shorter than l 0 is returned
* @author Agnieszka Budek
* @date 06/03/2020
*/
int Node::cluster_size(int l){

	int sum_tmp = 1;
	bool if_l = false;
	for(int s=0;s<b;s++) if(n[s]->x > x && p[s]->x==2) n[s]->check_cluster(l-1, sum_tmp, if_l);

	if(if_l) return sum_tmp;
	else     return 0;
}


void Node::check_cluster(int l_tmp, int& sum_tmp, bool& if_l){
    sum_tmp ++;
	if(l_tmp>0) {for(int s=0;s<b;s++) if(n[s]->x>x && p[s]->x==2) 	n[s]->check_cluster(l_tmp-1, sum_tmp, if_l);}
	else if_l = true;

}


/**
* This function returns true if this node is connected to the node n_tmp.
* @param n_tmp the second node
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Node::is_connected_to_Node(Node *n_tmp){
	for (int bb=0;bb<b;bb++) if(n[bb]==n_tmp) return true;
	return false;
}

/**
* This function returns true if this node is connected to the pore p_tmp
* @param p_tmp the pore
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Node::is_connected_to_Pore(Pore *p_tmp){
	for (int bb=0;bb<b;bb++) if(p[bb]==p_tmp) return true;
	return false;
}

/**
* This function switches node neighbors n_old to n_new in the list of neighboring nodes.
* @param n_old old node
* @param n_new new node
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::change_node_neighbours(Node *n_old,Node *n_new, Pore *p_tmp){
	for (int j=0; j<b;j++) if(n[j]==n_old && p[j]==p_tmp) n[j] = n_new;
}


int  Node::neighbor_number(Node *n_tmp){

	for(int s=0; s<b;s++) if(n[s]==n_tmp) return s;
	cerr<<"WARNING: The neighbor nr has not been found!!!"<<endl;
	exit(666);
	return -1;
}

/**
* This function adds a grain to the list of neighboring grains.
* @param g_tmp the grain to be added
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::add_Grain (Grain *g_tmp){
	//if grain already in the list do nothing
	for(int i=0; i<bG; i++) if(g[i]==g_tmp) return;

	Grain ** g_new = new Grain* [bG+1];
	for(int i=0; i<bG; i++) g_new[i]=g[i];
	g_new[bG] = g_tmp;
	delete [] g; g = g_new;
	bG++;
}

/**
* This function removes a grain from the list on neighboring grains.
* @param g_tmp the grain to be removed
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::remove_Grain (Grain *g_tmp){

	Grain ** g_new = new Grain * [bG];
	int b_new = 0;
	for(int bb=0;bb<bG;bb++) if(g[bb]!=g_tmp) 	g_new[b_new++] = g[bb];

	delete [] g; g = g_new;
	bG = b_new;
}

/**
* This function adds a node and a pore to the list of neighbors.
* @param n_tmp the node to be added
* @param p_tmp the pore to be added
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::add_neighbor (Node *n_tmp, Pore *p_tmp){
	//if neighbor already in the list do nothing
	for(int i=0; i<b; i++) if(n[i]==n_tmp || p[i]==p_tmp) return;

	Pore ** p_new = new Pore* [b+1];
	Node ** n_new = new Node* [b+1];

	for(int i=0; i<b; i++){ n_new[i]=n[i]; p_new[i]=p[i];}
	n_new[b] = n_tmp;
	p_new[b] = p_tmp;
	delete [] n; n = n_new;
	delete [] p; p = p_new;
	b++;
}

/**
* This function removes a node from the list of neighbors.
* @param n_tmp the node to be removed
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::remove_neighbor (Node *n_tmp){

	Pore ** p_new = new Pore * [b];
	Node ** n_new = new Node * [b];
	int b_new = 0;

	for(int bb=0;bb<b;bb++) if(n[bb]!=n_tmp) {
		n_new[b_new]  =n[bb];
		p_new[b_new++]=p[bb];
	}
	delete [] n; n = n_new;
	delete [] p; p = p_new;
	b = b_new;
}

/**
* This function removes a pore from the list of neighbors.
* @param p_tmp the pore to be removed
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::remove_neighbor (Pore *p_tmp){

	Pore ** p_new = new Pore * [b];
	Node ** n_new = new Node * [b];
	int b_new = 0;

	for(int bb=0;bb<b;bb++) if(p[bb]!=p_tmp) {
		n_new[b_new]  =n[bb];
		p_new[b_new++]=p[bb];
	}
	delete [] n; n = n_new;
	delete [] p; p = p_new;
	b = b_new;
}

/**
* This function adds a new grain to the list of neighboring grains.
* @param g_tmp the grain to be added
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Node::add_new_Grain(Grain *g_tmp){
	for(int i=0; i<bG; i++) if(g[i]==g_tmp) return;
	g[bG] = g_tmp;
	bG++;
}

bool Node::is_clogged() {
    for (int i=0;i<b;i++)
        if (p[i]->d>0)
            return false;
    return true;
}

ostream & operator << (ostream & stream, Node &n){

	stream<<"Node("<<n.a<<"): t = "<<int(n.t)<< " b = "<<n.b<<" bG = "<<n.bG<<" n = (";
	for(int k=0;k<n.b;k++){
		stream<<n.n[k]->a;
		if(k<n.b-1) stream<<",";}
	stream<<")   p = (";
	for(int k=0;k<n.b;k++) {
		stream<<n.p[k]->a;
		if(k<n.b-1) stream<<",";}
	stream<<")   g = (";
	for(int k=0;k<n.bG;k++) {
		stream<<n.g[k]->a;
		if(k<n.bG-1) stream<<",";}
	stream<<")";
  stream<<  "  u = "<<n.u << "  is-f = "<<n.is_fracture;

	return stream;
}

ofstream_txt & operator << (ofstream_txt & stream, Node &n){

	stream <<setw(12)<<n.a<<setw(12)<<setw(6)<<n.t<<n.u<<setw(12)<<n.cb<<setw(12)<<n.cc<<endl;
	return stream;
}



bool Node::can_be_calculated(){

    for (int k=0; k<b; k++)
        if(p[k]->d>0  and fabs(p[k]->q)>epsilon_for_c and n[k]->u > u and n[k]->tmp<2)
            return false;
    return true;
}

void Node::set_new_concentration(Network *N, SPECIES_NAME species){
	cerr<< "set_new_concentartion is called"<<endl;

    double Q = 0;
    double QB = 0;

    for (int k=0; k<b; k++)
    	// cerr<< "k="<<k<<endl;
        if(p[k]->d>0 and fabs(p[k]->q)>epsilon_for_c and n[k]->u > u ) {
	        Q+=fabs(p[k]->q);
        	if (species == SPECIES_NAME::B){
        	cerr<< "Calculating QB in set_new_concentration"<<endl;
        	QB += fabs(p[k]->q) * n[k]->cb * N->outlet_c_b_coeff(p[k]); // To do need to modify
        	cerr<< "QB="<<QB<<", k="<<k<<endl;
			}
            else if (species == SPECIES_NAME::C)
                QB += N->outlet_c_c_1(p[k]) + fabs(p[k]->q) * n[k]->cc * N->outlet_c_c_2_coeff(p[k]); // To do need to modify
            else std::cerr << "Unknown SPECIES_NAME." << std::endl;

        }

    if(Q!=0 and QB>0){
        if (species == SPECIES_NAME::B) cb = QB/Q;
        else if (species == SPECIES_NAME::C) cc = QB/Q;
        else std::cerr << "Unknown SPECIES_NAME." << std::endl;
    }

}


