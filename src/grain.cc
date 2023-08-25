

#include "grain.h"
#include "network.h"




/**
* @brief Constructor for defoult case.
* @author Agnieszka Budek
* @date 13/03/2020
*/
Grain::Grain (){

	Va = 0; 	Ve=0;     Vx=0;
	tmp = -1; a = -1; tmp2=-1;
	bN=0; bP=0;

	n = NULL;
	p = NULL;
}

/**
* @brief Constructor for triangular case.
* This is a constructor for triangular grain.
* @author Agnieszka Budek
* @param name initial name
* @param nn0 first node
* @param nn1 second node
* @param nn2 third node
* @date 25/09/2019
*/
Grain::Grain (float name, Node* nn0, Node* nn1, Node* nn2){
		
	Va = 0; 	Ve=0;     Vx=0;
	tmp = name; a = name; tmp2=0;
	bN=3; bP=3;

	n = new Node*[3];
	p = new Pore*[3];

	n[0]=nn0; n[1]=nn1; n[2]=nn2;	
		
	p[0]=findPore(n[0],n[1]);
	p[1]=findPore(n[1],n[2]);
	p[2]=findPore(n[2],n[0]);
}


/**
* @brief Constructor for diamondal case.
* This is a constructor for square grain.
* @author Agnieszka Budek
* @param name initial name
* @param nn0 first node
* @param nn1 second node
* @param nn2 third node
* @param nn3 forth node
* @date 13/03/2020
*/
Grain::Grain (float name, Node* nn0, Node* nn1, Node* nn2,Node *nn3){

	Va = 0; 	Ve=0;      Vx=0;
	tmp = name; a = name; tmp2=0;
	bN=4; bP=4;

	n = new Node*[4];
	p = new Pore*[4];

	n[0]=nn0; n[1]=nn1; n[2]=nn2; n[3]=nn3;

	p[0]=findPore(n[0],n[1]);
	p[1]=findPore(n[1],n[2]);
	p[2]=findPore(n[2],n[3]);
	p[3]=findPore(n[3],n[0]);
}

/**
* This is a basic constructor for general grain.
* @brief Constructor for general case.
* @author Agnieszka Budek
* @param g parent grain
* @date 25/09/2019
*/
Grain::Grain (Grain &g){

	Va  = g.Va;  Ve = g.Ve;   Vx = g.Vx;
	tmp = g.tmp; a  = g.a; tmp2=0;
	bN  = g.bN;  bP = g.bP;

	n = new Node*[bN];
	p = new Pore*[bP];

	for (int i=0;i<bN;i++) n[i]=g.n[i];
	for (int i=0;i<bP;i++) p[i]=g.p[i];

}


Grain& Grain::operator = (Grain &g){


	Va  = g.Va;  Ve = g.Ve;  Vx = g.Vx;
	tmp = g.tmp; a  = g.a; tmp2=g.tmp2;
	bN  = g.bN;  bP = g.bP;

	n = new Node*[bN];
	p = new Pore*[bP];

	for (int i=0;i<bN;i++) n[i]=g.n[i];
	for (int i=0;i<bP;i++) p[i]=g.p[i];
	return *this;

}

/**
* This is a basic constructor for general grain
* @author Agnieszka Budek
* @date 25/09/2019
*/
Grain::Grain (float name, int bbP, int bbN, Node** nn0, Pore** pp0){

	Va = 0; 	Ve=0;     Vx=0;
	tmp = name; a = name; tmp2=0;
	bN=bbN; bP=bbP;

	n = new Node*[bN];
	p = new Pore*[bP];

	for(int i=0;i<bN;i++) n[i] = nn0[i];
	for(int i=0;i<bP;i++) p[i] = pp0[i];
}

Grain::Grain (float name, double V_a_tmp, double V_e_tmp, double V_x_tmp, int bb_N, int bb_P){
	Va = V_a_tmp; 	Ve = V_e_tmp;   Vx = V_x_tmp;
	tmp = name; a = name; tmp2=0;
	bN=bb_N; bP=bb_P;

	n = new Node*[bN];
	p = new Pore*[bP];

	for(int i=0;i<bN;i++) n[i] = NULL;
	for(int i=0;i<bP;i++) p[i] = NULL;

}

/**
* In this function an initial volume of a grain is calculated.
* WARNING: Up to now only the trianglular volume is implemented.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Grain::calculate_initial_volume (Network *S){


	if (bN==3 ){ //version for 2D, triangular network with Heron formula
		double x = S->node_distance(n[0], n[1]);
		double y = S->node_distance(n[1], n[2]);
		double z = S->node_distance(n[2], n[0]);
		double P = (x+y+z)/2;

		Va = sqrt(P*(P-x)*(P-y)*(P-z));
		for (int i=0;i<bP;i++) Va -= (M_PI*p[i]->d*p[i]->d*p[i]->l/4)/2;
	}
	//WARNING: the general formula should be implemented for cubic network with added random node positions
	else if(bN==8){
		Va = 1;
		for (int i=0;i<bP;i++) Va -= (M_PI*p[i]->d*p[i]->d*p[i]->l/4)/4;
	}

	else if (bN==4){

		double x = S->node_distance(n[0], n[1]);
		double y = S->node_distance(n[1], n[2]);
		double z = S->node_distance(n[2], n[0]);
		double P = (x+y+z)/2;

		Va = sqrt(P*(P-x)*(P-y)*(P-z));

		x = S->node_distance(n[0], n[3]);
		y = S->node_distance(n[3], n[2]);
		z = S->node_distance(n[2], n[0]);
		P = (x+y+z)/2;

		Va += sqrt(P*(P-x)*(P-y)*(P-z));

		for (int i=0;i<bP;i++) Va -= (M_PI*p[i]->d*p[i]->d*p[i]->l/4)/2;
	}


	else if (bN==2){
		Va = 0;
	}



	else {
		cerr<<"WARNING: Initial grain volume for general network has to be implemented!!!"<<endl;
		Va = 1;
		cerr<<*this<<endl;
	}

	if(Va<0) Va = 0;
	if(!(Va>=0)){
		Va=sqrt(3)/4;
		//cerr<<"WARNING: Problem with calculating initial volume for grain"<<*this<<" Va = "<<Va<<endl;
		//cerr<<"Problematic pores: "<<" p = ("<<p[0]->l<<","<<p[1]->l<<","<<p[2]->l<<")"<<endl;
	}

	//updating Vx percentage
	if(S->Vx_perc > 0){
		Vx = S->Vx_perc*Va;
		Va = (1-S->Vx_perc)*Va;
	}
}

/**
* In this function maximal possible volume of a grain is calculated.
* This volume corresponds to the case in which the surrounding pores has d=0.
* WARNING: Up to now only the trianglular volume is implemented.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Grain::calculate_maximal_volume (Network *S){


	//cerr<<"Calculating maximal volume..."<<endl;
	double Va_tmp;

	if(S->initial_xy == NULL) {cerr<<"ERROR: table initial_xy has not been initialized!"<<endl; exit(666);}
	Point * pp = S->initial_xy[a];


	if (bN==3 ){ //version for 2D, triangular network with Heron formula
		double x = S->point_distance(pp[0], pp[1]);
		double y = S->point_distance(pp[1], pp[2]);
		double z = S->point_distance(pp[2], pp[0]);
		double P = (x+y+z)/2;

		Va_tmp = sqrt(P*(P-x)*(P-y)*(P-z));
	}

	//WARNING: the general formula should be implemented for cubic network with added random node positions
	else if(bN==8){
		Va_tmp = 1;
	}

	else if (bN==2){
		Va_tmp = 0;
	}

	else if (bN==4){



		double x = S->point_distance(pp[0], pp[1]);
		double y = S->point_distance(pp[1], pp[2]);
		double z = S->point_distance(pp[2], pp[0]);
		double P = (x+y+z)/2;

		Va_tmp = sqrt(P*(P-x)*(P-y)*(P-z));

		x = S->point_distance(pp[0], pp[3]);
		y = S->point_distance(pp[3], pp[2]);
		z = S->point_distance(pp[2], pp[0]);
		P = (x+y+z)/2;

		Va_tmp += sqrt(P*(P-x)*(P-y)*(P-z));
	}

	else {
		cerr<<"WARNING: Initial grain volume for general network has to be implemented!!!"<<endl;
		Va_tmp = 1;
		cerr<<*this<<endl;
	}

	if(Va_tmp<0) Va_tmp = 0;
	if(!(Va_tmp>=0)){
		Va_tmp=sqrt(3)/4;
	}

	return Va_tmp;
}

/**
* In this function the conditions for grain merging is checked.
* Now the condition checks if there is a material in the grain.
* This condition can be changed later.
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Grain::to_be_merge(){
	if	(Va+Ve+Vx<=0)   	 		return true;
	else						return false;
}


/**
* In this function it is checked either the grain is pathologiacal or not.
* Some grains are merged due to the topology problems, i.e. those with no pores.
* This condition can be changed later.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Grain::is_pathological (){		//condition for pathological grain, this kind of grain is merged automatically
	if (bN==2 && bP>1) 	return true;
	if (bP==0)			return true;
	return false;
}

/**
* This function checks if the node is connected to the grain.
*
* @param n_tmp node to be checked
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Grain::do_contain_Node(Node *n_tmp){
	for (int b=0;b<bN;b++) if(n[b]==n_tmp) return true;
	return false;
}


/**
* This function checks if the pore is connected to the grain.
*
* @param p_tmp pore to be checked
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Grain::do_contain_Pore(Pore *p_tmp){
	for (int b=0;b<bP;b++) if(p[b]==p_tmp) return true;
	return false;
}


/**
* This function removes pore form pore list.
*
* @param p_tmp pore to be removed
* @author Agnieszka Budek
* @date 25/09/2019
*/
void  Grain::remove_Pore(Pore *pp){

	Pore ** p_new = new Pore * [bP];
	int b_new = 0;
	for(int bb=0;bb<bP;bb++) if(p[bb]!=pp) 	p_new[b_new++] = p[bb];

	delete [] p; p = p_new;
	bP = b_new;

}



/**
* This function adds pore to the pore list
*
* @param p_tmp pore to be added
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Grain::add_Pore   (Pore *p_tmp){
	//if pore already in the list do nothing
	for(int i=0; i<bP; i++) if(p[i]==p_tmp) return;

	Pore ** p_new = new Pore* [bP+1];
	for(int i=0; i<bP; i++) p_new[i]=p[i];
	p_new[bP] = p_tmp;
	delete [] p; p = p_new;
	bP++;
}


/**
* This function removes node form the node list.
*
* @param n_tmp node to be removed
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Grain::remove_Node(Node *n_tmp){
	Node ** n_new = new Node * [bN];
	int b_new = 0;
	for(int bb=0;bb<bN;bb++) if(n[bb]!=n_tmp) 	n_new[b_new++] = n[bb];

	delete [] n; n = n_new;
	bN = b_new;
}

/**
* This function adds node to the node list
*
* @param n_tmp node to be added
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Grain::add_Node   (Node *n_tmp){
	//if node already in the list do nothing
	for(int i=0; i<bN; i++) if(n[i]==n_tmp) return;

	Node ** n_new = new Node* [bN+1];
	for(int i=0; i<bN; i++) n_new[i]=n[i];
	n_new[bN] = n_tmp;
	delete [] n; n = n_new;
	bN++;
}


/**
* This function changes node n_old to node n_new in the node list.
*
* @param n_old node to be removed
* @param n_new node to be added
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Grain::change_nodes(Node *n_old,Node * n_new){

	for (int i=0; i<bN;i++) if(n[i]==n_new) {  //if new node already in the list of nodes
		for (int j=0; j<bN;j++) if(n[j]==n_old){
			n[j] = n[bN-1];
			n[bN-1]=NULL;
			bN--;
			break;
		}
		return;
	}
	for (int j=0; j<bN;j++) if(n[j]==n_old) {n[j] = n_new; return;}   //if new node not in the list of nodes
	cerr<<"WARNING: Problem with changing nodes in a grain: n = "<<*n_old<<endl<<"  g = "<<*this<<endl;
}

/**
* This function changes pore p_old to pore p_new in the pore list.
* @param p_old pore to be removed
* @param p_new pore to be added
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Grain::change_pores(Pore *p_old,Pore * p_new){

	for (int i=0; i<bP;i++) if(p[i]==p_new) {  //if new pore already in the list of nodes
		for (int j=0; j<bP;j++) if(p[j]==p_old){
			p[j] = p[bP-1];
			p[bP-1]=NULL;
			bP--;
			break;
		}
		return;
	}
	for (int j=0; j<bP;j++) if(p[j]==p_old) {p[j] = p_new; return;}    //if new pore not in the list of nodes
	cerr<<"WARNING: Problem with changing pores in a grain: p = "<<*p_old<<endl<<"  g = "<<*this<<endl;
}


/**
* This function sets new d and l for a, so called, master pore, the one that remains after merging.
* New d and l are calculated to conserve total reaction area and resistance of a system after merging.
* The first one can be calculated when the second one is more problematic and can be improved later.
*
* @p_master pore to have d and l updated
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Grain::set_effective_d_and_l(Pore *p_master,Network *S){
	double s = 0; 			//reaction surface in a grain
	double r = 0;			//effective resistance in a grain
	double u_min = 1e10; 	//the smallest pressure in a grain
	double u_max = 0;		//the largest pressure in a grain
	double q_max = 0;		//flow in a grain (calculated for u_max)
	double q_min = 0;		//flow in a grain (calculated for u_min)
	Node *n_min=NULL;		//node with u_min
	Node *n_max=NULL;		//node with u_max

	//choosing n and u max/min
	for(int b=0;b<bN;b++){
		if (n[b]->u <= u_min) {n_min = n[b]; u_min = n[b]->u;}
		if (n[b]->u >= u_max) {n_max = n[b]; u_max = n[b]->u;}
	}
	//calculating s and q_max/min
	for(int b=0;b<bP;b++){
		s+=p[b]->l*p[b]->d*M_PI;
		if(p[b]->n[0]==n_max) q_max += fabs(p[b]->q);
		if(p[b]->n[0]==n_min) q_min += fabs(p[b]->q);
	}

	if(q_max+q_min <=0) {
		cerr<<"ERROR: Problem with calculating q_max and q_min in setting effective d and l."<<endl<<*this<<endl;
		cerr<<"q_max = "<<q_max<<endl;
		cerr<<"q_min = "<<q_min<<endl<<endl;
		for(int b=0;b<bP;b++) cerr<<*p[b]<<endl;

	}

	else   r = (u_max - u_min)/(q_max+q_min)/2;

	if(s<0 || r<0) cerr<<"ERROR: Problem with calculating r and s in setting effective d and l."<<endl<<*this<<endl;
	//calculating new d and l
	if(r>0 && s>0) {
		p_master->d = 2.*pow((4.*S->mu_0*s)/(M_PI*M_PI*r),0.2);
		p_master->l = 0.5*pow( (r*pow(s,4))/(4.*S->mu_0*pow(M_PI,3)),0.2);
	}

}

/**
* This function prints grain description.
*
* @param stream output stream
* @param g grain to be printed
* @author Agnieszka Budek
* @date 25/09/2019
*/
ostream& operator << (ostream & stream, Grain &g){
	

	stream<<"Grain("<<g.a<<"): n = (";
	for(int b=0;b<g.bN;b++) {
		stream<<g.n[b]->a;
		if(b<g.bN-1) stream<<",";
	}
	stream<<")  p = (";
	for(int b=0;b<g.bP;b++){
		stream<<g.p[b]->a;
		if(b<g.bP-1) stream<<",";
	}
	stream<<")";
	stream<<"  tmp = " << g.tmp;
	stream<<"  Va = "  << g.Va;
	stream<<"  Ve = "  << g.Ve;
	stream<<"  Vx = "  << g.Vx;
	return stream;
}

/**
* This function prints grain propoerties in one line
*
* @param stream output stream
* @param g grain to be printed
* @author Agnieszka Budek
* @date 25/09/2019
*/
ofstream_txt & operator << (ofstream_txt & stream, Grain &g){
	stream <<setw(12)<<g.a<<setw(12)<<g.Va<<setw(12)<<g.Ve<<setw(12)<<g.Vx<<endl;
	return stream;
}


