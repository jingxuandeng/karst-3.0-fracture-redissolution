

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
	tmp = name; a = name; tmp2=0;  x=1;
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
	tmp = name; a = name; tmp2=0; x=1; is_lhs=false;
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
	tmp = g.tmp; a  = g.a; tmp2=0; x=g.x; is_lhs=g.is_lhs;
	bN  = g.bN;  bP = g.bP;

	n = new Node*[bN];
	p = new Pore*[bP];

	for (int i=0;i<bN;i++) n[i]=g.n[i];
	for (int i=0;i<bP;i++) p[i]=g.p[i];

}


Grain& Grain::operator = (Grain &g){


	Va  = g.Va;  Ve = g.Ve;  Vx = g.Vx; is_lhs=g.is_lhs;
	tmp = g.tmp; a  = g.a; tmp2=g.tmp2; x=g.x;
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
	tmp = name; a = name; tmp2=0; x=1; is_lhs=false;
	bN=bbN; bP=bbP;

	n = new Node*[bN];
	p = new Pore*[bP];

	for(int i=0;i<bN;i++) n[i] = nn0[i];
	for(int i=0;i<bP;i++) p[i] = pp0[i];
}

Grain::Grain (float name, double V_a_tmp, double V_e_tmp, double V_x_tmp, int bb_N, int bb_P){
	Va = V_a_tmp; 	Ve = V_e_tmp;   Vx = V_x_tmp;
	tmp = name; a = name; tmp2=0; x=1; is_lhs=false;
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

		V0 = sqrt(P*(P-x)*(P-y)*(P-z))*S->H_z;    //additional parameter
        double d_mean = 0;
        for (int i=0;i<bP;i++)
            if (!(p[i]->is_fracture))  d_mean += p[i]->d / bP;
            else                       d_mean += S->d0 / bP;

        //Va = (pow(3. - (pow(S->d0,2)*M_PI)/Va_0,1.5)*Va_0)/(3.*sqrt(3));   //Felerny błąd ze wzorem \Delta d/ 2 : Wzór stary, dla dwuch ruszających się poprzeczek!!! (dla rombów i 1D itp.)
        //Va = (pow(2. - (pow(d_mean,2)*M_PI)/Va_0,1.5)*Va_0)/(2.*sqrt(2.)); //Felerny błąd ze wzorem \Delta d/ 2 : Poprawiony wzór na wszystkie poprzeczki ruszające się.



			//Va = (pow(6. - (pow(S->d0,2)*M_PI)/V0,1.5)*V0)/(6.*sqrt(6));   // the best 1D version, for half diamonds
			Va = (pow(4. - (pow(d_mean, 2) * M_PI) / V0, 1.5) * V0) / (8); //generic 2D version, all pores are changing

	}

	//WARNING: the general formula should be implemented for cubic network with added random node positions
	else if(bN==8){
		Va = 1; V0=1;
		for (int i=0;i<bP;i++) Va -= (M_PI*p[i]->d*p[i]->d*p[i]->l/4)/4;
	}

	else if (bN==4){

		double x = S->node_distance(n[0], n[1]);
		double y = S->node_distance(n[1], n[2]);
		double z = S->node_distance(n[2], n[0]);
		double P = (x+y+z)/2;

		Va = sqrt(P*(P-x)*(P-y)*(P-z))*S->H_z;

		x = S->node_distance(n[0], n[3]);
		y = S->node_distance(n[3], n[2]);
		z = S->node_distance(n[2], n[0]);
		P = (x+y+z)/2;

		Va += sqrt(P*(P-x)*(P-y)*(P-z))*S->H_z;
        V0=Va;


		double d_mean = 0;
		for (int i=0;i<bP;i++)
			if (!(p[i]->is_fracture))  d_mean += p[i]->d / bP;
			else                       d_mean += S->d0 / bP;


		Va  = 2* (pow(6. - (pow(S->d0,2)*M_PI)/(V0/2),1.5)*(V0/2))/(6.*sqrt(6));
		//cerr<<"I calculated initial volume for diamonds."<<endl;
	}


	else if (bN==2){
		Va = 0; V0=0;
	}



	else {
		cerr<<"WARNING: Initial grain volume for general network has to be implemented!!!"<<endl;
		Va = 1; V0=1;
		cerr<<*this<<endl;
	}

	if(Va<0) Va = 0;
    if(!V0>0) V0 = sqrt(3)/4*S->H_z;
	if(!(Va>=0)){
		Va=sqrt(3)/4*S->H_z;
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

		Va_tmp = sqrt(P*(P-x)*(P-y)*(P-z))*S->H_z;
	}

	//WARNING: the general formula should be implemented for cubic network with added random node positions
	else if(bN==8){
		Va_tmp = 1;
	}

	else if (bN<3){
        Va_tmp =  2*(S->l0*S->H_z*sqrt(3.)/4.);//*(n[0]->xy - n[1]->xy);

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
int Grain::to_be_merge(Network *S){

    int if_fracture = 0;
    for (int i=0; i<bP;i++)
            if (p[i]->is_fracture) if_fracture++;

    //if (if_fracture==0)               return false;   // this was hindering horizontal channels
    if((Va+Ve+Vx)/V0 < S->merge_factor or Va+Ve+Vx<=0) {
        if (if_fracture > 1) return 2;
        else return 1;
    }
    else
        return 0;
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
	if (bN<2) 	    return true; // fixme: Before (bN==2&& bP>1)  PROBLEM: if bN<2 (if we allow two nodes grains the mass conservation is violated, maybe some problem with calculating u or c filed)
	if (bP<2)		return true; //fixme: Before bP ==0
	if(bN == 2) if(n[0]->xy-n[1]->xy > 3.) return true;     //avoid long thin grains
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

    cerr<<"Setting effective d and l for "<<*p_master<<endl;

	double s = 0; 			//reaction surface in a grain
	double r = 0;			//effective resistance in a grain

	double q_max = 0;		//flow in a grain (calculated for u_max)
	double q_min = 0;		//flow in a grain (calculated for u_min)
    int   is_there_a_fracture = 0;
	Node *n_min=NULL;		//node with u_min
	Node *n_max=NULL;		//node with u_max

	//choosing n and u max/min

    if (p_master->n[0]->u < p_master->n[1]->u) {n_min = p_master->n[0]; n_max = p_master->n[1];}
    else                                       {n_min = p_master->n[1]; n_max = p_master->n[0];}

	//calculating s and q_max/min
	for(int b=0;b<bP;b++){
        bool pipe_formula = (!(S->sandwich_pores and p[b]->is_fracture) and (p[b]->d<S->H_z or S->no_max_z));
		if(pipe_formula)
            s+=p[b]->l*p[b]->d*M_PI;  //FIXME: tu był błąd bo nie obsługiwałam wcześńiej pora - kanapki, teraz jest ok
		else
            s+=p[b]->l*M_PI;
        if(p[b]->n[0]==n_max||p[b]->n[1]==n_max) q_max += fabs(p[b]->q);
		if(p[b]->n[0]==n_min||p[b]->n[1]==n_min) q_min += fabs(p[b]->q);
        if(p[b]->is_fracture) is_there_a_fracture++;
	}

	if(q_max+q_min <=0) {
		cerr<<"ERROR: Problem with calculating q_max and q_min in setting effective d and l."<<endl<<*this<<endl;
		cerr<<"p_master : "<<*p_master<<endl;
        cerr<<"q_max = "<<q_max<<endl;
		cerr<<"q_min = "<<q_min<<endl<<endl;
        cerr<<"n_min = "<<*n_min<<endl;
        cerr<<"n_max = "<<*n_max<<endl;
		for(int b=0;b<bP;b++) cerr<<*p[b]<<endl;
        //exit(876);
	}
    double delta_P=0;
    if(n_max->u - n_min->u > 0) delta_P = n_max->u - n_min->u;
    else{
        for(int b=0;b<bN;b++)
            if(n[b]->u != n_max->u)
                delta_P = fabs(n[b]->u - n_max->u);
    }
    if(delta_P==0) {cerr<<"ERROR: problem with pressure drop in the grain while merging."<<endl; }//exit(123321);}
    if(delta_P>0 and q_max+q_min >0)
        r = delta_P/((q_max+q_min)/2.);

	if(s<=0 || r<=0) cerr<<"ERROR: Problem with calculating r and s in setting effective d and l."<<endl<<*this<<endl;
	//calculating new d and l
	if(r>0 && s>0) {
        bool pipe_formula = (!(S->sandwich_pores and p_master->is_fracture) and (r>(128*S->mu_0*p_master->l)/(M_PI*pow(S->H_z,4)) or S->no_max_z));
        if(pipe_formula) {
            //old formulas for cylinger
            p_master->d = 2. * pow((4. * S->mu_0 * s) / (M_PI * M_PI * r), 0.2);
            p_master->l = 0.5 * pow((r * pow(s, 4)) / (4. * S->mu_0 * pow(M_PI, 3)), 0.2);
        }
        else {
            if(p_master->d<S->H_z){
                //new formulas for a thin fracture
//                p_master->d = 4. * pow((2. * S->mu_0 ) / ( M_PI * r), 1./3);
//                p_master->l = 1;
                p_master->l = fmax(fmin(s/M_PI,1.0),p_master->l);
                p_master->d = 4. * pow((2. * S->mu_0 * p_master->l) / (  M_PI * r), 1./3);

            }
            else{
                //new formulas for a fracture
//                p_master->d = (128 * S->mu_0) / (M_PI * r);
//                p_master->l = 1;
                p_master->l = fmax(fmin(s/M_PI,1.0),p_master->l);
                p_master->d = (128 *  p_master->l * S->mu_0) / ( M_PI * r);



            }
        }

        if(p_master->d <= S->d_min || p_master->l <= S->l_min) {
            cerr << "WARNING: P_master has problem with d = " << p_master->d << "   l = " << p_master->d << endl;
            p_master->d=0;
            //S->clear_unconeccted_pores();
        }
    }

    if(S->if_verbose) cerr<<"New  d and l for "<<*p_master<<endl;
    if(is_there_a_fracture==1) p_master->is_fracture = true;           //make sure the pore is tagged as a part of a fracture but do not create false fracture pores
    //else p_master->is_fracture = false;
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
* This function prints grain properties in one line
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


