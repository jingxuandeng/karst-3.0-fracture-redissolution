#include <random>
#include "network.h"
#include "triangulation.h"

//#define VISUALIZATION //only for debugging

#ifdef VISUALIZATION
#include <SFML/Graphics.hpp>
#endif


/**
* This function set the values of global mean diameter (d0) and mean pore length (l0)
*  taking into account the initial topology of the network and pore properties.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_initial_d0_and_l0 (){

	double d_tmp = 0, l_tmp = 0;
	int    i_tmp = 0;
	for(int i=0;i<NP;i++) if(p[i]->d>0) {d_tmp+=p[i]->d; l_tmp+=p[i]->l; i_tmp++;}
	d0 = d_tmp/i_tmp;
	l0 = l_tmp/i_tmp;
	mu_0 = M_PI*pow(d0,4)/(128.*l0);
	cerr<< "Initial d0 = "<<d0<<endl;
	cerr<< "Initial l0 = "<<l0<<endl;
	cerr<< "Setting mu_0 = "<<mu_0<<endl;
}

/**
* This function calculates initial mean flow through the system.
* It is important to calculate local Da_eff and G correctly.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_initial_mean_flow(){   //important for physical parameter calculations

	if (Q_tot!=0){    //constant flow through the system
		int N_inlet_pores=0;
		for(int i=0;i<N_wi;i++){
			Node* n_tmp = wi[i];
			for (int j=0; j<n_tmp->b;j++) if (n_tmp->p[j]->d>0) N_inlet_pores++;
		}
		q_in_0=Q_tot/N_inlet_pores;
	}
	else{			  //constant pressure in the system
		double Q_tot_tmp = 0;
		int N_inlet_pores=0;
		calculate_pressures();
		calculate_flows();
		for(int i=0;i<N_wi;i++){
			Node* n_tmp = wi[i];
			for (int j=0; j<n_tmp->b;j++) if (n_tmp->p[j]->d>0) {N_inlet_pores++; Q_tot_tmp += fabs(n_tmp->p[j]->q);}
		}
		q_in_0=Q_tot_tmp/N_inlet_pores;
	}
	cerr<<"Initial mean Q in inlet pore has been calculated, q_in_0 = "<<q_in_0<<endl;
}


/**
* This function calculates the initial total volume of species A
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network :: calculate_initial_total_Va(){
	Va_tot = 0;   Vx_tot = 0;
	for(int i=0;i<NG;i++) Va_tot+=g[i]->Va;
	for(int i=0;i<NG;i++) Vx_tot+=g[i]->Vx;

}

/**
* This function calculates the initial total volume of species E
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network :: calculate_initial_total_Ve(){
	Ve_tot = 0;
	for(int i=0;i<NG;i++) Ve_tot+=g[i]->Ve;

}



/**
* This function generates the 2D hexagonal network.
*
* @param N network width
* @param M network length
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::createHexagonalNetwork(int N, int M){

	if(M%2==1) {cerr<<"ERROR: N_y must be even in hexagonal network!"<<endl; exit(123);}

	cerr<<"Creating hexagonal network..."<<endl;

	if (random_seed==-1) srand (time(NULL));
	else                 srand (random_seed);
	cerr<<"Random seed: "<<random_seed<<endl;

	NN = N*M;
	NP = 3*N*M;
	NG = 2*N*M;

	b_max = 6;

	n = new Node*[NN];	p = new Pore*[NP];	g = new Grain*[NG];

	//random distribution
	default_random_engine gen;
	normal_distribution<double> N_dist(d0,d0/1); //distribution of initial diameters

	//pore and node generation
	for(int j=0;j<M;j++) for(int i=0;i<N;i++) 						n[N*j+i] 		= new Node(6,N*j+i);
	for(int j=0;j<M;j++) for(int i=0;i<N;i++) for(int k=0;k<3;k++) 	p[(N*j+i)*3+k] 	= new Pore(d0,l0,(N*j+i)*3+k);
	//for(int j=0;j<M;j++) for(int i=0;i<N;i++) for(int k=0;k<3;k++) 	p[(N*j+i)*3+k] 	= new Pore(N_dist(gen),l0,(N*j+i)*3+k);

	//pore and node connections
	for(int j=0;j<M;j++) for(int i=0;i<N;i++){
		Node* nn = n[N*j+i];

		nn->n[0] = n[N*((j+1)%M) 	+ (i-j%2+N)%N]; 	nn->p[0] = p[(N*j+i)*3+0];
		nn->n[1] = n[N*((j+1)%M) 	+ (i-j%2+1)%N]; 	nn->p[1] = p[(N*j+i)*3+1];
		nn->n[2] = n[N*j  			+ (i+1+N)%N]; 		nn->p[2] = p[(N*j+i)*3+2];
		nn->n[3] = n[N*((j-1+M)%M)  + (i-j%2+1)%N];		nn->p[3] = p[(N*((j-1+M)%M) + (i-j%2+1)%N)*3+0];
		nn->n[4] = n[N*((j-1+M)%M)  + (i-j%2+N)%N]; 	nn->p[4] = p[(N*((j-1+M)%M) + (i-j%2+N)%N)*3+1];
		nn->n[5] = n[N*j  			+ (i-1+N)%N];		nn->p[5] = p[(N*j+(i-1+N)%N)*3+2];


		p[(N*j+i)*3+0]->n[0] = nn; 	p[(N*j+i)*3+0]->n[1] = nn->n[0];
		p[(N*j+i)*3+1]->n[0] = nn; 	p[(N*j+i)*3+1]->n[1] = nn->n[1];
		p[(N*j+i)*3+2]->n[0] = nn; 	p[(N*j+i)*3+2]->n[1] = nn->n[2];
	}


	//grain generation
	for(int i=0;i<N;i++) for(int j=0;j<M;j++){
		Node* nn = n[N*j+i];
		g[(N*j+i)*2+0] = new Grain((N*j+i)*2+0,nn,nn->n[0],nn->n[1]);
		g[(N*j+i)*2+1] = new Grain((N*j+i)*2+1,nn,nn->n[1],nn->n[2]);
	}


	//connections between pores and grains
	if(if_verbose) cerr<<"Adding info about connections between pores and grains."<<endl;
	for(int i=0;i<NG;i++) for(int k=0;k<3;k++){
		if(g[i]->p[k]==NULL) {cerr<<"ERROR: Problem with generation hex network, g["<<i<<"] -> p["<<k<<"] = NULL"<<endl; continue;};
		if(g[i]->p[k]->g[0]==g[i] || g[i]->p[k]->g[1]==g[i])  continue;
		if(g[i]->p[k]->g[0]==NULL) { g[i]->p[k]->g[0]=g[i];   continue; }
		if(g[i]->p[k]->g[1]==NULL) { g[i]->p[k]->g[1]=g[i];   continue; }
		cerr<<"Problem with matching nodes and grains."<<endl;
	}

	if(if_verbose){  //for debugging
		cerr<<endl<<"List of  Pores:"<<endl;
		for(int i=0;i<NP;i++) cerr<<*p[i]<<endl;
		cerr<<endl<<"List of Nodes:"<<endl;
		for(int i=0;i<NN;i++) cerr<<*n[i]<<endl;
		cerr<<endl<<"List of Grains:"<<endl;
		for(int i=0;i<NG;i++) cerr<<*g[i]<<endl;
	}

	//adding nodes position
	if(if_verbose) cerr<<"Adding info about node positions..."<<endl;
	for(int i=0;i<N;i++) for(int j=0;j<M;j++){  n[N*j+i]->xy = Point(i-j%2*0.5,j*sqrt(3)/2);}

	//inlet and outlet
	if(if_verbose) cerr<<"Initializing input and output..."<<endl;

	wi = new Node* [N];
	wo = new Node* [N];
	N_wi=N; N_wo=N;

	for(int i=0;i<N;i++) {
		wi[i] = n[i];		wi[i]->t = 1;	wi[i]->tmp = - wi[i]->tmp-1;
		wo[i] = n[N*M-1-i];	wo[i]->t =-1;	wo[i]->tmp = - wo[i]->tmp-1;
	}
	for(int i=0;i<N;i++) for(int k=0;k<wi[i]->b;k++) if(wi[i]->p[k]->n[0]->tmp<0 && wi[i]->p[k]->n[1]->tmp<0) wi[i]->p[k]->d=0;
	for(int i=0;i<N;i++) for(int k=0;k<wo[i]->b;k++) if(wo[i]->p[k]->n[0]->tmp<0 && wo[i]->p[k]->n[1]->tmp<0) wo[i]->p[k]->d=0;


	//setting initial pore length and optionally adding randomness and
	if(if_randomness_in_regular_net) add_randomness_to_regular_network(gauss_sigma_d,max_rand_shift_xy);
	save_info_about_initial_node_positions ();
	//updating pore lengths
	for(int i=0; i<NP;i++) p[i]->l = point_distance(p[i]->n[0]->xy, p[i]->n[1]->xy);


	//optional barrier (resignation form periodic boundary conditions)
	if  (!if_periodic_bc)  for(int i=0; i<NP; i++) if(p[i]->n[0]->xy - p[i]->n[1]->xy > 5*l0) p[i]->d=0; //old version, leaving unnecessary grains


	for(int i=0;i<NN;i++) if(n[i]->a != i) cerr<<"Warning: problem with nodes names:  " <<i<<"  "<<n[i]->a<<endl;
	for(int i=0;i<NP;i++) if(p[i]->a != i) cerr<<"Warning: problem with pores names:  " <<i<<"  "<<p[i]->a<<endl;
	for(int i=0;i<NG;i++) if(g[i]->a != i) cerr<<"Warning: problem with grains names: " <<i<<"  "<<g[i]->a<<endl;



//adding information about grains to nodes
	if(if_track_grains) add_information_about_grains_in_nodes();


	//saving network topology
	export_topology_file ("topology_befor_cleaning.out");

	NN_max = NN; NP_max = NP; NG_max = NG;


	check_network_connections();
	if(if_clear_unused_pores){
		//Attention: This step before info of grains in nodes in generated
		//clearing nodes and pores
		cerr<<"Clearing unused pores:"<<endl;
		clear_unused_pores();
		cerr<<"Clearing unused nodes:"<<endl;
		clear_unused_nodes();
		cerr<<"Clearing unused grains:"<<endl;
		clear_unused_grains();
		cerr<<endl;
		check_network_connections();
		}


	cerr<<"Hexagonal network has been created."<<endl;

}

/**
* This function shifts the position of nodes and/or randomize pore diameters to add some randomness to the network.
*
* @param d_sigma sigma in Gauss distribution for pore diameters (mean always is equal d0)
* @param d_sigma max_nodes_shift maximal shift of node positions (in l0 units)
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::add_randomness_to_regular_network(double d_sigma, double max_nodes_shift){

	cerr<<"Adding randomness to hexagonal network..."<<endl;

	if(fabs(d_sigma)>0){  //calculating new, random diameters
		int i=NP-1;
		while(i>=0){
			 double X_1 = (double)rand()/RAND_MAX;
			 double X_2 = (double)rand()/RAND_MAX;
			 double tmp = sqrt(-2*log(X_1))*cos(2*M_PI*X_2);
			 if(p[i]->d!=0){
				 if(d_sigma<0) p[i]->d = d0*exp(fabs(d_sigma)*tmp);   //log normal distribution
				 else          p[i]->d = d_sigma*tmp + d0;            //gaussian distribution
			 }

			 if(p[i]->d>=0 && p[i]->d<10*d0) i--;}  //WARNING: The pore diameter can not exceed 10xd0
	}
	if(max_nodes_shift>0){ //adding shift to node positions
		for(int i=NN-1; i>=0;i--){ //new nodes position
			double dx_tmp = 0.9*max_nodes_shift*(double)rand()/RAND_MAX-0.45;
			double dy_tmp = max_nodes_shift*(sqrt(3)/2)*(0.9*(double)rand()/RAND_MAX-0.45);
			n[i]->xy.x+=dx_tmp;  n[i]->xy.y+=dy_tmp;
			if(type_of_topology == "cubic"){
				double dz_tmp = max_nodes_shift*(sqrt(3)/2)*(0.9*(double)rand()/RAND_MAX-0.45);
				n[i]->xy.z+=dz_tmp;
			}
		}
	}

	if(if_debugging_printing){
		for(int i=0;i<NN;i++) n[i]->tmp = n[i]->u;
		for(int i=0;i<NP;i++) p[i]->tmp = p[i]->l;
		description_note = "After adding randomness to the network: s = " + to_string(tot_steps);
		net_ps<<*this;}
}

/**
* This function generates random network using deloney traingulation.
* It is at some level comple since we want deloney propoerties and periodic boundary conditions.
*
* @param N_x network width
* @param N_y network length
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network:: createRandomTrianglesNetwork(int N_x, int N_y){

	double eps =1.;  //parameter for regular network threshold

	bool if_regular_points = false;

	N_wi = 0; N_wo = 0;
	NN = N_x*N_y; //total number of nodes

	cerr<<endl<<"Generating network from triangulation:"<<endl;

	// for communication with traingulation package
	std::vector<Vector2 <double> >  points;				//all points (repeting ones to periodic b.c.)
	std::vector<Vector2 <double> >  points_tmp;			//correct points (rest to set periodic b. c.)
	std::vector<Triangle<double> >  triangles; 			//all triangles
	std::vector<Edge    <double> >  edges;				//all edges


	triangulation(N_x, N_y, points, points_tmp, edges, triangles, if_regular_points, if_periodic_bc,random_seed);
	cerr<<"Tiangulation for all points done."<<endl;


	//export data form points,triangles and edges

	NN = N_x*N_y; NP = 0; NG = 0;
	n = new Node *[NN];
	cerr<<"Number of nodes: NN = "<<NN<<endl;

	cerr<<"Getting rid of unnecessary edges."<<endl;
	vector<Edge<double> > edges2;
	for(auto e : edges){
		bool if_new = true;
		for(auto e2 : edges2) if(e==e2)  if_new = false;
		if(if_new) if((e.p1.a<e.p2.a && e.p1.a >=0 && e.p1.a < NN)) edges2.push_back(e);
	}
	edges2.swap(edges);


	//counting nr of neighbors
	cerr<<"Counting nr of neighbors..."<<endl;
	int * tmp1 = new int [NN]; for(int i=0;i<NN;i++) tmp1[i] = 0;

	int NP_tmp = 0;
	for (const auto &e : edges) {
		tmp1[e.p1.a%NN]++;
		tmp1[e.p2.a%NN]++;
		NP_tmp++;
	}

	//creating nodes
	cerr<<"Creating nodes..."<<endl;
	for(const auto &p : points_tmp){
		int i=p.a;
		int b =tmp1[i];
		NP+=b;
		n[i] = new Node(i,int(b),0,Point(p.x,p.y,0));   //creating a node
	}


	cerr<<"Setting pores and neighbors..."<<endl;
	cerr<<"Numbers of pores "<<NP<<endl;
	if(NP%2==1) cerr<<"ERROR: Problem with NP = "<<NP<<endl;
	NP = NP/2;
	if(NP_tmp != NP) cerr<<"Cos nie tak w liczeniu ilosci porow"<<endl;
	p = new Pore*[NP];
	for (int i=0;i<NP;i++) p[i] = new Pore (d0, 1, i, 2); //creating pores

	//filling neighbors and pores
	int j=0; //pore number
	for(int i=0;i<NN;i++)   n[i]->tmp=0; //counting temporal nr of neighbours
	for (const auto &e : edges) {
		Node * n1 = n[e.p1.a%NN]; Node * n2 = n[e.p2.a%NN];
		Pore * p_tmp = findPore(n1,n2);
		if(p_tmp==NULL){  //to find pore only once
			if (j>=NP) cerr<<"ERROR: Problem with pore generation: j = "<<j<<endl<<flush;
			p[j]->n[0]  = n1;
			p[j]->n[1]  = n2;
			p[j]->l = point_distance (p[j]->n[0]->xy,p[j]->n[1]->xy);
			n1->n[(int)n1->tmp] = n2;  n1->p[(int)n1->tmp++] = p[j];
			n2->n[(int)n2->tmp] = n1;  n2->p[(int)n2->tmp++] = p[j];
			j++;
			}
		else cerr<<"WARNING: Double checking of pore!!!"<<e<<endl<<e.p1.a%NN<<" "<<e.p2.a%NN<<endl;
	}
	cerr<<"Checking number of connections..."<<endl<<flush;
	for(int i=0;i<NN;i++)   if(n[i]->tmp!=n[i]->b) cerr<<"ERROR: Problem in triangulation: i = "<< i <<" tmp = "<<n[i]->tmp << "  b = "<<n[i]->b<<endl<<flush;

	if(j!=NP) cerr<<"ERROR: Problem in triangulation: j = "<< j << "  NP = "<< NP <<endl;
	cerr<<"Calculating grains."<<endl;

	//creating grains

	if(NP%3!=0) cerr<<"ERROR: Problem with number of pores !!! NP/3 = "<<NP/3.<<endl;
	NG = int(NP*2./3. + 0.9);
	g = new Grain*[NG];
	cerr<<"Numbers of grain "<<NG<<endl;
	for(int i=0;i<NG;i++) g[i] = new Grain(i,0,0,0,3,3);
	int i=0;
	for (const auto &t : triangles) if(t.p1.a >=0 && t.p2.a >=0 && t.p3.a >=0){
			if((t.p1.a < NN )||(t.p2.a < NN )||(t.p3.a < NN))
				if(findGrain_T(n[t.p1.a%NN],n[t.p2.a%NN],n[t.p3.a%NN])==NULL){
					if (i>=NG) cerr<<"ERROR: Problem with pore generation: i = "<<i<<endl<<flush;
					int a1 = (t.p1.a+10*NN)%NN; int a2 = (t.p2.a+10*NN)%NN; int a3 = (t.p3.a+10*NN)%NN;
					g[i]->n[0] = n[a1];
					g[i]->n[1] = n[a2];
					g[i]->n[2] = n[a3];
					Pore * p_tmp1 = findPore(n[a1],n[a2]);
					Pore * p_tmp2 = findPore(n[a2],n[a3]);
					Pore * p_tmp3 = findPore(n[a3],n[a1]);
					if (p_tmp1 == NULL || p_tmp2==NULL || p_tmp3==NULL){
						cerr<<"WARNING: In setting grains p_tmp = NULL"<<endl<<flush; continue;}

					g[i]->p[0] = p_tmp1;   g[i]->p[1] = p_tmp2;   g[i]->p[2] = p_tmp3;
					for(int bb=0; bb<p_tmp1->bG; bb++) if(p_tmp1->g[bb]==NULL) {p_tmp1->g[bb]=g[i]; break;}
					for(int bb=0; bb<p_tmp2->bG; bb++) if(p_tmp2->g[bb]==NULL) {p_tmp2->g[bb]=g[i]; break;}
					for(int bb=0; bb<p_tmp3->bG; bb++) if(p_tmp3->g[bb]==NULL) {p_tmp3->g[bb]=g[i]; break;}
					i++;
			}
	}
	if(i!=NG) cerr<<"ERROR: Problem in triangulation: i = "<< i << "  NG = "<< NG <<endl;
	NG = i;

	//cheking if each pore has the propoer number of grains
	for(int i=0;i<NP;i++) for(int j=0;j<p[i]->bG;j++) if(p[i]->g[j]==NULL){
		if(j==1) p[i]->g[j]=g[j-1];
		else     {cerr<<"WARNING: Problem with pore "<<i<<"; None grains!!!"<<endl;}
		cerr<<"WARNING: Problem with filling info about grains: pore = "<<*p[i]<<endl;}


	//setting initial pore length and optionally adding randomness and
	if(if_randomness_in_regular_net) add_randomness_to_regular_network(gauss_sigma_d,0);
	save_info_about_initial_node_positions ();
	//updating pore lengths
	for(int i=0; i<NP;i++) p[i]->l = point_distance(p[i]->n[0]->xy, p[i]->n[1]->xy);





	cerr<<"Setting inlet and outlet..."<<endl;
	//inlet and outlet pores
	for (int i=0; i<NN; i++)
		if(int(n[i]->xy.y)==0)     N_wi++;
		else break;
	for (int i=NN-1; i>=0; i--)
		if(int(n[i]->xy.y)==N_y-1) N_wo++;
		else break;

	wi = new Node* [N_wi];
	wo = new Node* [N_wo];

	for(int i=0;i<N_wi;i++) {wi[i] = n[i];      wi[i]->t =  1;}
	for(int i=0;i<N_wo;i++) {wo[i] = n[NN-i-1]; wo[i]->t = -1;}

	//cutting vertical boundary conditions
	for (int i=0;i<NP;i++) if( abs(p[i]->n[0]->xy.y-p[i]->n[1]->xy.y)>N_y/2)          p[i]->d = 0;
	for (int i=0;i<NP;i++) if( abs(p[i]->n[0]->t) == 1 &&  abs(p[i]->n[1]->t) == 1)   p[i]->d = 0;

	//check_network_connections();

	//Check if memory should be freed
	points.clear(); points_tmp.clear(); triangles.clear(); edges.clear();


	if(if_debugging_printing){
		for(int i=0;i<NN;i++) n[i]->tmp = i;
		for(int i=0;i<NP;i++) p[i]->tmp = i;
		description_note = "Before clearing : s = " + to_string(tot_steps);
		net_ps<<*this;}


//adding information about grains to nodes
	if(if_track_grains) add_information_about_grains_in_nodes();

	NN_max = NN; NP_max = NP; NG_max = NG;


	check_network_connections();
	if(if_clear_unused_pores){
		//Attention: This step before info of grains in nodes in generated
		//clearing nodes and pores
		cerr<<"Clearing unused pores:"<<endl;
		clear_unused_pores();
		cerr<<"Clearing unused nodes:"<<endl;
		clear_unused_nodes();
		cerr<<"Clearing unused grains:"<<endl;
		clear_unused_grains();
		cerr<<endl;
		check_network_connections();
		}

	if(if_debugging_printing){
		for(int i=0;i<NN;i++) n[i]->tmp = i;
		for(int i=0;i<NP;i++) p[i]->tmp = i;
		description_note = "Just after clearing: s = " + to_string(tot_steps);
		net_ps<<*this;}


	cerr<<"Random triangles network has been created (using triangulation)."<<endl;

	return;
}


/**
* This function generates the 3D cubic network
*
* @param N network width
* @param M network length
* @param O network depth
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::createCubicNetwork(int N, int M, int O){


	cerr<<"Creating cubical network."<<endl;

	if (random_seed==-1) srand (time(NULL));
	else                 srand (random_seed);
	cerr<<"Random seed: "<<random_seed<<endl;

	NN = N*M*O;
	NP = 3*N*M*O;
	NG = N*M*O;

	b_max = 6;

	n = new Node*[NN];	p = new Pore*[NP];	g = new Grain*[NG];

	//random distribution
	default_random_engine gen;
	normal_distribution<double> N_dist(d0,d0/1); //distribution of initial diameters

	//pore and node generation
	for(int i=0;i<NN;i++) n[i] 	= new Node(6,i);
	for(int i=0;i<NP;i++) p[i] 	= new Pore(d0,1,i,0);

	//pore and node connections
	for(int k=0;k<O;k++) for(int j=0;j<M;j++) for(int i=0;i<N;i++){
		Node* nn = n[M*N*k + N*j + i];

		nn->n[0] = n[M*N*k           + N*j           + (i+1)  %N];   nn->p[0] = p[(M*N*k           + N*j           +  i        )*3 + 0];
		nn->n[1] = n[M*N*k           + N*j           + (i-1+N)%N];   nn->p[1] = p[(M*N*k           + N*j           + (i-1+N)%N )*3 + 0];
		nn->n[2] = n[M*N*k           + N*((j+1  )%M) +  i       ];   nn->p[2] = p[(M*N*k           + N*j           +  i        )*3 + 1];
		nn->n[3] = n[M*N*k           + N*((j-1+M)%M) +  i       ];   nn->p[3] = p[(M*N*k           + N*((j-1+M)%M) +  i        )*3 + 1];
		nn->n[4] = n[M*N*((k+1  )%O) + N*j           +  i       ];   nn->p[4] = p[(M*N*k           + N*j           +  i        )*3 + 2];
		nn->n[5] = n[M*N*((k-1+O)%O) + N*j           +  i       ];   nn->p[5] = p[(M*N*((k-1+O)%O) + N*j           +  i        )*3 + 2];

		p[(M*N*k + N*j + i)*3+0]->n[0] = nn; 	p[(M*N*k + N*j + i)*3+0]->n[1] = nn->n[0];
		p[(M*N*k + N*j + i)*3+1]->n[0] = nn; 	p[(M*N*k + N*j + i)*3+1]->n[1] = nn->n[2];
		p[(M*N*k + N*j + i)*3+2]->n[0] = nn; 	p[(M*N*k + N*j + i)*3+2]->n[1] = nn->n[4];

	}

	//grain generation
	for(int i=0;i<N;i++) for(int j=0;j<M;j++) for(int k=0;k<O;k++){
		Node * nn = n[M*N*k + N*j + i];
		Node * n1 = nn->n[0]->n[2]->n[4];
		Node * n_list [] = {nn,nn->n[0],nn->n[2],nn->n[4],n1,n1->n[1],n1->n[3],n1->n[5]};
		Pore * p_list [] = {nn->p[0],nn->p[2],nn->p[4],\
				          n1->p[1],n1->p[3],n1->p[5],\
						  nn->n[0]->p[2],nn->n[0]->p[4],nn->n[2]->p[4],nn->n[2]->p[0],nn->n[4]->p[2],nn->n[4]->p[0]};
		g[M*N*k + N*j + i] = new Grain (M*N*k + N*j + i, 12, 8, n_list, p_list);
	}


	if(if_verbose){  //for debugging
		cerr<<endl<<"List of  Pores:"<<endl;
		for(int i=0;i<NP;i++) cerr<<*p[i]<<endl;
		cerr<<endl<<"List of Nodes:"<<endl;
		for(int i=0;i<NN;i++) cerr<<*n[i]<<endl;
		cerr<<endl<<"List of Grains:"<<endl;
		for(int i=0;i<NG;i++) cerr<<*g[i]<<endl;
	}

	//adding nodes position (important for visualization)
	if(if_verbose) cerr<<"Adding info about node positions."<<endl;
	for(int i=0;i<N;i++) for(int j=0;j<M;j++) for(int k=0;k<O;k++) {  n[M*N*k + N*j + i]->xy = Point(i,j,k);}


	//inlet and outlet
	if(if_verbose) cerr<<"Initializing input and output..."<<endl;

	wi = new Node* [N*O];
	wo = new Node* [N*O];
	N_wi=N*O; N_wo=N*O;


	for(int i=0;i<N;i++) for(int k=0;k<O;k++){
		int ii = k*N+i;
		wi[ii] = n[k*N*M + N*0    +i];	wi[ii]->t = 1;	wi[ii]->tmp = - wi[ii]->tmp-1;
		wo[ii] = n[k*N*M + N*(M-1)+i];	wo[ii]->t =-1;	wo[ii]->tmp = - wo[ii]->tmp-1;
	}
	for(int i=0;i<N_wi;i++) for(int k=0;k<wi[i]->b;k++) if(wi[i]->p[k]->n[0]->tmp<0 && wi[i]->p[k]->n[1]->tmp<0) wi[i]->p[k]->d=0;
	for(int i=0;i<N_wo;i++) for(int k=0;k<wo[i]->b;k++) if(wo[i]->p[k]->n[0]->tmp<0 && wo[i]->p[k]->n[1]->tmp<0) wo[i]->p[k]->d=0;

	//setting initial pore length and optionally adding randomness and
	if(if_randomness_in_regular_net) add_randomness_to_regular_network(gauss_sigma_d,max_rand_shift_xy);
	save_info_about_initial_node_positions ();
	//updating pore lengths
	for(int i=0; i<NP;i++) p[i]->l = point_distance(p[i]->n[0]->xy, p[i]->n[1]->xy);




	for(int i=0;i<NN;i++) if(n[i]->a != i) cerr<<"Warning: problem with nodes names: " <<i<<"  "<<n[i]->a<<endl;
	for(int i=0;i<NP;i++) if(p[i]->a != i) cerr<<"Warning: problem with pores names: " <<i<<"  "<<p[i]->a<<endl;
	for(int i=0;i<NG;i++) if(g[i]->a != i) cerr<<"Warning: problem with grains names: "<<i<<"  "<<g[i]->a<<endl;


//adding information about grains to nodes
	if(if_track_grains) add_information_about_grains_in_nodes();

//adding information about grains to nodes
	if(if_track_grains) add_information_about_grains_in_pores();

//optional barrier (resignation form periodic boundary conditions)
	if  (!if_periodic_bc){
		//for(int k=0;k<O;k++) for(int j=0;j<M;j++)  for(int b=0;b<3;b++) p[(M*N*k + N*j + 0)*3+0]->d = 0;
		//for(int i=0;i<N;i++) for(int j=0;j<M;j++)  for(int b=0;b<3;b++) p[(M*N*0 + N*j + i)*3+0]->d = 0;
		for(int i=0; i<NP; i++) if(p[i]->n[0]->xy - p[i]->n[1]->xy > 5*l0) p[i]->d=0;
	}

	NN_max = NN; NP_max = NP; NG_max = NG;


	//saving network topology
	export_topology_file ("topology_befor_cleaning.out");

	check_network_connections();
	if(if_clear_unused_pores){
		//Attention: This step before info of grains in nodes in generated
		//clearing nodes and pores
		cerr<<"Clearing unused pores:"<<endl;
		clear_unused_pores();
		cerr<<"Clearing unused nodes:"<<endl;
		clear_unused_nodes();
		cerr<<"Clearing unused grains:"<<endl;
		clear_unused_grains();
		cerr<<endl;
		check_network_connections();
		}

	cerr<<"Cubic network has been created."<<endl;

}





/**
* This function creates a 2D network of square lattice.
*
*
* @author Rishabh
* @date 25/09/2019
*/


void Network::createSquareNetwork(int N, int M){


	cerr<<"Creating Square network..."<<endl<<flush;

	if (random_seed==-1) srand (time(NULL));
	else                 srand (random_seed);
	cerr<<"Random seed: "<<random_seed<<endl;

	NN = N*M;
	NP = 2*N*M;
	NG = N*M;
	b_max = 4;

	n = new Node*[NN];	p = new Pore*[NP]; g = new Grain*[NG];;



	//pore and node generation
	double tmp_do  = 0;
	for(int i=0;i<NN;i++) n[i] 	= new Node (4,i);
	for(int i=0;i<NP;i++) p[i] 	= new Pore (d0,l0,i,2);




	//pore and node connections
	for(int j=0;j<M;j++) for(int i=0;i<N;i++){
		Node* nn = n[N*j + i];

		nn->n[0] = n[N*j           + (i+1)  %N];   nn->p[0] = p[(N*j           +  i        )*2 + 0];
		nn->n[1] = n[N*j           + (i-1+N)%N];   nn->p[1] = p[(N*j           + (i-1+N)%N )*2 + 0];
		nn->n[2] = n[N*((j+1  )%M) +  i       ];   nn->p[2] = p[(N*j           +  i        )*2 + 1];
		nn->n[3] = n[N*((j-1+M)%M) +  i       ];   nn->p[3] = p[(N*((j-1+M)%M) +  i        )*2 + 1];


		p[(N*j + i)*2+0]->n[0] = nn; 	p[(N*j + i)*2+0]->n[1] = nn->n[0];
		p[(N*j + i)*2+1]->n[0] = nn; 	p[(N*j + i)*2+1]->n[1] = nn->n[2];
	}

	//grain generation
	for(int i=0;i<N;i++) for(int j=0;j<M;j++){
		Node* nn = n[N*j+i];
		g[N*j+i] = new Grain(N*j+i,nn,nn->n[0],nn->n[0]->n[2],nn->n[2]);
	}



	//connections between pores and grains
	if(if_verbose) cerr<<"Adding info about connections between pores and grains."<<endl;
	for(int i=0;i<NG;i++) for(int k=0;k<4;k++){
		if(g[i]->p[k]==NULL) {cerr<<"ERROR: Problem with generation hex network, g["<<i<<"] -> p["<<k<<"] = NULL"<<endl; continue;};
		if(g[i]->p[k]->g[0]==g[i] || g[i]->p[k]->g[1]==g[i])  continue;
		if(g[i]->p[k]->g[0]==NULL) { g[i]->p[k]->g[0]=g[i];   continue; }
		if(g[i]->p[k]->g[1]==NULL) { g[i]->p[k]->g[1]=g[i];   continue; }
		cerr<<"Problem with matching nodes and grains."<<endl;
	}


	if(if_verbose){  //for debugging
		cerr<<endl<<"List of  Pores:"<<endl;
		for(int i=0;i<NP;i++) cerr<<*p[i]<<endl;
		cerr<<endl<<"List of Nodes:"<<endl;
		for(int i=0;i<NN;i++) cerr<<*n[i]<<endl;
	}

	//adding nodes position (important for visualization)
	if(if_verbose) cerr<<"Adding info about node positions."<<endl;
	for(int i=0;i<N;i++) for(int j=0;j<M;j++) {  n[N*j + i]->xy = Point(i,j);}


	//inlet and outlet
	if(if_verbose) cerr<<"Initializing input and output..."<<endl;

	wi = new Node* [N];
	wo = new Node* [N];
	N_wi = N; N_wo = N;


	for(int i=0;i<N;i++) {
		wi[i] = n[N*0    +i];	wi[i]->t = 1;	wi[i]->tmp = - wi[i]->tmp-1;
		wo[i] = n[N*(M-1)+i];	wo[i]->t =-1;	wo[i]->tmp = - wo[i]->tmp-1;
	}

	//setting initial pore length and optionally adding randomness and
	if(if_randomness_in_regular_net) add_randomness_to_regular_network(gauss_sigma_d,max_rand_shift_xy);
	save_info_about_initial_node_positions ();
	//updating pore lengths
	for(int i=0; i<NP;i++) p[i]->l = point_distance(p[i]->n[0]->xy, p[i]->n[1]->xy);




	for(int i=0;i<N_wi;i++) for(int k=0;k<wi[i]->b;k++) if(wi[i]->p[k]->n[0]->tmp<0 && wi[i]->p[k]->n[1]->tmp<0) wi[i]->p[k]->d=0;
	for(int i=0;i<N_wo;i++) for(int k=0;k<wo[i]->b;k++) if(wo[i]->p[k]->n[0]->tmp<0 && wo[i]->p[k]->n[1]->tmp<0) wo[i]->p[k]->d=0;



	for(int i=0;i<NN;i++) if(n[i]->a != i) cerr<<"Warning: problem with nodes names:  " <<i<<"  "<<n[i]->a<<endl;
	for(int i=0;i<NP;i++) if(p[i]->a != i) cerr<<"Warning: problem with pores names:  " <<i<<"  "<<p[i]->a<<endl;
	for(int i=0;i<NG;i++) if(g[i]->a != i) cerr<<"Warning: problem with grains names: " <<i<<"  "<<g[i]->a<<endl;



//adding information about grains to nodes
	if(if_track_grains) add_information_about_grains_in_nodes();


//optional barrier (resignation form periodic boundary conditions)
	if  (!if_periodic_bc){
		for(int i=0; i<NP; i++) if(p[i]->n[0]->xy - p[i]->n[1]->xy > 5*l0) p[i]->d=0;
	}


	//saving network topology
	export_topology_file ("topology_befor_cleaning.out");

	check_network_connections();
	if(if_clear_unused_pores){
		//clearing nodes and pores
		cerr<<"Clearing unused pores:"<<endl;
		clear_unused_pores();
		cerr<<"Clearing unused nodes:"<<endl;
		clear_unused_nodes();
		cerr<<endl;
		check_network_connections();
		}


	cerr<<"Square network has been created."<<endl;

}


/**
* This function generates the 2D diamondal network.
* It is a square network but rotated.
* It is similar to hexagonal one but missing horizontal pores.
*
* @param N network width
* @param M network length
*
* @author Agnieszka Budek
* @date 13/03/2020
*/
void Network::createDiamondNetwork(int N, int M){

	if(M%2==1) {cerr<<"ERROR: N_y must be even in hexagonal network!"<<endl; exit(123);}

	cerr<<"Creating diamond network."<<endl;

	if (random_seed==-1) srand (time(NULL));
	else                 srand (random_seed);
	cerr<<"Random seed: "<<random_seed<<endl;

	NN = N*M;
	NP = 2*N*M;
	NG = N*M;

	b_max = 4;

	n = new Node*[NN];	p = new Pore*[NP];	g = new Grain*[NG];

	//random distribution
	default_random_engine gen;
	normal_distribution<double> N_dist(d0,d0/1.); //distribution of initial diameters

	//pore and node generation
	for(int j=0;j<M;j++) for(int i=0;i<N;i++) 						n[N*j+i] 		= new Node(4,N*j+i);
	for(int j=0;j<M;j++) for(int i=0;i<N;i++) for(int k=0;k<2;k++) 	p[(N*j+i)*2+k] 	= new Pore(d0,l0,(N*j+i)*2+k);


	//pore and node connections
	for(int j=0;j<M;j++) for(int i=0;i<N;i++){
		Node* nn = n[N*j+i];

		nn->n[0] = n[N*((j+1)%M) 	+ (i-j%2+N)%N]; 	nn->p[0] = p[(N*j+i)*2+0];
		nn->n[1] = n[N*((j+1)%M) 	+ (i-j%2+1)%N]; 	nn->p[1] = p[(N*j+i)*2+1];
		nn->n[2] = n[N*((j-1+M)%M)  + (i-j%2+1)%N];		nn->p[2] = p[(N*((j-1+M)%M) + (i-j%2+1)%N)*2+0];
		nn->n[3] = n[N*((j-1+M)%M)  + (i-j%2+N)%N]; 	nn->p[3] = p[(N*((j-1+M)%M) + (i-j%2+N)%N)*2+1];


		p[(N*j+i)*2+0]->n[0] = nn; 	p[(N*j+i)*2+0]->n[1] = nn->n[0];
		p[(N*j+i)*2+1]->n[0] = nn; 	p[(N*j+i)*2+1]->n[1] = nn->n[1];
	}


	//grain generation
	for(int i=0;i<N;i++) for(int j=0;j<M;j++){
		Node* nn = n[N*j+i];
		g[N*j+i] = new Grain(N*j+i,nn,nn->n[0],nn->n[0]->n[1],nn->n[1]);
	}


	//connections between pores and grains
	if(if_verbose) cerr<<"Adding info about connections between pores and grains."<<endl;
	for(int i=0;i<NG;i++) for(int k=0;k<4;k++){
		if(g[i]->p[k]==NULL) {cerr<<"ERROR: Problem with generation hex network, g["<<i<<"] -> p["<<k<<"] = NULL"<<endl; continue;};
		if(g[i]->p[k]->g[0]==g[i] || g[i]->p[k]->g[1]==g[i])  continue;
		if(g[i]->p[k]->g[0]==NULL) { g[i]->p[k]->g[0]=g[i];   continue; }
		if(g[i]->p[k]->g[1]==NULL) { g[i]->p[k]->g[1]=g[i];   continue; }
		cerr<<"Problem with matching nodes and grains."<<endl;
	}

	if(if_verbose){  //for debugging
		cerr<<endl<<"List of  Pores:"<<endl;
		for(int i=0;i<NP;i++) cerr<<*p[i]<<endl;
		cerr<<endl<<"List of Nodes:"<<endl;
		for(int i=0;i<NN;i++) cerr<<*n[i]<<endl;
		cerr<<endl<<"List of Grains:"<<endl;
		for(int i=0;i<NG;i++) cerr<<*g[i]<<endl;
	}

	//adding nodes position (important for visualization)
	if(if_verbose) cerr<<"Adding info about node positions."<<endl;
	for(int i=0;i<N;i++) for(int j=0;j<M;j++){  n[N*j+i]->xy = Point(i-j%2*0.5,j*sqrt(3)/2);}


	//inlet and outlet
	if(if_verbose) cerr<<"Initializing input and output..."<<endl;

	wi = new Node* [N];
	wo = new Node* [N];
	N_wi=N; N_wo=N;

	for(int i=0;i<N;i++) {
		wi[i] = n[i];		wi[i]->t = 1;	wi[i]->tmp = - wi[i]->tmp-1;
		wo[i] = n[N*M-1-i];	wo[i]->t =-1;	wo[i]->tmp = - wo[i]->tmp-1;
	}
	for(int i=0;i<N;i++) for(int k=0;k<wi[i]->b;k++) if(wi[i]->p[k]->n[0]->tmp<0 && wi[i]->p[k]->n[1]->tmp<0) wi[i]->p[k]->d=0;
	for(int i=0;i<N;i++) for(int k=0;k<wo[i]->b;k++) if(wo[i]->p[k]->n[0]->tmp<0 && wo[i]->p[k]->n[1]->tmp<0) wo[i]->p[k]->d=0;

	//setting initial pore length and optionally adding randomness and
	if(if_randomness_in_regular_net) add_randomness_to_regular_network(gauss_sigma_d,max_rand_shift_xy);
	save_info_about_initial_node_positions ();
	//updating pore lengths
	for(int i=0; i<NP;i++) p[i]->l = point_distance(p[i]->n[0]->xy, p[i]->n[1]->xy);



	//optional barrier (resignation form periodic boundary conditions)
	if  (!if_periodic_bc)  for(int i=0; i<NP; i++) if(p[i]->n[0]->xy - p[i]->n[1]->xy > 5*l0) p[i]->d=0; //old version, leaving unnecessary grains


	for(int i=0;i<NN;i++) if(n[i]->a != i) cerr<<"Warning: problem with nodes names:  " <<i<<"  "<<n[i]->a<<endl;
	for(int i=0;i<NP;i++) if(p[i]->a != i) cerr<<"Warning: problem with pores names:  " <<i<<"  "<<p[i]->a<<endl;
	for(int i=0;i<NG;i++) if(g[i]->a != i) cerr<<"Warning: problem with grains names: " <<i<<"  "<<g[i]->a<<endl;



//adding information about grains to nodes
	if(if_track_grains) add_information_about_grains_in_nodes();


	//saving network topology
	export_topology_file ("topology_befor_cleaning.out");
	NN_max = NN; NP_max = NP; NG_max = NG;


	check_network_connections();
	if(if_clear_unused_pores){
		//Attention: This step before info of grains in nodes in generated
		//clearing nodes and pores
		cerr<<"Clearing unused pores:"<<endl;
		clear_unused_pores();
		cerr<<"Clearing unused nodes:"<<endl;
		clear_unused_nodes();
		cerr<<"Clearing unused grains:"<<endl;
		clear_unused_grains();
		cerr<<endl;
		check_network_connections();
		}

	cerr<<"Hexagonal network has been created."<<endl;

}



/**
* This function adds information about grains into pores. Important if we want to track grains.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::add_information_about_grains_in_pores(){

	int *nrG = new int [NP];
	for(int i=0;i<NP;i++) nrG[i] = 0;
	for(int i=0;i<NG;i++) for(int j=0;j<g[i]->bP;j++) nrG[int(g[i]->p[j]->a)]++;
	for(int i=0;i<NP;i++) p[i]->g = new Grain *[nrG[i]];
	for(int i=0;i<NG;i++) for(int j=0;j<g[i]->bP;j++) g[i]->p[j]->add_grain(g[i]);

	delete [] nrG;

}

/**
* This function adds information about grains into nodes. Important if we want to track grains.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::add_information_about_grains_in_nodes(){

	int *nrG = new int [NN];
	for(int i=0;i<NN;i++) nrG[i] = 0;
	for(int i=0;i<NG;i++) for(int j=0;j<g[i]->bN;j++) nrG[int(g[i]->n[j]->a)]++;
	for(int i=0;i<NN;i++) n[i]->g = new Grain *[nrG[i]];
	for(int i=0;i<NG;i++) for(int j=0;j<g[i]->bN;j++) g[i]->n[j]->add_new_Grain(g[i]);

	delete [] nrG;

}
