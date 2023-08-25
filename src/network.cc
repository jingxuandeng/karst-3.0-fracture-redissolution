#include "network.h"


//functions conected with creation of the network

Network::Network (string input_file_name) {

	cerr<<endl<<"Network is being created."<<endl<<endl;

//default values

	N_x = 10;		//size of regular network
	N_y = 10;		//size of regular network

	P_in   = N_y-1;   //pressure at the inlet (must by positive)
	P_out  = 0;	      //pressure at the outlet, always should be set to zero
	Q_tot  = 2*N_x;   //total flow through the system (if == 0 the constant pressure is kept)
	Va_tot = 0;       //total volume of dissolving species
	Ve_tot = 0;       //total volume of precipitating species
	Vx_tot = 0;       //total amount of non reacting species
	Vx_perc= 0;       // percentage of non reacting species in the system


	//dimenssionless parameters describing evolution of the system

	l0    = 1;		 //initial characteristic pore length (should be always equal to one!!!!)
	d0    = 0.1*l0;	 //initial characteristic pore diameter
	Da    = 1;		 //effective Damkohler number for dissolution
	Da2   = 1;		 //effective Damkohler number for precipitation
	G1    = 1;		 //DaPe for dissolution
	G2    = 1;		 //DaPe for precipitation
	Pe1   = 1;		 //Peclet number for dissolution (D along pore)
	Pe2   = 1;		 //Peclet number for precipitation (D along pore)
	gamma = 1;	     //ratio of acid capacity numbers between dissolution and precipitation (gamma1/gamma2)
	kappa = 1;	     //ratio of Da_1/Da_2 = ratio of reaction rates
	theta = 1;	     //ratio of G_1/G_2
	d_min = d0/100.; //minimal possible pore diameter (important in precipitation)
	l_min = l0*1e-10;//minimal possible pore length (must be >0 for numerical reasons)



	//physical parameters -> should be set after choosing dimenssionless one
	
	k1	= 10e-8;		//reaction rate for dissolution
	k2	= 1;			//reaction rate for precipitation
	D1	= 0;			//diffusion coefficient for dissolution
	D2	= 1;			//diffusion coefficient for precipitation
	DD1 =-1;			//transversal diffusion coefficient for dissolution
	DD2 = 1;			//transversal diffusion coefficient for precipitation
	Sh  = 4;			//Sherwood number for pipe
	gamma_1	= 1;		//capacity number for dissolution   (c_sol = 1 by default)
	gamma_2 = 1;		//capacity number for precipitation (c_sol = 1 by default)
	Cb_0 	= 1;		//acid inlet concentration 
	Cc_0	= 0;		//precipitating species inlet concentration
	mu_0    = M_PI*pow(d0,4)/(128*l0);		//viscosity  always set to M_PI*pow(d0,4)/(128*l0)
	dt_unit = 2*k1 * gamma_1/d0;            //(in dimensionless units [2 k1 * gamma_1/d0])

	//evolution parameters
	T_max       = 10;     	  //maximal number of time steps
	tot_steps   = 0;          //total nr of steps in simulation
	tot_time    = 0;          //total time in simulation
	dt          = 1;      	  //time step (in dimensionless units [2 k1 * gamma_1/d0])
	d_d_max     = 0.1;	      //maximal change of pore diameter in one step (in %; if obtained then dt = 2/3 dt)
	d_d_min     = 0.01;	      //minimal change of pore diameter in one step (in %; if not obtained then dt = 1.2dt)
	d_d_dis     = 4;          //minimal dissolution at the outlet to finish simulation (unit d0)
	u_min       = 0.001;      ///< minimal pressure drop to finish simulation.
	d_V_min     = 0.01;       //(only for precipitation)  minimal change of pore volume in one step (if not obtained the dt = 1.2dt)
	d_V_max     = 0.1;        //(only for precipitation)  maximal change of pore volume in one step (if obtained the dt = 2/3dt)
	set_new_dt  =-1;          //0-don't change dt; 1 - increase dt by factor 1.2; -1 - decrease dt by factor 0.75
	d_max_for_u = 1000;       //maximal diameter that consume pressure, for d>d_max_for_u * d0 delta u in pore is zero

	//merging parameters
	type_of_merging = "merge_empty_grains"; //type of merging: "none", "merge_empty_grains", not implemented yet: "merge_pores"


	//printing pictures
	L_out               = 1;			          //printing scale
	pages_tot           = 100;		          //total nr of pages in the pictures (should be recalculated in future)
	pages_saved         = 0;	                  //number of printed pages
	description_note    = "Default network.";   //title note in ps/pdf files
	printing_mode       = "debugging";            //printing network style
	s_save_data         = 50;                    //how often save txt and ps files (later will be automatized)
	print_diss_factor   = 4;					   //definition of dissolution pattern for printing; only pores with d>d0*print_diss_factor are printed
	initial_xy          = NULL;                   //initial position of nodes: for printing in grains style
	pattern_anal_factor = 3;

//CONTROL PARAMETERS
	//type of network
	type_of_topology                   = "from_file";    //three options: "hexagonal", "from_file", "triangularization"
	in_topology_file_name              = "net_0.out";    //file name with input topology of the network
	in_topology_file_name_g            = "net_g_0.out";  //file name with input topology of the network
	in_pore_size_file_name             = "pores_0.out";  //file name with input pore sizes
	if_randomness_in_regular_net       = true;  	     //if true randomness is added to hexagonal network (working for hexagonal net)
	if_clear_unused_pores              = true;           //if true unused pores and nodes and grains are deleted
	if_track_grains                    = true;           //if true grains and their volume are tracked (important for merging and precipitation)
	if_periodic_bc					   = true;	         //if periodic boundary condition
	random_seed                        = -1;
	gauss_sigma_d                      = 0;    			 //if randomness is on this give information about width of the initial diameter distribution (log normal used here)
	max_rand_shift_xy                  = 1;       	     //if randomness is on this give information about max shift in positions



	//dynamics
	if_leapfrog                          = false;        //if true frog leap instead of Euler algorithm is used in evolution (not implemented yet)
	if_full_dissolution                  = true;         //if true evolution stops when system is fully dissolved
	if_system_dissolved                  = false;        //check if system is dissolved (fulfilling condition given by d_d_diss)
	if_adaptive_dt                       = true; 	     //adapting dt according to d_d_max and d_d_min;
	if_recalculate_physical_parameters   = true;         //if true recalculate physical parameters according to dimensionless one
	if_smarter_calculation_of_pressure   = true;         //if true pressure and flow is calculate in two steps
	if_precipitation				     = false;		 //if true apart form dissolution the precipitation in on
	if_dynamical_length				     = true;		 //if true length of pore is changing according to dissolution and precipitation
	if_streamtube_mixing                 = false;        // if true the stream-tube mixing is perform while calculation the concentration (works only for dissolution now)

	//output
	if_save_ps            = true;     //if true ps pictures are created
	if_save_txt           = true;     //if true data about network (nodes, pores and grains) is saved in text file
	if_save_table         = true;     //if true save tables with diameters, flow and concentration
	if_save_topology      = true; 	  //if true topology is saved in each save_all
	if_verbose            = false;    //if true verbose version for debugging
	if_debugging_printing = false;	  //if true debugging printing is done after each calculation

	//addition inlet cut
	inlet_cut_factor = 1;      //factor of an inlet cut (in a cut: d = d*factor)
	inlet_cut_w = 0;           //width of an inlet cut
	inlet_cut_l = 0;		   //length of aGn inlet cut


//Reading from the config. file

	conf_in.open (input_file_name,	ios_base::in);
	if(conf_in.is_open() == false) 	cerr<<"Problem with reading file "<<input_file_name<<endl;
	else 							read_setup_file(conf_in);

//input/output files
	net_ps    .open("net.ps",	      ios_base::out | ios_base::trunc );
	pores_out .open("pores.out",	  ios_base::out | ios_base::trunc );
	grains_out.open("grains.out",	  ios_base::out | ios_base::trunc );
	nodes_out .open("nodes.out",	  ios_base::out | ios_base::trunc );
	net_out   .open("net.out",   	  ios_base::out | ios_base::trunc );
	net_g_out .open("net_g.out",      ios_base::out | ios_base::trunc );

	time_evolution_out     .open("time_evolution.out"     ,ios_base::out | ios_base::trunc );
	pattern_analysis_out   .open("pattern_analysis.out"   ,ios_base::out | ios_base::trunc );
	child_distribution_out .open("child_distribution.out" ,ios_base::out | ios_base::trunc );
	fork_distribution_out  .open("fork_distribution.out"  ,ios_base::out | ios_base::trunc );
	cluster_size_out       .open("cluster_size.out"       ,ios_base::out | ios_base::trunc );

	if (if_save_table){
		diameters_out     .open("d.out",	      ios_base::out | ios_base::trunc );
		flow_out          .open("q.out",	      ios_base::out | ios_base::trunc );
		pressure_out      .open("u.out",	      ios_base::out | ios_base::trunc );
		concentration_out .open("c.out",	      ios_base::out | ios_base::trunc );
		concentration2_out.open("c2.out",	      ios_base::out | ios_base::trunc );
		VA_out			  .open("VA.out",	      ios_base::out | ios_base::trunc );
		VE_out			  .open("VE.out",	      ios_base::out | ios_base::trunc );
		VX_out			  .open("VX.out",	      ios_base::out | ios_base::trunc );
		lengths_out		  .open("l.out",	      ios_base::out | ios_base::trunc );
	}



	cerr<<endl<<"Initialization of the network."<<endl<<endl;

//creating network of nodes
	//srand(0); //provisional setting of random seed
	if     (type_of_topology == "hexagonal")     	createHexagonalNetwork 			(N_x,N_y);
	else if(type_of_topology == "cubic")         	createCubicNetwork     			(N_x,N_y,N_x);
	else if(type_of_topology == "diamond")          createDiamondNetwork     		(N_x,N_y);
	else if(type_of_topology == "triangulation")  	createRandomTrianglesNetwork 	(N_x,N_y);
	else if(type_of_topology == "square") 			createSquareNetwork				(N_x,N_y);
	else if(type_of_topology == "from_file"){
		import_topology_from_file  (in_topology_file_name);
		import_pore_size_from_file (in_pore_size_file_name);
		if(if_track_grains) {
			import_grains_from_file    (in_topology_file_name_g);
			add_information_about_grains_in_nodes();}
		else                 NG = 0;
		if(if_clear_unused_pores){
			cerr<<"Clearing unused pores:"<<endl;
			clear_unused_pores();
			cerr<<"Clearing unused nodes:"<<endl;
			clear_unused_nodes();
			cerr<<"Clearing unused grains:"<<endl;
			clear_unused_grains();
			cerr<<endl;
			check_network_connections();
		}
	}
	else{
		cerr<<"ERROR: wrong value of type_of_topology. "
				<<"Possible values: \n \"hexagonal\" \n \"diamond\" \n \"square\" \n \"cubic\""
						"\"from_file\" \n \"triangulation\"."<<endl;
		exit(1);}

//Check node/pore sequence (is n[s] connected through p[s])
	for(int i=0;i<NN;i++)
		for(int s=0;s<n[i]->b;s++) if(findPore(n[i],n[i]->n[s])!= n[i]->p[s]) {cerr<<"WARNING: Problem with pore/node sequence."<<endl; exit(123);}


//check if to track grains, if not NG=0
	if (!if_track_grains) NG = 0;

//additional options
	//inlets cuts
	if(inlet_cut_factor!=1) create_an_inlet_cut (inlet_cut_w, inlet_cut_l, inlet_cut_factor);
	create_an_initial_pattern();

//calculating initial Va volume for grains
	if(if_track_grains){
		cerr<<"Calculating initial grain volume..."<<endl;
		if(type_of_topology != "from_file") for(int i=0;i<NG;i++) g[i]->calculate_initial_volume(this);
		calculate_initial_total_Va();
		calculate_initial_total_Ve();}


//updating pore lengths
	if(if_dynamical_length && type_of_topology != "from_file") for(int i=0; i<NP;i++) p[i]->calculate_actual_length(this);

//calculating physical parameters based on dimensionless one
	calculate_initial_d0_and_l0 ();
	calculate_initial_mean_flow();
	if(if_recalculate_physical_parameters) {
		recalculate_k1 ();
		recalculate_DD1();
		//Add precipitation reaction parameters: k2 i DD2
	}


//saving network topology
	export_topology_file ("");
	if(if_track_grains) export_topology_file_with_grains ("");

//printing network before dissolving
	description_note = "At the beginning of the simulation.";
	for(int i=0;i<NN;i++) n[i]->tmp = n[i]->a;
	for(int i=0;i<NP;i++) p[i]->tmp = p[i]->a;
	for(int i=0;i<NG;i++) g[i]->tmp = g[i]->a;
	if(if_save_ps)   net_ps<< *this;
	if(if_save_txt)  print_net_txt();




	//Optionally I print info about the system at the beginning of the simulation
	if(if_verbose){
		for(int i=0; i<NN; i++) cerr<<"Node:   "<<*n[i]<<endl;
		for(int i=0; i<NP; i++) cerr<<"Pore:   "<<*p[i]<<endl;
		for(int i=0; i<NG; i++) cerr<<"Grains: "<<*g[i]<<endl;
	}


	cerr<<endl<<endl<<"Network has been generated."<<endl<<"NN = "<<NN<<"   NP = "<<NP<<"   NG = "<<NG<<endl<<endl<<endl;
}


Network:: ~Network (){

	cerr<<"Network is being deleted."<<endl;

	for(int i=0;i<NN;i++) delete n[i];
	for(int i=0;i<NP;i++) delete p[i];
	for(int i=0;i<NG;i++) delete g[i];

	delete[] n; n=NULL;
	delete[] p; p=NULL;
	if(if_track_grains) delete[] g; g=NULL;

//file closing
	
	net_ps    .close();
	net_out   .close();
	pores_out .close();
	grains_out.close();
	nodes_out .close();
	net_out   .close();
	net_g_out .close();

	time_evolution_out     .close();
	pattern_analysis_out   .close();
	child_distribution_out .close();
	fork_distribution_out  .close();

	if (if_save_table){
		diameters_out      .close();
		flow_out           .close();
		pressure_out       .close();
		concentration_out  .close();
		concentration2_out .close();
		VA_out			   .close();
		VE_out 			   .close();
		VX_out 			   .close();
		lengths_out		   .close();
	}

}

