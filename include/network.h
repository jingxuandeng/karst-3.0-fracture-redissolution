/**
 * @file
 * @author  Agnieszka Budek
 * @date 25/09/2019
 * @class Network
 *
 * @section DESCRIPTION
 *
 * Class Network represents the porous material as a network of pores, nodes and grains.
 * Flow of reagents through the network results in dissolution and optionally precipitation reactions.
 * As a results of those reactions the pores properties will change.
 *
 * Class Network consist of list of all pores, nodes and grains and all physical parameters governing the reaction.
 * Additionally in the class Network there are control parameters for the entire simulation.
 *
 */


#ifndef NETWORK_H
#define NETWORK_H 

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <list>

#include "pore.h"
#include "grain.h"
#include "node.h"
#include "constants.h"
#include "printing.h"

using namespace std;

class Punkt;

/// main class, represents the simulated network and includes all physical properties and  simulation parameters
class Network  
{
	public:
		
		Node** n;			///< node list
		Pore** p;			///< pore list
		Grain** g;			///< grain list
		int b_max;			///< maximal number of neighbors
			
		Node** wi;			///< inlet list
		Node** wo;			///< outlet list

		int NN;				///< number of nodes
		int NP;				///< number of pores
		int NG;				///< number of grains
		
		int NN_max;			///< maximal number of nodes (before merging)
		int NP_max;			///< maximal number of pores (before merging)
		int NG_max;			///< maximal number of grains (before merging)

		int N_x;			///< size of regular network (for hexagonal and triangularization but should be set (the average value) for printing ps)
		int N_y;			///< size of regular network (for hexagonal and triangularization but should be set (the average value) for printing ps)
		int N_wi;			///< number of inlet nodes
		int N_wo;			///< number of outlet nodes

		// physical parameters
		double P_in;		///< pressure at the inlet (by default equal to N_y, don't change when Q_tot!=0) (must by positive)
		double P_out;		///< pressure at the outlet (by default equal to N_y, don't change when Q_tot!=0) (should be set to zero)
		double Q_tot;		///< total flow through the system
		double Va_tot;		///< total amount of volume of dissolving species
		double Ve_tot;		///< total amount of volume of precipitating species
		double Vx_tot;      ///< total amount of non reacting
		double Vx_perc;     ///< percentage of non reacting species

		double q_in_0;      ///< initial mean flow through pores (by definition initial flow through inlet pores)

		double k1;			///< reaction rate of dissolution
		double k2;			///< reaction rate of precipitation
		double D1;			///< diffusion coefficient for dissolution (not implemented yet, we assume that flow is faster then diffusion along pore)
		double D2;			///< diffusion coefficient for precipitation (not implemented yet)
		double DD1;			///< transversal diffusion coefficient for dissolution
		double DD2;			///< transversal diffusion coefficient for precipitation
		double Sh;			///< Sherwood number for pipe
		double gamma_1;		///< capacity number for dissolution (important only for unit of dt, generally should be set to 1)
		double gamma_2;		///< capacity number for precipitation
		double Cb_0;		///< acid inlet concentration
		double Cc_0;		///< inlet concentration of precipitating species
		double mu_0;        ///< viscosity  (by default set to M_PI*pow(d0,4)/(128*l0))
		double dt_unit;     ///< unit of time step (in dimensionless units [2 k1 * gamma_1/d0] or in diffusion limited case 2 DD1*Sh*gamma/d0^2)

		// dimenssionless parameters
		double d0;		///< initial characteristic pore diameter [natural unit is l0]
		double l0;		///< initial characteristic pore length (l0 is by defoult set to one)
		double Da;		///< effective Damkohler number for dissolution
		double Da2;		///< effective Damkohler number for precipitation
		double G1;		///< DaPe for dissolution
		double G2;		///< DaPe for precipitation
		double Pe1;		///< Peclet number for dissolution (D along pore) (not used now)
		double Pe2;		///< Peclet number for precipitation (D along pore)  (not used now)
		double gamma;	///< ratio of acid capacity numbers between dissolution and precipitation
		double kappa;	///< ratio of Da_2/Da_1 of reaction rates (dissolution vs precipitation)
		double theta;   ///< ratio of G2/G1 (dissolution vs precipitation)
		double d_min;	///< minimal possible pore diameter (important for precipitation)
		double l_min;   ///< minimal pore length (for numerical reason)


		// evolution parameters
		int    T_max;	  	///< maximal number of time steps in simulation
		int    tot_steps; 	///< total nr of steps in simulation
		double tot_time;  	///< total time in simulation
		double dt;			///< time step (in dimensionless units [2 k1 * gamma_1/d0] or in diffusion limited case 2 DD1*Sh*gamma/d0^2)
		double d_d_max;		///< maximal change of pore volume in one step (if obtained the dt = 2/3dt)
		double d_max_for_u; ///< maximal diameter that consume pressure, for d>d_max_for_u*d0 delta u in pore is zero
		double d_d_min;		///< minimal change of pore diameter in one step (if not obtained the dt = 1.2dt)
		double d_V_max;     ///< (only for precipitation)  maximal change of pore volume in one step (if obtained the dt = 2/3dt)
		double d_V_min;     ///< (only for precipitation)  minimal change of pore volume in one step (if not obtained the dt = 1.2dt)
		double d_d_dis; 	///< minimal dissolution at the outlet to finish simulation (unit d0)
		double u_min;       ///< minimal pressure drop to finish simulation.
		int    set_new_dt;  ///< 0-don't change dt; 1 - increase dt by factor 1.2; -1 - decrease dt by factor 0.75
		
		// merging parameters
		string type_of_merging; ///< type of merging: "none", "merge_empty_grains", "empty_grains_to_pores", "merge_pores"

		// input and output files
		ofstream_ps 	net_ps;
		ofstream_txt    pores_out, nodes_out, grains_out, net_out, net_g_out;
		ofstream_txt    time_evolution_out, pattern_analysis_out, child_distribution_out, fork_distribution_out, cluster_size_out;
		ofstream_txt    diameters_out, flow_out, pressure_out, concentration_out, concentration2_out, VA_out, VE_out, VX_out, lengths_out;
	   	ifstream 	    conf_in, net_in, net_g_in, pores_in, grains_in;
	

	//CONTROL PARAMETERS
	   	// type of network
	   	string type_of_topology;        ///< three options: "hexagonal", "form_file", "triangularization"
	   	string in_topology_file_name;   ///< file with initial topology of the network
	   	string in_pore_size_file_name;  ///< file with initial pore size of the network
	   	string in_topology_file_name_g; ///< file with initial topology of grains
		bool if_clear_unused_pores;     ///< if true unused pores and nodes and grains are deleted
	    bool if_track_grains;           ///<  if true grains are tracked (important for merging and precipitation)
		bool if_periodic_bc;            ///< if false - network without boundary conditions is created, works for hexagonal, cubic or triangulation (default true)
	    double random_seed;             ///< random seed for generating random networks
	   	bool if_randomness_in_regular_net;  ///< if true randomness is added to hexagonal network (working for hexagonal net)
	   	double gauss_sigma_d;           ///< if randomness is on this give information about width of the initial diameter distribution
	   	double max_rand_shift_xy;       ///< if randomness is on this give information about max shift in positions

		// dynamics
		bool if_leapfrog;               ///< if true frog leap instead of Euler algorithm is used in evolution (not implemented yet)
		bool if_full_dissolution;       ///< if true evolution stops when system is fully dissolved
		bool if_system_dissolved;       ///< check if system is dissolved (fulfilling condition given by d_d_diss)
		bool if_adaptive_dt;			///< adapting dt according to d_d_max and d_d_min;
		bool if_recalculate_physical_parameters;    ///< if true, recalculate physical parameters according to dimensionless one, by now must be true
		bool if_smarter_calculation_of_pressure;    ///< if true pressure and flow is calculate in two steps
		bool if_precipitation;					 	///< if true apart form dissolution the precipitation in on
		bool if_dynamical_length;					///< if true length of pore is changing according to dissolution and precipitation
		bool if_streamtube_mixing;					///< if true and we have square lattice  stream-tube mixing is performed while calculating the species B concentration


		// output
		bool if_save_ps;                          ///< if true ps pictures are created
		bool if_save_txt;                         ///< if true data about network (nodes, pores, grains) in text file
		bool if_save_table; 					  ///< if true tables with diameters, flow and concentration are being saved
		bool if_save_topology; 					  ///< if true topology is saved in each save_all
		bool if_verbose;                          ///< if true verbose version for debugging
		bool if_debugging_printing;				  ///< if true debugging printing is done after each calculation

		// printing parameters
		// Point * xy;			///< list of nodes positions, filled when printing_ps is ON
		double L_out;				///< scale factor for printing pores
		int pages_tot;				///< total nr of pages in the pictures
		int pages_saved;			///< number of printed pages
		string printing_mode;       ///<  printing network style
		string description_note;    ///< note to be printed in the title of the page, description of current state of simulation
		double s_save_data;         ///< how often save txt and ps files if <0 save after breaking every n-th level (row)
		double print_diss_factor;   ///< definition of dissolution pattern for printing; only pores with d>d0*print_diss_factor are printed
		Point ** initial_xy;        ///< initial position of nodes: for printing in grains style
		double pattern_anal_factor; ///< threshold for pattern analysis (pores with d> d0*diss_pattern_factor are treated as pattern when looking for a spanning tree)

		// addition inlet cut
		double inlet_cut_factor;       ///< factor of an inlet cut (in a cut: d = d*factor)
		int    inlet_cut_w;            ///< width of an inlet cut
		int    inlet_cut_l;			   ///< length of an inlet cut

	public:
		Network  (string input_file_name);
		~Network();

	private:	
	
	public:

// evolution of the system
		void evolution(long int T = 0);				///< make whole evolution
		void do_one_leapfrog_step();				///< do one leap-frog step, WARNING: not implemented yet
		void do_one_euler_step();					///< do one euler step, simplest numerical method for integrating differential equation of the system evolution
		void calculate_pressures();					///< calculate pressure field for the whose network
		void calculate_flows();						///< calculate flow field for the whole network
		void calculate_concentrations();			///< calculate concentration profile for species b
		void calculate_concentrations_c();			///< calculate concentration profile for species c
		void calculate_concentrations_streamtube_mixing(); 		     ///< new fancy mixing method where particles prefer to go straight through the crossing
 		void dissolve();											 ///< change the pore sizes due to the dissolution
		void dissolve_and_precipitate();							 ///< change the pore sizes due to both dissolution and precipitation
		void calculate_pressures_and_flows_smarter(double d_max);    ///< alternative way of calculating pressure and flow field, using d_max: pores larger then d_max do not consume pressure drop but can consume acid :)
		void calculate_pressures_for_small_d(double d_max);			 ///< part of an alternative way of calculating pressure and flow field, WARNING: must be tested more carefully
		void calculate_pressures_for_large_d(double d_max);			 ///< part of an alternative way of calculating pressure and flow field, WARNING: must be tested more carefully
		void calculate_flows_for_large_d(double d_max);				 ///< part of an alternative way of calculating pressure and flow field, WARNING: must be tested more carefully


		void adapt_dt();												///< part of adapting the time step, for the dissolution to be not too slow not too fast
		void recalculate_flows_to_keep_Q_tot(string type_of_nodes);     ///< when total flow through the system is kept the flow field must be rescaled in each time step
		void set_adaptive_dt(double dd, double dV);                     ///< the time dt will be adapted in each time step according to the speed of dissolution and precipitation
		void check_if_dissolved();										///< checks if the system is dissolved, if yes the simulation stops

		//Properties of particular pores
		double outlet_c_b       (Pore* p);    	///< returns the outlet concentration of species b in the pore p as a function of c0_b
		double outlet_c_c_1     (Pore* p);		///< returns the outlet concentration of species c in the pore p as a function of c0_b
		double outlet_c_c_2     (Pore* p);		///< returns the outlet concentration of species c in the pore p as a function of c0_c
		double k_eff_in_pore    (Pore* p);		///< returns the value of k_eff in the pore p
		double k_eff_2_in_pore  (Pore* p);		///< returns the value of k_eff_2 (see precipitation) in the pore p
		double G_in_pore        (Pore* p);		///< returns the value of G in the pore p
		double G2_in_pore       (Pore* p);		///< returns the value of g2 (see precipitation) in the pore p
		double Da_eff_in_pore   (Pore* p);		///< returns the value of Da_eff in the pore p
		double Da_eff_2_in_pore (Pore* p);		///< returns the value of Da_eff_2 (see precipitation) in the pore p

		void check_diss_pattern(double threshold);		///< check the position of the dissolution pattern - pores larger that threshold x d0
		bool check_diss_front  (double threshold, int row);  ///< check if the dissolution pattern -reaches particular row of pores

// calculating initial properties of the system
		void calculate_initial_mean_flow();		///< calculate initial mean flow through the system, important for setting Da_eff and G correctly
		void calculate_initial_d0_and_l0 ();    ///< calculate initial mean pore length and diameter, important for setting Da_eff and G correctly
		void calculate_initial_total_Va();		///< calculate initial amount of species A, important for mass balance
		void calculate_initial_total_Ve();		///< calculate initial amount of species E, important for mass balance
		void recalculate_k1 ();					///< recalculate reaction rate, k, has no impact on the simulation, just for potential curiosity
		void recalculate_DD1();					///< recalculate diffusion coefficient, has no impact on the simulation, just for potential curiosity
		void recalculate_k2 ();					///< recalculate reaction rate for second reaction, has no impact on the simulation, just for potential curiosity
		void recalculate_DD2();					///< recalculate diffusion coefficient, has no impact on the simulation, just for potential curiosity

// creating network
		void createHexagonalNetwork(int N, int M);	            		///< creating hexagonal network of nodes/pores
		void createDiamondNetwork(int N, int M);						///< creating diamondal network of nodes/pores (its like hexagonal but without horizontal pores)
		void read_network_of_nodes(string net_in_name);					///< reading network of nodes/pores from a file
		void createCubicNetwork(int N, int M, int O);					///<  creating 3D cubic network
		void createRandomTrianglesNetwork(int N1, int M2);			    ///< creating random network of pores;
		void createSquareNetwork(int N, int M);							///<  creating 2D square network (check if works)
		void create_an_inlet_cut(int cut_w, int cut_l, double factor);	///< adding one central cut
		void create_an_initial_pattern();                               ///< creating initial pattern for debugging purpose
		void add_randomness_to_regular_network(double d_sigma, double max_nodes_shift);   ///< Add random node position shift or random initial diameters for hexagonal network
		void add_information_about_grains_in_nodes();					///< if tracking grains is necessary we need info of grain in connected nodes
		void add_information_about_grains_in_pores();					///< if tracking grains is necessary we need info of grain in connected pores

// maintaining the network
		Pore*  findPore   (Node* n1, Node* n2);
		Grain* findGrain_T(Node* n1, Node* n2, Node *n3);
		Grain* findGrain  (Node* n1, Node* n2, Node *n3);
		void clear_unused_nodes();								///< delete nodes  unconnected to the rest of the network
		void clear_unused_pores();								///< delete pores  unconnected to the rest of the network
		void clear_unused_grains();								///< delete grains unconnected to the rest of the network
		void clear_unconeccted_pores();                    		///< deletes parts of the network that are not connected to the inlet/outlet
		void move_pore_to_the_end  (int i, int NP_tmp);			///< delete pore form the network but still preserve info about it in the end of pore list
		void move_grain_to_the_end (int i, int NG_tmp);			///< delete node form the network but still preserve info about it in the end of node list
		void move_node_to_the_end  (int i, int NN_tmp);			///< delete grain form the network but still preserve info about it in the end of grain list
		double node_distance (Node  *n1, Node   *n2);      		///< return distance between two nodes
		double point_distance(Point p1, Point  p2);	       		///< return distance between two point taking into account boundary conditions
		void save_info_about_initial_node_positions ();    		///< for printing in grain style only


// network merging
		void do_merging();										///< do standard merging
		void merge_empty_grains();								///< merge empty (Va + Ve <=0) grains
		void merge_nodes(Node *n1, Node *n2);					///< merge two nodes, this is a pert of merging one grain function
		void merge_one_grain(Grain *g);							///<  the main function responsible of merging grains
		void merge_one_pore_grain (Grain *g);					///< merging special type of grain, the one-pore one
		void debugging_for_merging(string text="In merging");	///< Extra debugging for merging


// I/O
		void print_net_txt   ();			///< print info about network in text mode
		void print_tables_txt();			///< print info about network in tables (works better for regular networks)
		void save_all_data                      (bool if_save_now = false);		    ///< save info about network, the format depends of if_dave_... flags
		void read_setup_file                    (ifstream& input_file);   		///< read setup file with initial parameters  for the whole simulation
		void export_topology_file               (string out_file_name = "");  	///< export basic info about topology: list of nodes with their neighboring nodes and pores
		void export_topology_file_with_grains   (string out_file_name = "");  	///< export extra info about topology: list of grains and their neighboring nodes and pores
		void import_topology_from_file          (string in_file_name);   		///< import basic info about topology: list of nodes with their neighboring nodes and pores
		void import_pore_size_from_file         (string in_file_name);			///< import info about pore sizes (l and d)
		void import_grains_from_file            (string in_file_name);			///< import extra info about grains: list of grains and their neighboring nodes and pores plus info avout Va and Ve

// verification
		void check_flow_balance();		  	///< Checks if the flow through the system is conserve (if not the pressure field is not calculated precisely enough)
		void check_acid_balance();			///< Checks if the amount of all species is conserved, if not some problem with merging occurs
		void check_precipitating_balance(); ///< Checks if the amount of all species is conserved, when both dissolution and precipitation occures
		void check_network_connections();   ///< Checks consistency of network connections, if false the initial network is badly defined or there is a problem with merging
		void check_GMash_connections();		///< FOr checking GMash only
		void print_network_for_debugging (string text = "Network",string type_n = "name", string type_p="name", string type_g = "name"); ///< Print all information about network in large pdf (WARNING: not to be used for large networks and large no of time steps)

// results analysis
		void analyze_diss_pattern(bool if_anal_pattern = false);
		void write_time_step_properties();
		void write_pattern_properties();
		void find_the_largest_tree(double threshold);
		void find_minimal_spanning_tree(double threshold, Node* root);
		int* find_child_distribution ();
		void find_forks();
		double find_reverse_forks();
		double find_cluster_size (int n);
		int* find_fork_distribution();
		double find_minimal_tree_length();
		double distance_to_root(Node *);  ///< Returns distance to the root of the tree, if node does't belongs to the tree -1 is returned.

//other output generation (for Rishabh)
		void write_vtk_data();
		void write_point_data(string file_name);
		void write_line_data(string file_name);
		void write_cell_data(string file_name);
		void write_diameter(string file_name);
		void write_flow_rate(string file_name);
		void write_concentration(string file_name);


};


#endif


