#include "network.h"
#include "printing.h"

//functions related to evolution


/**
* This function simulates the evolution of the system.
*
* @param T maximal no of time steps to be simulated
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::evolution(double T){

	cerr<<"\n\n\nEvolution:"<<endl;

	if (T!=0) T_max = T/dt_unit;
    else      T_max = T_max/dt_unit;

	if(if_save_txt)     print_net_txt();
	save_all_data       (true);
	analyze_diss_pattern(true);

	while(tot_time < T_max){
		cerr<<endl;
		cerr << tot_steps << ". step of evolution, tot_time = " << tot_time << " (t = "<<tot_time*dt_unit<< ")" << endl;

		if(if_leapfrog)  do_one_leapfrog_step();	//do frog leap version (not implemented yet)
		else			 do_one_euler_step();		//do normal Euler version

		tot_steps++;
        tot_time+=dt;

		save_all_data();
		analyze_diss_pattern();

		if(if_system_dissolved){
			save_all_data       (true);
			analyze_diss_pattern(true);
			break;
		}


	}
	save_all_data(true);
	//analyze_diss_pattern(true);
}

/**
* This function do one euler time step.
* The order of calculations is as follows:
*
* 1. Calculate pressure field
*
* 2. Calculate flow field
*
* 3. Calculate the concentration profiles.
*
* 4. Change pore sizes according to the dissolution (and precipitation).
*
* 5. Check the flow and mass balance.
*
* 6. Do merging (optionally).
*
* This is the simplest method for integrating of differential equations describing the system dynamics.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::do_one_euler_step(){

    static bool if_percolation = false;

	if(if_smarter_calculation_of_pressure)	{calculate_pressures_and_flows_smarter(d_max_for_u);}
	else 									{calculate_pressures();           calculate_flows();}

	if(if_precipitation){

		calculate_concentrations();
		calculate_concentrations_c();
		dissolve_and_precipitate();
		}
	else{
		if(if_streamtube_mixing) 	calculate_concentrations_streamtube_mixing();
		else						calculate_concentrations();

	    dissolve();
        cerr<<"Kocham Janka."<<endl;
	}

	//checking acid and flow balance
	check_flow_balance();
	if(if_track_grains) check_acid_balance();
	if(if_precipitation&&if_track_grains) check_precipitating_balance();
	//check_GMash_connections();  //additional checking for G_Mash network

	//topology changes
    do_merging();


	if(if_adaptive_dt)      adapt_dt();           //adapt dt
	if(if_full_dissolution) check_if_dissolved(); //check if the system is dissolved
    //checking percolation
    if (!if_percolation && find_percolation()>0 && tot_steps>10){
        cerr<<"\nSystem clogged due to d_min percolation.\nSimulation finished."<<endl;
        if_percolation = true;
        if_system_dissolved=true;
        T_max = tot_time;  //after clogging just end the simulation
        return;}

}
/**
* This function will do one leap-frog step.
* Will be implemented later.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::do_one_leapfrog_step(){

//trzeba będzie też znak pradu dyfuzji lepiej liczyć w leaf frogu!!!
//bo tam biore znak z wczesniejszego kroku (mozna by to iteracyjnie, ale moim zdaniem bez sensu strata czasu)

}

