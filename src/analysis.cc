#include "network.h"
#include "printing.h"


//functions related to evolution

/**
* This function analyzes dissolution pattern.
* For example: pattern length, ramification, fractal dim, etc.
* Additionally in each time step basic information about system is printed in write_time_step_properties();
* @author Agnieszka Budek
* @date 26/03/2020
*/
void Network::analyze_diss_pattern(bool if_anal_pattern){


	write_time_step_properties();   //time step properties are printed always




	//deciding either to analyze pattern or not
	//now we print pattern properties as we save other data, can be changed later
	static double Va_old = 0;

	if (s_save_data<0 and !if_anal_pattern)  if_anal_pattern = check_diss_front(print_diss_factor, pages_saved*abs(s_save_data));
	else if (s_save_data>=1 && tot_steps%int(s_save_data)==0) if_anal_pattern = true;
	else if (s_save_data>0 && s_save_data<1){ //check the volume condition!!!
		if (tot_steps==0) Va_old = 0;
		if (fabs(Va_old-tot_time*dt_unit)/(T_max*dt_unit)>s_save_data){
			if_anal_pattern = true;
			Va_old = tot_time*dt_unit;
		}
	}

	if(if_anal_pattern) write_pattern_properties  ();

}


/**
* This function analyzes dissolution pattern.
* For example: pattern length, ramification, fractal dim, etc.
* @author Agnieszka Budek
* @date 30/03/2020
*/
void Network::write_pattern_properties(){


	if(tot_steps==0){
			pattern_analysis_out   << "#"                   <<\
						setw(15)   << "step  "              <<\
						setw(15)   << "time  "              <<\
						setw(15)   << "time*dt_unit  "      <<\
						setw(15)   << "diss_pattern_l "     <<\
                        setw(15)   << "preci_pattern_l "    <<\
						setw(15)   << "horton_diss"         <<\
                        setw(15)   << "horton_preci"        <<\
						endl;
			pattern_analysis_out<<"#  ----------------------------------------------------------------------------------------------------------------"<<endl;
		}



    // 1. precipitation pattern analysis
	find_the_largest_tree(0.5,true);
    double max_horton_pre  = find_reverse_forks();
    double dist_pre        = find_minimal_tree_length();


    //finding forks distribution
    int * child_dist = find_child_distribution();
    int * fork_dist  = find_fork_distribution();

    child_distribution2_out<<setw(10)<<tot_steps;
    fork_distribution2_out <<setw(10)<<tot_steps;
    for(int i=0; i<=child_dist[0]; i++) child_distribution2_out<<setw(8)<<child_dist[i];
    for(int i=0; i<=fork_dist [0]; i++)  fork_distribution2_out<<setw(8)<< fork_dist[i];
    child_distribution2_out<<endl<<flush;
    fork_distribution2_out <<endl<<flush;

    delete [] child_dist;
    delete [] fork_dist;


    //finding cluster size (cluster growing method) (for fractal dimension approximation)
    cerr<<"Finding clusters size for diss pattern..."<<endl;
    cluster_size2_out<<setw(10)<<tot_steps;
    for(int i=3; i<2*int(pow(NN,0.25)); i++) cluster_size2_out<<setw(5)<<i<<setw(12)<<find_cluster_size(i);
    cluster_size2_out <<endl<<flush;



    //2. dissolution pattern analysis
    find_the_largest_tree(pattern_anal_factor);
	double max_horton_diss = find_reverse_forks();
	double dist_diss       = find_minimal_tree_length();


	//finding forks distribution
	child_dist = find_child_distribution();
	fork_dist  = find_fork_distribution();

	child_distribution_out<<setw(10)<<tot_steps;
	fork_distribution_out <<setw(10)<<tot_steps;
	for(int i=0; i<=child_dist[0]; i++) child_distribution_out<<setw(8)<<child_dist[i];
	for(int i=0; i<=fork_dist [0]; i++)  fork_distribution_out<<setw(8)<< fork_dist[i];
	child_distribution_out<<endl<<flush;
	fork_distribution_out <<endl<<flush;

	delete [] child_dist;
	delete [] fork_dist;

	//finding cluster size (cluster growing method) (for fractal dimension approximation)
	cerr<<"Finding clusters size for diss pattern..."<<endl;
	cluster_size_out<<setw(10)<<tot_steps;
	for(int i=3; i<2*int(pow(NN,0.25)); i++) cluster_size_out<<setw(5)<<i<<setw(12)<<find_cluster_size(i);
	cluster_size_out <<endl<<flush;



    //printing pattern properties for both diss/preci mode
    pattern_analysis_out        <<\
		setw(15) << tot_steps   <<\
		setw(15) << tot_time    <<\
		setw(15) << tot_time*dt_unit    <<\
		setw(15) << dist_diss   <<\
        setw(15) << dist_pre    <<\
		setw(15) << max_horton_diss  <<\
        setw(15) << max_horton_pre   <<\
		endl<<flush;



	if(printing_mode=="debugging" ){//for debugging only
				description_note = "At the end of write pattern properties: s = " + to_string(tot_steps);
				for (int i=0;i<NP;i++) p[i]->tmp = 666;
				printing_mode = "debugging";
				net_ps<< *this;
				//printing_mode = "grains";  //optionally one can switch to the grain printing_mode
		}
}


/**
* This function prints basic system properties:
* - step
* - total time
* - total flow
* - inlet pressure
* - total A volume
* - total E volume
* - maximal diameter at inlet
* - maximal diameter at outlet
* @author Agnieszka Budek
* @date 30/03/2020
*/
void Network::write_time_step_properties(){

	//printing first line
	if(tot_steps==0){
		time_evolution_out<<"#" <<setw(15)<<\
				"step"          <<setw(15)<<\
				"time"          <<setw(15)<<\
				"time*dt_unit"  <<setw(15)<<\
				"Q"             <<setw(15)<<\
				"P"             <<setw(15)<<\
				"VA_tot"        <<setw(15)<<\
				"VE_tot"        <<setw(15)<<\
				"VX_tot"        <<setw(15)<<\
                "plots"         <<setw(15)<<\
                "percolation"   <<setw(15)<<\
                "sim_state"     <<setw(15)<<\
				"d_in_max"      <<setw(15)<<\
				"d_out_max"     <<setw(15)<<endl;
		time_evolution_out<<"#  ----------------------------------------------------------------------------------------------------------------"<<endl;
	}

	//calculating maximal pore diameter at the outlet
	double d_out_max = 0;double d_in_max = 0;
	for(int i=0;i<N_wi;i++)
        for(int s=0;s<wi[i]->b;s++)
            if(d_in_max<wi[i]->p[s]->d)
                d_in_max=wi[i]->p[s]->d;

    for(int i=0;i<N_wo;i++)
        for(int s=0;s<wo[i]->b;s++)
            if(d_out_max<wo[i]->p[s]->d)
                d_out_max=wo[i]->p[s]->d;


    double percolation = find_percolation();



	time_evolution_out <<setw(15)<<\
			tot_steps  <<setw(15)<<\
			tot_time   <<setw(15)<<\
			tot_time*dt_unit   <<setw(15)<<\
			Q_tot_tmp  <<setw(15)<<\
			wi[0]->u   <<setw(15)<<\
			Va_tot     <<setw(15)<<\
			Ve_tot     <<setw(15)<<\
			Vx_tot     <<setw(15)<<\
            pages_saved<<setw(15)<<\
            percolation<<setw(15)<<\
            sim_state  <<setw(15)<<\
			d_in_max   <<setw(15)<<\
			d_out_max  <<endl<<flush;
}

