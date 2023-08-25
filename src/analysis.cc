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

	if (s_save_data<0)                                   if_anal_pattern = check_diss_front(print_diss_factor, pages_saved*abs(s_save_data));
	else if (s_save_data>=1 && tot_steps%int(s_save_data)==0) if_anal_pattern = true;
	else if (s_save_data>0 && s_save_data<1){ //check the volume condition!!!
		if (tot_steps==0) Va_old = Va_tot;
		if ((Va_old-Va_tot)/Va_tot>s_save_data){
			if_anal_pattern = true;
			Va_old = Va_tot;
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
			pattern_analysis_out   << "#"                 <<\
						setw(15)   << "step  "            <<\
						setw(15)   << "time  "            <<\
						setw(15)   << "pattern_length "   <<\
						setw(15)   << "horton "           <<\
						endl;
			pattern_analysis_out<<"#  ----------------------------------------------------------------------------------------------------------------"<<endl;
		}

	find_the_largest_tree(pattern_anal_factor);
	int max_horton = find_reverse_forks();
	double dist    = find_minimal_tree_length();


	//printing pattern properties
	pattern_analysis_out        <<\
		setw(15) << tot_steps   <<\
		setw(15) << tot_time    <<\
		setw(15) << dist        <<\
		setw(15) << max_horton  <<\
		endl<<flush;

	//finding forks distribution
	int * child_dist = find_child_distribution();
	int * fork_dist  = find_fork_distribution();

	child_distribution_out<<setw(10)<<tot_steps;
	fork_distribution_out <<setw(10)<<tot_steps;
	for(int i=0; i<=child_dist[0]; i++) child_distribution_out<<setw(8)<<child_dist[i];
	for(int i=0; i<=fork_dist [0]; i++)  fork_distribution_out<<setw(8)<< fork_dist[i];
	child_distribution_out<<endl<<flush;
	fork_distribution_out <<endl<<flush;

	delete [] child_dist;
	delete [] fork_dist;

	//finding cluster size (cluster growing method) (for fractal dimension approximation)
	cerr<<"Finding clusters size..."<<endl;
	cluster_size_out<<setw(10)<<tot_steps;
	for(int i=3; i<2*int(pow(NN,0.25)); i++) cluster_size_out<<setw(5)<<i<<setw(12)<<find_cluster_size(i);
	cluster_size_out <<endl<<flush;

	if(printing_mode=="debugging"){//for debugging only
				description_note = "At the end of write pattern properties: s = " + to_string(tot_steps);
				for (int i=0;i<NP;i++) p[i]->tmp = 666;
				printing_mode = "debugging";
				net_ps<< *this;
				printing_mode = "grains";
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
				"Q"             <<setw(15)<<\
				"P"             <<setw(15)<<\
				"VA_tot"        <<setw(15)<<\
				"VE_tot"        <<setw(15)<<\
				"VX_tot"        <<setw(15)<<\
				"d_in_max"      <<setw(15)<<\
				"d_out_max"     <<setw(15)<<endl;
		time_evolution_out<<"#  ----------------------------------------------------------------------------------------------------------------"<<endl;
	}

	//calculating maximal pore diameter at the outlet
	double d_out_max = 0;double d_in_max = 0;
	for(int i=0;i<N_wi;i++) for(int s=0;s<wi[i]->b;s++)
		if(d_in_max<wi[i]->p[s]->d) d_in_max=wi[i]->p[s]->d;
	for(int i=0;i<N_wo;i++) for(int s=0;s<wo[i]->b;s++)
		if(d_out_max<wo[i]->p[s]->d) d_out_max=wo[i]->p[s]->d;



	time_evolution_out <<setw(15)<<\
			tot_steps  <<setw(15)<<\
			tot_time   <<setw(15)<<\
			Q_tot      <<setw(15)<<\
			wi[0]->u   <<setw(15)<<\
			Va_tot     <<setw(15)<<\
			Ve_tot     <<setw(15)<<\
			Vx_tot     <<setw(15)<<\
			d_in_max   <<setw(15)<<\
			d_out_max  <<endl<<flush;
}

