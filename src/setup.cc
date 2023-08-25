#include "network.h"

#include <string>  
#include <iostream>
#include <sstream>


/**
* This function reads setup information form setup file.
* @param fp_setup setup file name
* @author Agnieszka Budek
* @date 25/09/2019
*/
void::Network::read_setup_file(ifstream& fp_setup){

	int i=0;
	string s;

	while(getline(fp_setup,s)){
		i++;
		
		if(s.size()<3) 	continue; //cerr<<"Line "<<i<<" is empty."<<endl; continue;}
		if(s[0]=='#' ) 	continue;

		stringstream line;
		line << s;
		
		string name,value, e;
		value = "";
		
		line >> name >> e >> value;
		
		if(e!="=" || value=="") {cerr<<"Problem with parsing line nr "<<i<<endl; continue;}



		if(name == "N_x"){
			N_x = stod(value);
			Q_tot = 2*N_x;
			cerr<< "Setting N_x = "<<N_x<<endl;
			cerr<< "Additionally setting Q_tot = "<<Q_tot<<endl;
		}

		else if(name == "N_y"){
			N_y = stod(value);
			P_in  = N_y-1;
			cerr<< "Setting N_y = " <<N_y <<endl;
			cerr<< "Additionally setting P_in = "<<P_in<<endl;
		}

		else if(name == "N_wi"){
			N_wi = stod(value);
			cerr<< "Setting N_wi = "<<N_wi<<endl;}

		else if(name == "N_wo"){
			N_wo = stod(value);
			cerr<< "Setting N_wo = "<<N_wo<<endl;}

		else if(name == "P_in"){
			P_in = stod(value);
			cerr<< "Setting P_in = "<<P_in<<endl;}

		else if(name == "P_out"){
			P_out = stod(value);
			cerr<< "Setting P_out = "<<P_out<<endl;}

		else if(name == "Q_tot"){
			Q_tot = stod(value);
			cerr<< "Setting Q_tot = "<<Q_tot<<endl;}

		else if(name == "q_in_0"){
			q_in_0 = stod(value);
			cerr<< "Setting q_in_0 = "<<q_in_0<<endl;}

		else if(name == "k1"){
			k1 = stod(value);
			cerr<< "Setting k1 = "<<k1<<endl;}

		else if(name == "k2"){
			k2 = stod(value);
			cerr<< "Setting k2 = "<<k2<<endl;}

		else if(name == "D1"){
			D1 = stod(value);
			cerr<< "Setting D1 = "<<D1<<endl;}

		else if(name == "D2"){
			D2 = stod(value);
			cerr<< "Setting D2 = "<<D2<<endl;}

		else if(name == "DD1"){
			DD1 = stod(value);
			cerr<< "Setting DD1 = "<<DD1<<endl;}

		else if(name == "DD2"){
			DD2 = stod(value);
			cerr<< "Setting DD2 = "<<DD2<<endl;}

		else if(name == "Sh"){
			Sh = stod(value);
			cerr<< "Setting Sh = "<<Sh<<endl;}

		else if(name == "gamma_1"){
			gamma_1 = stod(value);
			cerr<< "Setting gamma_1 = "<<gamma_1<<endl;}

		else if(name == "gamma_2"){
			gamma_2 = stod(value);
			cerr<< "Setting gamma_2 = "<<gamma_2<<endl;}

		else if(name == "Cb_0"){
			Cb_0 = stod(value);
			cerr<< "Setting Cb_0 = "<<Cb_0<<endl;}

		else if(name == "Cc_0"){
			Cc_0 = stod(value);
			cerr<< "Setting Cc_0 = "<<Cc_0<<endl;}

		else if(name == "mu_0"){
			mu_0 = stod(value);
			cerr<< "Setting mu_0 = "<<mu_0<<endl;}

		else if(name == "d0"){
			d0 = stod(value);
			cerr<< "Setting d0 = "<<d0<<endl;}

		else if(name == "l0"){
			l0 = stod(value);
			cerr<< "Setting l0 = "<<l0<<endl;}

		else if(name == "Da"){
			Da = stod(value);
			cerr<< "Setting Da = "<<Da<<endl;}

		else if(name == "Da2"){
			Da2 = stod(value);
			cerr<< "Setting Da2 = "<<Da2<<endl;}

		else if(name == "G1"){
			G1 = stod(value);
			cerr<< "Setting G1 = "<<G1<<endl;}

		else if(name == "G2"){
			G2 = stod(value);
			cerr<< "Setting G2 = "<<G2<<endl;}

		else if(name == "Pe1"){
			Pe1 = stod(value);
			cerr<< "Setting Pe1 = "<<Pe1<<endl;}

		else if(name == "Pe2"){
			Pe2 = stod(value);
			cerr<< "Setting Pe2 = "<<Pe2<<endl;}

		else if(name == "gamma"){
			gamma = stod(value);
			cerr<< "Setting gamma = "<<gamma<<endl;}

		else if(name == "kappa"){
			kappa = stod(value);
			cerr<< "Setting kappa = "<<kappa<<endl;}

		else if(name == "theta"){
			theta = stod(value);
			cerr<< "Setting theta = "<<theta<<endl;}

		else if(name == "d_min"){
			d_min = stod(value);
			cerr<< "Setting d_min = "<<d_min<<endl;}

		else if(name == "l_min"){
			l_min = stod(value);
			cerr<< "Setting l_min = "<<l_min<<endl;}

		else if(name == "Cb_0"){
			Cb_0 = stod(value);
			cerr<< "Setting Cb_0 = "<<Cb_0<<endl;}

		else if(name == "Cc_0"){
			Cc_0 = stod(value);
			cerr<< "Setting Cc_0 = "<<Cc_0<<endl;}

		else if(name == "Vx_perc"){
			Vx_perc = stod(value);
			cerr<< "Setting Vx_perc = "<<Vx_perc<<endl;}

		else if(name == "T_max"){
			T_max = stod(value);
			cerr<< "Setting T_max = "<<T_max<<endl;}

		else if(name == "tot_steps"){
			tot_steps = stod(value);
			cerr<< "Setting tot_steps = "<<tot_steps<<endl;}

		else if(name == "tot_time"){
			tot_time = stod(value);
			cerr<< "Setting tot_time = "<<tot_time<<endl;}

		else if(name == "dt"){
			dt = stod(value);
			cerr<< "Setting dt = "<<dt<<endl;}

		else if(name == "d_max_for_u"){
			d_max_for_u = stod(value);
			cerr<< "Setting d_max_for_u = "<<d_max_for_u<<endl;}

		else if(name == "d_d_max"){
			d_d_max = stod(value);
			cerr<< "Setting d_d_max = "<<d_d_max<<endl;}

		else if(name == "u_min"){
			u_min = stod(value);
			cerr<< "Setting u_min = "<<u_min<<endl;}

		else if(name == "d_d_min"){
			d_d_min = stod(value);
			cerr<< "Setting d_d_min = "<<d_d_min<<endl;}

		else if(name == "d_V_min"){
			d_V_min = stod(value);
			cerr<< "Setting d_V_min = "<<d_V_min<<endl;}

		else if(name == "d_V_max"){
			d_V_max = stod(value);
			cerr<< "Setting d_V_max = "<<d_V_max<<endl;}

		else if(name == "d_d_dis"){
			d_d_dis = stod(value);
			cerr<< "Setting d_d_dis = "<<d_d_dis<<endl;}

		else if(name == "set_new_dt"){
			set_new_dt = stod(value);
			cerr<< "Setting set_new_dt = "<<set_new_dt<<endl;}

		else if(name == "type_of_topology"){
			type_of_topology =   value;
			cerr<< "Setting type_of_topology = "<<type_of_topology<<endl;}

		else if(name == "type_of_merging"){
			type_of_merging =   value;
			cerr<< "Setting type_of_merging = "<<type_of_merging<<endl;}

		else if(name == "in_topology_file_name"){
			in_topology_file_name = value;
			cerr<< "Setting in_topology_file_name = "<<in_topology_file_name<<endl;}

		else if(name == "in_topology_file_name_g"){
			in_topology_file_name_g = value;
			cerr<< "Setting in_topology_file_name_g = "<<in_topology_file_name_g<<endl;}

		else if(name == "in_pore_size_file_name"){
			in_pore_size_file_name = value;
			cerr<< "Setting in_pore_size_file_name = "<<in_pore_size_file_name<<endl;}

		else if(name == "if_periodic_bc"){
			if      (value == "true" )   if_periodic_bc  = true;
			else if (value == "false")   if_periodic_bc  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_periodic_bc. Set true or false."<<endl;
			cerr<< "Setting if_periodic_bc = "<<if_periodic_bc<<endl;}


		else if(name == "if_randomness_in_regular_net"){
			if      (value == "true" )   if_randomness_in_regular_net  = true;
			else if (value == "false")   if_randomness_in_regular_net  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_randomness_in_regular_net. Set true or false."<<endl;
			cerr<< "Setting if_randomness_in_regular_net = "<<if_randomness_in_regular_net<<endl;}

		else if(name == "if_leapfrog"){
			if      (value == "true" )   if_leapfrog  = true;
			else if (value == "false")   if_leapfrog  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_leapfrog. Set true or false."<<endl;
			cerr<< "Setting if_leapfrog = "<<if_leapfrog<<endl;}

		else if(name == "if_full_dissolution"){
			if      (value == "true" )   if_full_dissolution  = true;
			else if (value == "false")   if_full_dissolution  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_full_dissolution. Set true or false."<<endl;
			cerr<< "Setting if_full_dissolution = "<<if_full_dissolution<<endl;}

		else if(name == "if_system_dissolved"){
			if      (value == "true" )   if_system_dissolved  = true;
			else if (value == "false")   if_system_dissolved  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_system_dissolved. Set true or false."<<endl;
			cerr<< "Setting if_system_dissolved = "<<if_system_dissolved<<endl;}

		else if(name == "if_adaptive_dt"){
			if      (value == "true" )   if_adaptive_dt  = true;
			else if (value == "false")   if_adaptive_dt  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_adaptive_dt. Set true or false."<<endl;
			cerr<< "Setting if_adaptive_dt = "<<if_adaptive_dt<<endl;}

		else if(name == "if_recalculate_physical_parameters"){
			if      (value == "true" )   if_recalculate_physical_parameters  = true;
			else if (value == "false")   if_recalculate_physical_parameters  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_recalculate_physical_parameters. Set true or false."<<endl;
			cerr<< "Setting if_recalculate_physical_parameters = "<<if_recalculate_physical_parameters<<endl;}

		else if(name == "if_smarter_calculation_of_pressure"){
			if      (value == "true" )   if_smarter_calculation_of_pressure  = true;
			else if (value == "false")   if_smarter_calculation_of_pressure  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_smarter_calculation_of_pressure. Set true or false."<<endl;
			cerr<< "Setting if_smarter_calculation_of_pressure = "<<if_smarter_calculation_of_pressure<<endl;}

		else if(name == "if_precipitation"){
			if      (value == "true" )   if_precipitation  = true;
			else if (value == "false")   if_precipitation  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_precipitation. Set true or false."<<endl;
			cerr<< "Setting if_precipitation = "<<if_precipitation<<endl;}

		else if(name == "if_dynamical_length"){
			if      (value == "true" )   if_dynamical_length  = true;
			else if (value == "false")   if_dynamical_length  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_dynamical_length. Set true or false."<<endl;
			cerr<< "Setting if_dynamical_length = "<<if_dynamical_length<<endl;}

		else if(name == "if_streamtube_mixing"){
			if      (value == "true" )   if_streamtube_mixing  = true;
			else if (value == "false")   if_streamtube_mixing  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_streamtube_mixing. Set true or false."<<endl;
			cerr<< "Setting if_streamtube_mixing = "<<if_streamtube_mixing<<endl;}

		else if(name == "if_save_ps"){
			if      (value == "true" )   if_save_ps  = true;
			else if (value == "false")   if_save_ps  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_save_ps. Set true or false."<<endl;
			cerr<< "Setting if_save_ps = "<<if_save_ps<<endl;}

		else if(name == "if_save_txt"){
			if      (value == "true" )   if_save_txt  = true;
			else if (value == "false")   if_save_txt  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_save_txt. Set true or false."<<endl;
			cerr<< "Setting if_save_txt = "<<if_save_txt<<endl;}

		else if(name == "if_save_topology"){
			if      (value == "true" )   if_save_topology  = true;
			else if (value == "false")   if_save_topology  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_save_topology. Set true or false."<<endl;
			cerr<< "Setting if_save_topology = "<<if_save_topology<<endl;}


		else if(name == "if_save_table"){
			if      (value == "true" )   if_save_table  = true;
			else if (value == "false")   if_save_table = false;
			else                     cerr<<"WARNING: Wrong value of variable if_save_table. Set true or false."<<endl;
			cerr<< "Setting if_save_table = "<<if_save_table<<endl;}


		else if(name == "if_verbose"){
			if      (value == "true" )   if_verbose  = true;
			else if (value == "false")   if_verbose  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_verbose. Set true or false."<<endl;
			cerr<< "Setting if_verbose = "<<if_verbose<<endl;}

		else if(name == "if_debugging_printing"){
			if      (value == "true" )   if_debugging_printing  = true;
			else if (value == "false")   if_debugging_printing  = false;
			else                     cerr<<"WARNING: Wrong value of variable if_debugging_printing. Set true or false."<<endl;
			cerr<< "Setting if_debugging_printing = "<<if_debugging_printing<<endl;}

		else if(name == "if_clear_unused_pores"){
					if      (value == "true" )   if_clear_unused_pores  = true;
					else if (value == "false")   if_clear_unused_pores  = false;
					else                     cerr<<"WARNING: Wrong value of variable if_clear_unused_pores. Set true or false."<<endl;
					cerr<< "Setting if_clear_unused_pores = "<<if_clear_unused_pores<<endl;}

		else if(name == "if_track_grains"){
					if      (value == "true" )   if_track_grains  = true;
					else if (value == "false")   if_track_grains  = false;
					else                     cerr<<"WARNING: Wrong value of variable if_track_grains. Set true or false."<<endl;
					cerr<< "Setting if_track_grains = "<<if_track_grains<<endl;}


		else if(name == "L_out"){
			L_out = stod(value);
			cerr<< "Setting L_out = "<<L_out<<endl;}

		else if(name == "s_save_data"){
			s_save_data = stod(value);
			cerr<< "Setting s_save_data = "<<s_save_data<<endl;}

		else if(name == "pages_tot"){
			pages_tot = stod(value);
			cerr<< "Setting pages_tot = "<<pages_tot<<endl;}

		else if(name == "pages_saved"){
			pages_saved = stod(value);
			cerr<< "Setting pages_saved = "<<pages_saved<<endl;}

		else if(name == "printing_mode"){
			printing_mode = value;
			cerr<< "Setting printing_mode = "<<printing_mode<<endl;}

		else if(name == "description_note"){
			description_note = value;
			cerr<< "Setting description_note = "<<description_note<<endl;}

		else if(name == "print_diss_factor"){
			print_diss_factor = stod(value);
			cerr<< "Setting print_diss_factor = "<<print_diss_factor<<endl;}

		else if(name == "pattern_anal_factor"){
			pattern_anal_factor = stod(value);
			cerr<< "Setting diss_pattern_factor = "<<pattern_anal_factor<<endl;}

		else if(name == "inlet_cut_factor"){
			inlet_cut_factor = stod(value);
			cerr<< "Setting inlet_cut_factor = "<<inlet_cut_factor<<endl;}

		else if(name == "inlet_cut_w"){
			inlet_cut_w = stod(value);
			cerr<< "Setting inlet_cut_w = "<<inlet_cut_w<<endl;}

		else if(name == "inlet_cut_l"){
			inlet_cut_l = stod(value);
			cerr<< "Setting inlet_cut_l = "<<inlet_cut_l<<endl;}

		else if(name == "random_seed"){
			random_seed = stod(value);
			cerr<< "Setting random_seed = "<<random_seed<<endl;}

		else if(name == "gauss_sigma_d"){
			gauss_sigma_d = stod(value);
			cerr<< "Setting gauss_sigma_d = "<<gauss_sigma_d<<endl;}

		else if(name == "max_rand_shift_xy"){
			max_rand_shift_xy = stod(value);
			cerr<< "Setting max_rand_shift_xy = "<<max_rand_shift_xy<<endl;}


		else {cerr<<"Problem with parsing name \""<<name<<"\" in line "<<i<<endl;}


		
	}
	cerr<<"I read "<<i<<" lines."<<endl;

//Checking consistency

	if(!if_track_grains){
		if(type_of_merging!="none") {
			cerr<<"WARNING: Performing merging requires grains tracking!"<<endl;
			cerr<<"WARNING: type_of_merging is set to none."<<endl;
			type_of_merging="none";
		}
		if(if_clear_unused_pores) {
			cerr<<"WARNING: Clearing unused pores requires grains tracking!"<<endl;
			cerr<<"WARNING: if_clear_unused_pores is set to false."<<endl;
			if_clear_unused_pores = false;
		}
	}

	if(if_streamtube_mixing == true) {
		if(type_of_topology!="square" && type_of_topology!="diamond") {
			cerr<<"WARNING: Stream-tube mixing works only for square topology. "<<endl;
			if_streamtube_mixing = false;
			cerr<<"if_streamtube_mixing = " <<if_streamtube_mixing<<endl;
		}
		if(if_precipitation) {
			cerr<<"WARNING: Stream-tube mixing works only for pure dissolution."<<endl;
			if_streamtube_mixing = false;
			cerr<<"if_streamtube_mixing = " <<if_streamtube_mixing<<endl;
		}

	}
	if(!if_precipitation) {d_V_min = 0; d_V_max = 0;}


//adjust additional parameters that depends on above
	mu_0  = M_PI*pow(d0,4)/(128.*l0); cerr<<"Setting mu_0 = "<<mu_0<<endl;
	P_out = 0;	     //pressure at the outlet, always should be set to zero
	NN    = N_x * N_y;
	Da2   = kappa * Da;
	G2    = theta * G1;


	cerr<<"At the beginning of the simulation:\nDa = "<<Da<<"\nG1 = "<<G1<<"\nDa2 = "<<Da2<<"\nG2 = "<<G2<<endl;

}
