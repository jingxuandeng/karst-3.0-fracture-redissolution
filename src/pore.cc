#include "pore.h"
#include "network.h"
#include "printing.h"
#include "grain.h"
#include <tuple>
// #define _USE_MATH_DEFINES
// #include <cmath>

Pore::Pore (double dd, double ll, float name, int bb){
	d = dd; l = ll; a=name; tmp=name; q=0; x=1; bG=bb; c_in=0;
    is_active = true; is_fracture = false; l0=0;
	n[0]=NULL; n[1]=NULL;	
	if(bG>0){
		g = new Grain*[bG];
		for (int i=0;i<bG;i++) g[i]=NULL;
	}
}

double Pore::perm(Network*S){
    bool pipe_formula = (!(S->sandwich_pores and is_fracture) and (d<S->H_z or S->no_max_z));

    if(!pipe_formula and d<S->H_z)
        return M_PI*pow(d,3)/(128*S->mu_0*l);   ///< permeability of a thin fracture
    if(pipe_formula)
        return M_PI*pow(d,4)/(128*S->mu_0*l);   ///< permeability of a particular pore
    else
        return M_PI*d/(128*S->mu_0*l);          ///WARNING: the tube can not have the diameter larger than H_z, later the formula for H_z = l_Z = 1
}



/**
* This function returns the inlet concentration of species B for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_inlet_cb(){
	if (q>0)  return n[0]->cb;
	else      return n[1]->cb;
}

/**
* This function returns the outlet concentration of species B for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_outlet_cb(){
	if (q>0) return n[1]->cb;
	else     return n[0]->cb;
}

/**
* This function returns the inlet concentration of species C for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_inlet_cc(){
	if (q>0)  return n[0]->cc;
	else      return n[1]->cc;
}

/**
* This function returns the outlet concentration of species C for the pore depending on flow direction through the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::calculate_outlet_cc(){
	if (q>0) return n[1]->cc;
	else     return n[0]->cc;
}


/**
* This function returns true if the species A is still present in neighboring grains.
* This information is important when deciding if the dissolution can still occur in the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Pore::is_Va_left(){
	double Va_tot=0;
	for (int b=0;b<bG;b++) Va_tot+=g[b]->Va;
	if(Va_tot>0)   return true;
	else		   return false;
}

/**
* This function returns true if the species E is still present in neighboring grains.
* This information is important when deciding if the redissolution can still occur in the pore.
* @author Jingxuan Deng
* @date 07/07/2025
*/
bool Pore::is_Ve_left(){
    double Ve_tot=0;
    for (int b=0;b<bG;b++) Ve_tot+=g[b]->Ve;
    if(Ve_tot>0)   return true;
    else		   return false;
}

/**
* This function returns true if the species E is being generated from precipitation step.
* This information is important when deciding if the redissolution can still occur in the pore.
* @author Jingxuan Deng
* @date 04/08/2025
*/

bool Pore::is_Ve_generated(Network *S){

	// The is_Ve_generated() function returns false only when there is full conversion of A to C.
	// Therefore, before computing Ve_prec, we can check certain conditions (ccin == 0 and Va == 0)
	// to avoid unnecessary computations and save time.
	// Va == 0 means there is no dissolution of A to generate C, and thus, no precipitation occurs.
	// Precipitation also depends on the inlet concentration of species C,
	// so if cc0 == 0, it indicates no precipitation can happen.
	double cc0 = calculate_inlet_cc();
	if (cc0==0||!is_Va_left()) return false;

	double cb0 = calculate_inlet_cb();
	double f1= local_Da_eff(S);
	double f2= local_Da_eff_2(S);
	double f3= local_Da_eff_3(S);
	double g= local_G(S);

	// Compute total volume of mineral E precipitated during this step.
	double Ve_prec_tmp = (S->dt*S->d0*M_PI*d*l/(2*f1*(1+g)))*S->gamma*(cc0*(1-exp(-f1-f3))+(cb0*((f1+f3)*(1-exp(-f2))-f2*(1-exp(-f1-f3)))/(f1+f3-f2)));
	if (Ve_prec_tmp>0) return true;
	else		   return false;
}

/**
* This function returns the actual length the pore if the if_dynamical_length flag is on.
*
* In this case the length of a pore changes in time depending on dissolution and precipitation processes.
* @param S pointer to the network
* @param l_max maximal length for the network, important for periodic boundary conditions
* @param l_0 characteristic pore's length, if there is a problem with calculation max length l_0 is returned
* @author Agnieszka Budek
* @date 26/09/2019
*/
void Pore::calculate_actual_length(Network *S, double l_max, double l_0){


	l = S->point_distance(n[0]->xy, n[1]->xy);  //l0;

    if (is_fracture and S->type_of_merging == "merge_for_fracture") return;        // fracture pore has always maximal length


	if (d == 0) return;

	if(!S->if_dynamical_length || !S->if_track_grains) return;

	double factor=0; double V_max=0; double V_act=0;
	for (int s=0; s<bG; s++){
		V_max = g[s]->V0;
		V_act = g[s]->Va + g[s]->Ve + g[s]->Vx;
		if(V_max>0 && V_act>0) factor += pow(V_act/V_max,1./3);  //One may ask about exponent for semi 2D network, maybe it should be two?
		}

	if(factor>0) 	l = l*factor/bG;
	else			l = S->l_min;

	if(l<=S->l_min) l = S->l_min;

	if(l<=S->l_min) {                                                   //FIXME: OPCJA: maleńkie ziarno powinnno miec opcje wtornego osadzania
		if(S->if_verbose) cerr<<"l = l_min for Pore:"<<*this<<endl;
		for (int s=0; s<bG; s++){
			if(S->if_verbose) cerr<<"The following gains is going to be empty."<<endl\
					<<"Grain:"<<*g[s]<<endl;
			S->Va_tot -= g[s]->Va; g[s]->Va =0;
            S->Ve_tot -= g[s]->Ve; g[s]->Ve =0;


		}
	}
}

/**
* This function returns the local value of G parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_G(Network* S){

    //double d_tmp = min(1,d);  //possible feature for a fracture

	if     (S->G1==0)	return	0;						         // reaction limited case, G = 0
	else if(S->G1>0 )	return  S->G1*d/S->d0;	             // mixed case: k1 ~ DD1
	else			    return  -1;							     // diffusion limited case, convention: G<0 => G = Inf

}

/**
* This function returns the local value of G2 (used in precipitation) parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_G_2(Network* S){

    //double d_tmp = min(1,d);   //possible feature for a fracture

	if     (S->G2==0)	return	0;						         // reaction limited case, G = 0
	else if(S->G2>0 )	return  S->G2*d/S->d0;	             // mixed case: k1 ~ DD1
	else			    return  -1;							     // diffusion limited case, convention: G<0 => G = Inf

}

/**
* This function returns the local value of G3 (used in redissolution) parameter for the pore.
* @param S pointer to the network
* @author Jingxuan Deng
* @date 16/06/2025
*/
double Pore::local_G_3(Network* S){

	//double d_tmp = min(1,d);   //possible feature for a fracture

	if     (S->G3==0)	return	0;						         // reaction limited case, G = 0
	else if(S->G3>0 )	return  S->G3*d/S->d0;	             // mixed case: k1 ~ DD1
	else			    return  -1;							     // diffusion limited case, convention: G<0 => G = Inf

}

/**
* This function returns the local value of Da_eff parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_Da_eff(Network* S){

    //if(!is_Va_left()) return 0; //0.000001;
	if (q==0) return -1;
	double G = this->local_G(S);

    //formula for aperture
    if((d>S->H_z and !S->no_max_z) or (S->sandwich_pores and is_fracture)) {
        if (G > 0)  return S->Da * (1 / S->d0) * (l / S->l0) * (S->q_in_0 / fabs(q)) * ((1 + S->G1) / (1 + G));
        if (G == 0) return S->Da * (1 / S->d0) * (l / S->l0) * (S->q_in_0 / fabs(q));
    }

    // old formula for a cylinder
	if      (G>0)    return S->Da*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G1)/(1+G));
	else if (G==0)   return S->Da*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	else             return S->Da*(l/S->l0)*(S->q_in_0/fabs(q));
}

double Pore::is_there_precipitation(Network *S){

    if (S->C_eq==0 and is_Ve_left()) return  1;  //if C_eq == 0 the precipitation is reversible

    double cc0 = calculate_inlet_cc(); // 20250708: This could be negative because the definition of R=k2theta(c_eq-c) for reversible precipitation
    double cb0 = calculate_inlet_cb();
    if(!is_Va_left()) cb0 = 0;

    if(cc0 >= 0)
        return 1;       //normal behaviour with precipitation

    double f1 = local_Da_eff(S);
	double f3= local_Da_eff_3(S);
    // double dcc = cb0*(1-exp(-f1));
	double dcc = cb0*(1-exp(-f1-f3));

    if(cc0+dcc < 0)
        return 0;

	// 20250708: because cb0 must be positive, so if cc0 is negative
    double alpha = 1./f1 * log(cb0/(cb0+cc0));
    if(alpha>1 or alpha<0) {cerr<<"ERROR: is_there_precipitation has wrong value."; return 0;}
    return (1-alpha);

}
// To do 20250708
// this function should be value between 0 and 1, value inbetween means we are redissolving that not comes from the equation
/**
* This function returns a value between 0 and 1 that control the redissolution process.
* @param S pointer to the network
* @author Jingxuan Deng
* @date 08/07/2025
*/
double Pore::is_there_redissolution(Network *S){
	// cerr<<"running is_there_redissolution..."<<endl;
	double cc0 = calculate_inlet_cc();
	double cb0 = calculate_inlet_cb();
	double f1 = local_Da_eff(S);
	double f2 = local_Da_eff_2(S);
	double f3 = local_Da_eff_3_tmp(S);
	double g = local_G(S);

	// 1. Compute initial total volume of mineral E in all grains around the pore
	double Ve_0= 0;
	for (int b = 0; b < bG; ++b) {
		double Ve_b=this->g[b]->Ve;
		double Ve_diss2_max_tot=0;
		for (int pi = 0; pi < this->g[b]->bP; ++pi) {
			auto* pore_i = this->g[b]->p[pi];
			double cb0_i = pore_i ->calculate_inlet_cb();
			double f1_i = pore_i ->local_Da_eff(S);
			double f3_i = pore_i ->local_Da_eff_3_tmp(S);
			double g_i = pore_i ->local_G(S);
			double d_i= pore_i ->d;
			double l_i= pore_i ->l;

			// total maximum amount of mineral E of a grain will be redissolved.
			double k_ratio=f3_i/(f1_i+f3_i);
			double time_factor=(S->dt*S->d0*M_PI*d_i*l_i/(2*f1_i*(1+g_i)));
			double exp_factor=1-exp(-f1_i-f3_i);
			Ve_diss2_max_tot += 0.5*k_ratio*time_factor*(S->gamma)*cb0_i*exp_factor;
			// cerr<<"k_ratio= "<<k_ratio<<". time factor= "<<time_factor<<". exp_factor= "<<exp_factor<<". gamma= "<<S->gamma<<". cb0_i= "<<cb0_i<<". Ve_diss2_max_tot= "<<Ve_diss2_max_tot<<endl;
		}
		// cerr<<"Ve_diss2_max_tot_final= "<<Ve_diss2_max_tot<<endl;
		double Ve_diss2_b =0.5*(f3/(f1+f3))*(S->dt*S->d0*M_PI*d*l/(2*f1*(1+g)))*S->gamma*cb0*(1-exp(-f1-f3));
		// cerr<<"Ve_diss2_b = "<<Ve_diss2_b<<endl;

		double w_b = 0;
		if (Ve_diss2_max_tot>0) w_b = Ve_diss2_b/Ve_diss2_max_tot;
		Ve_0 +=w_b*Ve_b;

		// cerr<<"w_b= "<<w_b<<" and Ve_b= "<<Ve_b<<endl;
		// cerr<<"Ve_0= "<<Ve_0<<endl;
	}

	// // Compute total available Ve in grain
	// double Ve_tot=0;
	// for (int b=0;b<bG;b++) Ve_tot+=this->g[b]->Ve;
	// cerr<<"Ve_tot= "<<Ve_tot<<endl;


	// 2. Compute total volume of mineral E precipitated during this step
	double Ve_prec_tmp = (S->dt*S->d0*M_PI*d*l/(2*f1*(1+g)))*S->gamma*(cc0*(1-exp(-f1-f3))+(cb0*((f1+f3)*(1-exp(-f2))-f2*(1-exp(-f1-f3)))/(f1+f3-f2)));

	// cerr<<"Ve_prec_tmp = "<<Ve_prec_tmp<<endl;

	// 3. Compute total redissolution of mineral E during this step
	double Ve_diss2_tmp =(f3/(f1+f3))*(S->dt*S->d0*M_PI*d*l/(2*f1*(1+g)))*S->gamma*cb0*(1-exp(-f1-f3));

	// cerr<<"Ve_diss2_tmp = "<<Ve_diss2_tmp<<endl;

	double alpha = 0;
	if (Ve_diss2_tmp>0) alpha = (Ve_0+Ve_prec_tmp) / Ve_diss2_tmp;

	// cerr<<"alpha= "<<alpha<<endl;

	// if (alpha>=1) return 1; // normal behavior
	// else if (alpha==0) return 0;
	// else if (alpha<0) {cerr<<"ERROR: is_there_redissolution has wrong value."; return 0;}
	// else return -f1-log(1-((f1+f3)/f3)*(alpha*Ve_diss2_tmp)*(2*f1*(1+g))/(S->dt*S->d0*M_PI*d*l*S->gamma)*(1/cb0));

	double Da_factor=0;
	// if (alpha>=1) Da_factor=1; // normal behavior
	// else if (alpha==0) Da_factor= 0;
	// else if (alpha<0) {cerr<<"ERROR: is_there_redissolution has wrong value."; return 0;}
	// else {
	// 	double k_factor=(f3+f1)/f3;
	// 	double scale_t=(2*f1*(1+g))/(S->dt*S->d0*M_PI*d*l*S->gamma*cb0);
	// 	Da_factor=-f1-log(1-k_factor*(alpha*Ve_diss2_tmp)*scale_t);
	// }
	if (alpha>=1) Da_factor=1; // normal behavior
	else if (alpha<1 and alpha>=0) Da_factor=0; // no redissolution
	else  {cerr<<"ERROR: is_there_redissolution has wrong value."; return 0;}

	// cerr<<"Da_factor= "<<Da_factor<<endl;
	return  Da_factor;
}

/**
* This function returns the local value of Da_eff_2 (used in precipitation) parameter for the pore.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Pore::local_Da_eff_2(Network* S){

	if (q==0) return -1;
	double G = this->local_G_2(S);
    double Da2local = S->Da2;

    //double d_tmp = min(1,d);     //possible feature for a fracture

	// precipitation depends on pH treshold
    if(S->if_dynamic_k2){

        if(abs(S->dyn_k2_alpha)>100){
            if(calculate_inlet_cb()*_sign(S->dyn_k2_alpha)>S->dyn_k2_c0*_sign(S->dyn_k2_alpha))
                return 0; // if inlet cB_in concentration of a given pore is larger than treshold value, then 0. Right now this is implemented, dont need to worry about the else{...}
        }
        else{
            double kappa = 1./(1+pow(calculate_inlet_cb()/S->dyn_k2_c0,S->dyn_k2_alpha));
            Da2local = Da2local*kappa;
        }
    }

    Da2local = Da2local * is_there_precipitation(S); // this could be useful form for redissolution as well

    //formula for an aperture
    if((d>S->H_z and !S->no_max_z) or (S->sandwich_pores and is_fracture)){
        if      (G>0)    return Da2local*(1/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G2)/(1+G));
        else if (G==0)   return Da2local*(1/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
    }

    //old formula for a cylinder
	if      (G>0)    return Da2local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G2)/(1+G));
	else if (G==0)   return Da2local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	else             return Da2local*(l/S->l0)*(S->q_in_0/fabs(q));
}

/**
* This function returns the local value of Da_eff_3 (used in redissolution) parameter for the pore.
* @param S pointer to the network
* @author Jingxuan Deng
* @date 16/06/2025
*/
// To do, modify this
double Pore::local_Da_eff_3(Network* S){
	// cerr<<"running local_Da_eff_3..."<<endl;
	// if (q==0) return -1; // the model will not reach this step. Just comment this out.
	double G = this->local_G_3(S);
	double Da3local = S->Da3;

	// cerr<<"Da3local= "<<Da3local<<endl;
	Da3local = Da3local * is_there_redissolution(S);
	// cerr<<"Da3local*is_there_redissolution= "<<Da3local<<endl;
	//
	// if      (G>0)    return Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G3)/(1+G));
	// else if (G==0)   return Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	// else             return Da3local*(l/S->l0)*(S->q_in_0/fabs(q));

	double f3_tmp=0.0;
	if      (G>0) f3_tmp=Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G3)/(1+G));
	else if (G==0) f3_tmp=Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	else f3_tmp=Da3local*(l/S->l0)*(S->q_in_0/fabs(q));
	// cerr<<"local_Da_eff_3 with constrain = "<<f3_tmp<<endl;

	return f3_tmp;
}

/**
* This function returns the local value of Da_eff_3 without constrain (used in determining maximal redissolution volume) parameter for the pore.
* @param S pointer to the network
* @author Jingxuan Deng
* @date 22/08/2025
*/
double Pore::local_Da_eff_3_tmp(Network* S){
	// cerr<<"running local_Da_eff_3_tmp..."<<endl;
	// if (q==0) return -1; // the model will not reach this step. Just comment this out.
	double G = this->local_G_3(S);
	double Da3local = S->Da3;

	if      (G>0)    return Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G3)/(1+G));
	else if (G==0)   return Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	else             return Da3local*(l/S->l0)*(S->q_in_0/fabs(q));

	// double f3_tmp=0.0;
	// if      (G>0)    f3_tmp=Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q))*((1+S->G3)/(1+G));
	// else if (G==0)   f3_tmp=Da3local*(d/S->d0)*(l/S->l0)*(S->q_in_0/fabs(q));
	// else             f3_tmp=Da3local*(l/S->l0)*(S->q_in_0/fabs(q));
	// // cerr<<"local_Da_eff_3_tmp= "<< f3_tmp <<endl;
	// return f3_tmp;
}

/**
* This function returns the change in diameter due to dissolution in one time step.
* @param S pointer to the network
* @author Agnieszka Budek
* @date 14/03/2020
*/
double Pore::default_dd_plus(Network*S){
	// To do: add check how much E is available: use a loop to get all the E ...
	if(S->if_track_grains && !is_Va_left())  return 0;   //no reaction if there is no A species available
	if(d==0 || q ==0)  return 0;   //pore with no flow
	if(l<=S->l_min)    return 0;   //no reaction in tiny grain
    if(!is_active)     return 0;

	//dissolution parameters
	double f1      = local_Da_eff(S);
	double g       = local_G(S);
	double c0;
	if(S->if_streamtube_mixing) c0 = c_in;
	else                        c0 = calculate_inlet_cb();

	double dd_plus = 0; 		//diameter change
    //double d_tmp = min(1,d);       //possible feature for a fracture

	//finding dissolution contribution
	if      (f1==0)      dd_plus = 0;
	else if (S->G1 >=0)  dd_plus = S->dt*c0*(1-exp(-f1))/(1+g)/f1;
	else        	     dd_plus = S->dt*c0*(1-exp(-f1))/f1/d;


	return dd_plus;
}

/**
* This function returns the change in diameter due to redissolution in one time step.
* @param S pointer to the network
* @author Jingxuan Deng
* @date 07/07/2025
*/

double Pore::default_dd_plus_rediss(Network*S){
	// even if E is present, check how much will be redissolved, make sure that this amount is not larger than Ve left, do not have negative volume. check dissolution.cc file, for tracking and checking mineral A
	if(S->if_track_grains && !is_Ve_left() && !is_Ve_generated(S))  return 0;   //no reaction if there is no E species available
	if(d==0 || q ==0)  return 0;   //pore with no flow
	if(l<=S->l_min)    return 0;   //no reaction in tiny grain

	//redissolution parameters
	double f1      = local_Da_eff(S);
	double f3      = local_Da_eff_3(S); // This Damkoehler number considered all the situation when Ve0+Ve_prec<Ve_diss2_max
	double g       = local_G(S);

	double c0;
	if(S->if_streamtube_mixing) c0 = c_in;
	else                        c0 = calculate_inlet_cb();

	double dd_plus = 0; 		//diameter change


	//finding redissolution contribution
	if      (-f1-f3==0)      dd_plus = 0;
	else if (S->G3 >=0)  dd_plus = S->gamma*S->dt*(f3/(f1+f3))*c0*(1-exp(-f1-f3))/(1+g)/f1;
	else        	     dd_plus = S->gamma*S->dt*(f3/(f1+f3))*c0*(1-exp(-f1-f3))/f1/d;


	// //finding redissolution contribution
	// if      (-f1-f3==0)      dd_plus = 0;
	// else if (S->G3 >=0)  dd_plus = S->gamma*S->dt*(f3/(f1+f3))*(1-exp(-f1-f3))/(1+g)/f1;
	// else        	     dd_plus = S->gamma*S->dt*(f3/(f1+f3))*(1-exp(-f1-f3))/f1/d;

	return dd_plus;
}


/**
* This function returns the change in diameter due to precipitation in one time step (no condition for left space is checked here).
* @param S pointer to the network
* @author Agnieszka Budek
* @date 14/03/2020
*/
double Pore::default_dd_minus(Network*S){


	if(d==0 || q ==0)  return 0;   //pore with no flow
	if(l==S->l_min)    return 0;   //no reaction in tiny grain
	if(d<=S->d_min && (!is_Va_left())) return 0;
    if(!is_active)      return 0;

	//dissolution parameters
	double f1      = local_Da_eff(S);
	double g       = local_G(S);
	double c0;
	if(S->if_streamtube_mixing) c0 = c_in;
	else                        c0 = calculate_inlet_cb();


	//precipitation parameters
    double dd_minus = 0; 		//diameter change
    double f2       = local_Da_eff_2(S);
    double c0_c     = calculate_inlet_cc();
    if(S->C_eq!=0)   c0_c     = fmax(0,calculate_inlet_cc());   //irreversible reaction



	//finding precipitation contribution
	if      (f2==0)         dd_minus = 0;
	else if (f1==f2 && is_Va_left())        dd_minus = S->gamma*S->dt/(1+g)/f1*((c0_c + c0)*(1-exp(-f1)) -c0*exp(-f1)*f1);
	else if (!is_Va_left()) 				dd_minus = S->gamma*S->dt/(1+g)/f1*  c0_c*      (1-exp(-f2));
	else                    				dd_minus = S->gamma*S->dt/(1+g)/f1*(\
								       	   	   	   c0  * (f1*(1-exp(-f2)) - f2*(1-exp(-f1)))/(f1-f2)+\
												   c0_c*  (1-exp(-f2)) );


    if(!is_Ve_left() and dd_minus<0)  return 0;         //no E dissolution if there is no E left
    else                              return dd_minus ; //* is_there_precipitation(S);

}

/**
* This function removes the pore form the list of neighboring nodes.
* This function is used when deleting the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::remove_info_from_attached_nodes(){

	for(int j=0;j<2;j++) n[j]->remove_neighbor(this);
}

/**
* This function removes the pore form the list of neighboring grains.
* This function is used when deleting the pore.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::remove_info_from_attached_grains(){

	for(int j=0;j<bG;j++) g[j]->remove_Pore(this);
}

/**
* This function adds the g_tmp grain to the pore's list of neighboring grains.
* Important for merging and dynamic topology.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::add_grain (Grain *g_tmp){
	//if grain already in the list do nothing
	for(int i=0; i<bG; i++) if(g[i]==g_tmp) return;

	Grain ** g_new = new Grain* [bG+1];
	for(int i=0; i<bG; i++) g_new[i]=g[i];
	g_new[bG] = g_tmp;
	delete [] g; g = g_new;
	bG++;
}

/**
* This function removes the grain from the pore's list of neighboring grains.
* Important for merging and dynamic topology.
* @param gg grain to be removed
* @author Agnieszka Budek
* @date 25/09/2019
*/
void  Pore::remove_grain(Grain *gg){
	Grain ** g_new = new Grain * [bG];
	int b_new = 0;
	for(int bb=0;bb<bG;bb++) if(g[bb]!=gg) 	g_new[b_new++] = g[bb];

	delete [] g; g = g_new;
	bG = b_new;
}

/**
* This function returns true if a node belongs to the list of neighboring nodes.
* @param n_tmp node to be checked
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Pore::is_contain_Node  (Node *n_tmp){
	if(n[0]==n_tmp || n[1]==n_tmp) 	return true;
	else                			return false;
}

/**
* This function returns true if a grain belongs to the list of neighboring grains.
* @param g_tmp grain to be checked
* @author Agnieszka Budek
* @date 25/09/2019
*/
bool Pore::is_contain_Grain (Grain *g_tmp){
	for (int bb=0;bb<bG;bb++) if(g[bb]==g_tmp) return true;
	return false;
}

/**
* This function returns the common node of this pore and pore p1.
* @param p1 second pore
* @author Agnieszka Budek
* @date 25/09/2019
*/
Node*  Pore::find_common_node(Pore *p1){
	for(int i=0;i<2;i++) for(int j=0;j<2;j++) if(n[i]==p1->n[j]) return n[i];
	cerr<<"WARNING: Pores"<<*p1<<" and "<<*this<<" don't have common node."<<endl;
	return NULL;
}

/**
* This function returns the closest of two neighboring nodes.
* @param n_tmp second node
* @author Agnieszka Budek
* @date 25/09/2019
*/
Node*  Pore::find_closest_node(Node *n_tmp){

	double u_0 = fabs(n[0]->u - n_tmp->u);    //distance to node 0
	double u_1 = fabs(n[1]->u - n_tmp->u);    //distance to node 1
	if(u_1>u_0) return  n[0];
	else        return  n[1];

}

/**
* This function switches n_old to n_new in the list of neighboring nodes.
* @param n_old old node
* @param n_new new node
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Pore::change_pore_neighbours(Node * n_old, Node *n_new){
	//changing info about nodes in grains
	for (int i=0;i<bG;i++) g[i]->change_nodes(n_old,n_new);
	//changing info about nodes in pores
	for (int j=0; j<2;j++) if(n[j]==n_old) {
		n[j] = n_new;
		if(n[(j+1)%2]==n_new) {
			cerr<<"WARNING: In changing pore neighbors generating pore with two identical nodes."<<endl<<*this<<endl;
			//d=0;
		}
		return;}
	cerr<<"WARNING: Problem with changing neighbors in pore "<<*this<<endl;
}




double Pore::calculate_d_nbr() {

    double sum_tmp=0;
    int n_tmp=0;

    for(int i=0; i<n[0]->b;i++)
        if(n[0]->p[i]!=this) {sum_tmp+=n[0]->p[i]->d; n_tmp++;}
    for(int i=0; i<n[1]->b;i++)
        if(n[1]->p[i]!=this) {sum_tmp+=n[1]->p[i]->d; n_tmp++;}

    return sum_tmp/n_tmp;

}

double Pore::calculate_sin_angle() {
    double x0=n[0]->xy.x;
    double x1=n[1]->xy.x;
    double y0=n[0]->xy.y;
    double y1=n[1]->xy.y;

    double sqrt_tmp=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
    return abs(y0-y1)/sqrt_tmp;

}


double Pore::calculate_cos_angle() {
    double x0=n[0]->xy.x;
    double x1=n[1]->xy.x;
    double y0=n[0]->xy.y;
    double y1=n[1]->xy.y;

    double sqrt_tmp=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
    return abs(x0-x1)/sqrt_tmp;

}


std::tuple<double,double> Pore::calculate_d_nbr_direction() {

    double sum_tmp_hor=0;
    double sum_tmp_ver=0;
    int n_tmp=0;

    for(int i=0; i<n[0]->b;i++)
        if(n[0]->p[i]!=this) {
            sum_tmp_hor+=n[0]->p[i]->d*n[0]->p[i]->calculate_sin_angle();
            sum_tmp_ver+=n[0]->p[i]->d*n[0]->p[i]->calculate_cos_angle();
            n_tmp++;}
    for(int i=0; i<n[1]->b;i++)
        if(n[1]->p[i]!=this) {
            sum_tmp_hor+=n[1]->p[i]->d*n[1]->p[i]->calculate_sin_angle();
            sum_tmp_ver+=n[1]->p[i]->d*n[1]->p[i]->calculate_cos_angle();
            n_tmp++;
        }

    return {sum_tmp_hor/n_tmp,sum_tmp_ver};

}

double Pore::calculate_l_nbr() {

    double sum_tmp=0;
    int n_tmp=0;

    for(int i=0; i<n[0]->b;i++)
        if(n[0]->p[i]!=this) { sum_tmp+=n[0]->p[i]->l; n_tmp++;}
    for(int i=0; i<n[1]->b;i++)
        if(n[1]->p[i]!=this) { sum_tmp+=n[1]->p[i]->l; n_tmp++;}

    return sum_tmp / n_tmp;

}





ofstream_txt & operator << (ofstream_txt & stream, Pore &p){

	stream <<setw(12)<<p.a<<setw(12)<<p.d<<setw(12)<<p.l<<setw(8)<<p.bG<<setw(12)<<p.q;
	return stream;
}



ostream & operator << (ostream & stream, Pore &p){

	stream<<"Pore("<<p.a<<"):  tmp = "<<p.tmp<<"  n = ("<<p.n[0]->a<<","<<p.n[1]->a<<")  g = (";
	for(int b=0;b<p.bG;b++){
		stream<<p.g[b]->a;
		if(b<p.bG-1) stream<<",";
	}
	stream<<")";
    stream << "   Properties: ("<<" d = "<<p.d<<" l = "<<p.l<<" q = "<<p.q<<" is-f = "<<p.is_fracture<<")";
	return stream;

}
