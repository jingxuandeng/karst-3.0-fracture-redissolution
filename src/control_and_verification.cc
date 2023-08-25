#include "network.h"

/**
* This function checks flow balance in the system.
* If the total inlet flow (flow coming through inlet pores) is not equal total outlet flow the warning is printed.
* The problem with flow balance can be a result of numerical problems with solving linear equations for pressure.
* This was a reason to try to calculate pressure in two steps.
*
* @author Agnieszka Budek
* @date 26/09/2019
*/
void Network::check_flow_balance(){

	double eps = 10e-2;

	double Q_in_tmp  = 0;
	double Q_out_tmp = 0;


	if (if_smarter_calculation_of_pressure && false){ //unused, to be deleted
		check_diss_pattern(d_max_for_u);
		//calculate total input flow
		for(int i=0;i<NP;i++){
			if((p[i]->n[0]->x>0 && (p[i]->n[1]->x==0||p[i]->n[1]->t==-1)) || (p[i]->n[1]->x>0 && (p[i]->n[0]->x==0||p[i]->n[0]->t==-1)))
				if(p[i]->d>0)Q_in_tmp+=fabs(p[i]->q);
		}
	}
	else{
		//calculate total input flow
		for(int i=0;i<N_wi;i++){
			Node* n_tmp = wi[i];
			for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) Q_in_tmp+=fabs(n_tmp->p[j]->q);
		}
	}
	//calculate total output flow
	for(int i=0;i<N_wo;i++){
		Node* n_tmp = wo[i];
		for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) Q_out_tmp+=fabs(n_tmp->p[j]->q);
	}

	//print_network_for_debugging ("In calculating flow balance ","x","diameter",  "name");


	if(fabs(Q_in_tmp-Q_out_tmp)/(Q_in_tmp+Q_out_tmp) > eps)
		cerr<<"WARNING: Flow is not conserved: Q_in = "<<Q_in_tmp <<" Q_out = "<<Q_out_tmp<<"."<<endl;
	else
		cerr<<"Flow is conserved: Q_in = "<<Q_in_tmp <<" Q_out = "<<Q_out_tmp<<"."<<endl;
	if(fabs(Q_in_tmp-Q_out_tmp)/(Q_in_tmp+Q_out_tmp) > 1000*eps) exit(888);
}

/**
* This function checks acid balance.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::check_acid_balance(){

	double eps = 10e-7;
	double VB_in=0, VB_out =0;

	if(if_streamtube_mixing){

		//calculate total input of acid
		for(int i=0;i<N_wi;i++){
			for (int j=0; j<wi[i]->b;j++){
				Pore * pp = wi[i]->p[j];
				if(pp->d>0) VB_in+=fabs(pp->q)*pp->c_in*dt/dt_unit;
			}
		}

		//calculate total output of acid
		for(int i=0;i<N_wo;i++){
			for (int j=0; j<wo[i]->b;j++){
				Pore * pp = wo[i]->p[j];
				if(pp->d>0) VB_out+=fabs(pp->q)*pp->c_in*outlet_c_b(pp)*dt/dt_unit;
			}
		}
	}


	else{

		//calculate total input of acid
		for(int i=0;i<N_wi;i++){
			Node* n_tmp = wi[i];
			for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) VB_in+=fabs(n_tmp->p[j]->q)*n_tmp->cb*dt/dt_unit;
		}

		//calculate total output of acid
		for(int i=0;i<N_wo;i++){
			Node* n_tmp = wo[i];
			for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) VB_out+=fabs(n_tmp->p[j]->q)*n_tmp->cb*dt/dt_unit;
			}
	}

	//calculate consumption of acid
	double Va_tot_tmp = 0;
	for(int i=0;i<NG;i++) Va_tot_tmp+= g[i]->Va;

	double Va_delta = Va_tot - Va_tot_tmp;
	double Vb_delta = VB_in - VB_out;

	if(fabs(Vb_delta - Va_delta)/fabs(Vb_delta + Va_delta) > eps)
		{cerr<<"WARNING: Mass is not conserved: Vb_delta = "<<setprecision(10)<<Vb_delta <<" Va_delta = "<<Va_delta<<setprecision(10)<<"."<<endl;}
	else
		cerr<<"Mass is conserved: Vb_delta = "<<Vb_delta <<" Va_delta = "<<Va_delta<<"."<<endl;
	if(!if_precipitation) Va_tot=Va_tot_tmp;

}


void Network::check_precipitating_balance(){


	double eps = 10e-3;
	double VC_in=0, VC_out =0;

	//calculate total production of C
	double Va_tot_tmp = 0;
	for(int i=0;i<NG;i++) Va_tot_tmp+= g[i]->Va;

	double Va_delta = Va_tot - Va_tot_tmp;
	VC_in = Va_delta*gamma; //cerr<<"VC_in = "<<VC_in<<endl;

	//calculate total input of precipitatin species
	for(int i=0;i<N_wi;i++){
		Node* n_tmp = wi[i];
		for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) VC_in+=fabs(n_tmp->p[j]->q)*n_tmp->cc*dt/dt_unit;
	}

	//calculate total output of acid
	for(int i=0;i<N_wo;i++){
		Node* n_tmp = wo[i];
		for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) VC_out+=(1.*gamma)*fabs(n_tmp->p[j]->q)*n_tmp->cc*dt/dt_unit; //UWAGA: to ewentualnie pomnozym przez gamma i kappa trzeba lub nie
		}

	//calculate consumption of C
	double Ve_tot_tmp = 0;
	for(int i=0;i<NG;i++) Ve_tot_tmp+= g[i]->Ve;

	double Ve_delta =  Ve_tot_tmp - Ve_tot;

	//tymczasowe smieci
	double d_min = 1;
	for(int i=0;i<NP;i++)   if (p[i]->d<d_min && p[i]->d>0) d_min = p[i]->d;
//	cerr<<"d_min = "    << d_min<<endl;
//	cerr<<"VC_in = "    << VC_in<<endl;
//	cerr<<"VC_out = "   << VC_out<<endl;
//	cerr<<"Ve_stare = " << Ve_tot<<endl;
//	cerr<<"Ve_nowe =  " << Ve_tot_tmp<<endl;



	double Vc_delta =  VC_in - VC_out;
	if(fabs(Vc_delta - Ve_delta)/fabs(Vc_delta + Ve_delta) > eps)
		{cerr<<"WARNING: Mass is not conserved in precipitation: Vc_delta = "<<setprecision(10)<<Vc_delta <<" Ve_delta = "<<Ve_delta<<setprecision(10)<<"."<<endl;}
	else
		cerr<<"Mass is conserved in precipitation: Vc_delta = "<<Vc_delta <<" Ve_delta = "<<Ve_delta<<"."<<endl;
	Ve_tot=Ve_tot_tmp;
	Va_tot=Va_tot_tmp;

}


/**
* This function checks network connections consistency.
* It is useful for debugging and merging.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::check_network_connections(){

	cerr<<"Checking network structure..."<<endl;
	bool if_ok=true;

	//checking nodes
	//cerr<<"Checking nodes:"<<endl;
	for (int i=0;i<NN;i++){
		//cerr<<"Checking: "<<*n[i]<<endl<<flush;
		Node *nn = n[i];
		if (nn == NULL) {cerr<<"ERROR Node does not exist: i = "<<i<<endl<<flush; if_ok=false;}
		if (nn->b<1)    {cerr<<"ERROR Node is not connected to the network: "<<*nn<<endl<<flush; if_ok=false;}
		if (nn->a<0)    {cerr<<"ERROR: Negative node in the main list of grains: "<<*nn<<endl; if_ok=false;}
		for(int b=0;b<nn->b;b++) {
			if (nn->n[b]==NULL) {cerr<<"ERROR: in checking network's nodes: n["<<i<<"]->n["<<b<<"] = NULL."<<endl; if_ok=false;}
			if (nn->p[b]==NULL) {cerr<<"ERROR: in checking network's nodes: n["<<i<<"]->p["<<b<<"] = NULL."<<endl; if_ok=false;}
			if (nn->n[b]->a<0)  {cerr<<"ERROR: in checking network's nodes: "<<*nn<<endl<<*nn->n[b]<<endl; if_ok=false;}
			if (nn->p[b]->a<0)  {cerr<<"ERROR: in checking network's nodes: "<<*nn<<endl<<*nn->p[b]<<endl; if_ok=false;}
			if (!nn->p[b]->is_contain_Node(nn)) {cerr<<"ERROR: in checking network's nodes: "<<*nn<<endl<<*nn->p[b]<<endl;if_ok=false;}
			if (!nn->n[b]->is_connected_to_Node(nn)) {cerr<<"ERROR: in checking network's nodes: "<<*nn<<endl<<*nn->n[b]<<endl;if_ok=false;}
			if (!((nn->p[b]->n[0]==nn && nn->p[b]->n[1]==nn->n[b]) || (nn->p[b]->n[1]==nn && nn->p[b]->n[0]==nn->n[b])))
				{cerr<<"ERROR: in checking nodes: Wrong pores/nodes structure for: "<<*nn<<" and "<<*nn->p[b]<<endl;if_ok=false;}
		}
		if (if_track_grains) for(int b=0;b<nn->bG;b++) {
			if (nn->g[b]==NULL) {cerr<<"ERROR: in checking network's nodes: n["<<i<<"]->g["<<b<<"] = NULL."<<endl;if_ok=false;}
			if (nn->g[b]->a<0)  {cerr<<"ERROR: in checking network's nodes: "<<*nn<<endl<<*nn->n[b]<<endl;if_ok=false;}
		}
	}



	//checking pores
	//cerr<<"Checking pores:"<<endl;
	for (int i=0;i<NP;i++){
		Pore *pp = p[i];
		if (pp == NULL) {cerr<<"ERROR: Pore does not exist: i = "<<i<<endl<<flush;if_ok=false;}
		if (pp->a<0)    {cerr<<"ERROR: Negative pore in the main list of pores: "<<*pp<<endl;if_ok=false;}
		if (pp->n[0]==pp->n[1]) {cerr<<"ERROR: Pathological pore: "<<*pp<<endl;if_ok=false;}
		for(int b=0;b<2;b++){
			if (pp->n[b]==NULL) {cerr<<"ERROR: in checking network's pores: p["<<i<<"]->n["<<b<<"] = NULL."<<endl<<flush;if_ok=false;}
			if (pp->n[b]->a<0)  {cerr<<"ERROR: in checking network's pores: "<<*pp<<endl<<flush;if_ok=false;}
		}
		if (if_track_grains) for(int b=0;b<pp->bG;b++){
			if (pp->g[b]==NULL) {cerr<<"ERROR: in checking network's pores: p["<<i<<"]->g["<<b<<"] = NULL."<<endl<<flush;if_ok=false;}
			if (pp->g[b]->a<0)  {cerr<<"ERROR: in checking network's pores: "<<*pp<<endl<<flush;if_ok=false;}
			int check=0;
			Grain * gg=pp->g[b];
			for(int bb=0; bb<gg->bP;bb++) if (gg->p[bb]==pp) check++;
			if(check!=1) {cerr<<"ERROR: in checking network's pores: Wrong pores/grains structure for : "<<*pp<<endl<<*gg<<endl<<flush;if_ok=false;}
			if(pp->bG<1) cerr<<"WARNING: Pore "<<i<<" is connected to unusual nr of grains: "<<*pp<<endl<<flush;
		}
		//optional checking for duplicate pores:
		if(type_of_merging=="none" and false) for(int j=0;j<NP;j++) if(p[j]!=pp)
			if(((pp->n[0]==p[j]->n[0])&&(pp->n[1]==p[j]->n[1])) ||\
			   ((pp->n[1]==p[j]->n[0])&&(pp->n[0]==p[j]->n[1])))
				{cerr<<"ERROR: in checking network's pores: Duplicate pores: "<<i<<" and"<<j<<endl;if_ok=false;}
	}

	//checking grains
	//cerr<<"Checking grains:"<<endl;
	if (if_track_grains) for (int i=0;i<NG;i++){
		Grain *gg = g[i];
		if (gg == NULL) {cerr<<"ERROR: Grain does not exist: i = "<<i<<endl<<flush;if_ok=false;}
		if (gg->a<0)    {cerr<<"ERROR: Negative grain in the main list of grains:"<<*gg<<endl;if_ok=false;}
		for(int b=0;b<gg->bN;b++){

			if (gg->n[b]==NULL) {cerr<<"ERROR: in checking network's grains: "<<*gg<<endl;if_ok=false;}
			if (gg->n[b]->a<0)  {cerr<<"ERROR: in checking network's grains: "<<*gg<<endl<<*gg->n[b]<<endl;if_ok=false;}
		}
		for(int b=0;b<gg->bP;b++){
			if (gg->p[b]==NULL) {cerr<<"ERROR: in checking network's grains: "<<*gg<<endl;if_ok=false;}
			if (gg->p[b]->a<0)  {cerr<<"ERROR: in checking network's grains: "<<*gg<<endl<<*gg->p[b]<<endl;if_ok=false;}
			int check=0;
			for(int bb=0; bb<gg->p[b]->bG;bb++)
				if (gg->p[b]->g[bb]==gg) check++;
			if(check!=1) {cerr<<"ERROR: in checking network's grains : " <<*gg<<endl<<*gg->p[b]<<endl<<"Check = "<<check<<endl;if_ok=false;}
			if(!gg->do_contain_Node(gg->p[b]->n[0]) || !gg->do_contain_Node(gg->p[b]->n[1]))
					{cerr<<"ERROR: in checking network's grains: "<<*gg<<endl<<*gg->p[b]<<endl;if_ok=false;}
		}
		//optionally chacking if we don't have duplicates (for grain = triangles)
		if(false) for(int j=0;j<NG;j++) if(g[j]!=gg) if(g[j]->bN==gg->bN){
			Grain* t1 = gg; Grain* t2 = g[j];
			bool if_tmp=true;
			for(int b1=0;b1<t1->bN;b1++){
				bool if_tmp_in=false;
				for(int b2=0;b2<t2->bN;b2++) if(t1->n[b1] == t2->n[b2]) if_tmp_in=true;
				if(!if_tmp_in) if_tmp=false;}
			bool if_tmp2=true;
			for(int b1=0;b1<t1->bP;b1++){
				bool if_tmp_in=false;
				for(int b2=0;b2<t2->bP;b2++) if(t1->p[b1] == t2->p[b2]) if_tmp_in=true;
				if(!if_tmp_in) if_tmp2=false;}
			if(if_tmp && if_tmp2)  {cerr<<"ERROR: in checking network's grains: Duplicate grains: "<<endl<<*t1<<endl<<*t2<<endl;if_ok=false;}
		}
	}

	if(if_ok)   cerr<< "Network structure is ok."     <<endl;
	else      { cerr<< "Network structure is NOT ok." <<endl; exit(666);}

}


void Network::check_GMash_connections(){
//	double threshold = 2;
//	check_diss_pattern(threshold);
//	check_if_dissolved();
//
//	for(int i=0;i<NP;i++) if(p[i]->d > threshold*d0 && p[i]->x==1){
//		cerr<<"Problematic pore: " <<*p[i]<<endl;
//
//		for (int b=0;b<2;b++){
//			Node * nn = p[i]->n[b];
//			cerr<<"Problematic node nr "<<b<<":  "<<*nn<<endl;
//			cerr<<"Node properties: "<<setw(8)<<nn->a<<setw(12)<<nn->xy.x<<setw(12)<<nn->xy.y<<setw(12)<<nn->xy.z<<setw(8)<<int(nn->t)<<setw(8)<< nn->b<<endl;
//
//		}
//
//		exit(666);
//	}

	for(int i=0;i<NN;i++) if(n[i]->b>10) if(n[i]->xy.x>10 && n[i]->xy.x<90 &&n[i]->xy.y>10&&n[i]->xy.y<90) {

		//Problematic pore: 2129

		cerr<<"Probelm with node "<<i<<endl;
		Node *nn=n[i];

		for(int i=0;i<NP;i++) p[i]->tmp = p[i]->d;
		for(int bb=0;bb<nn->b;bb++){
			//nn->n[bb]->t=6;
			//cerr<<"Node nr "<<bb<<" : "<<*nn->n[bb]<<endl;
			cerr<<"Pore nr "<<bb<<" : "<<*nn->p[bb]<<endl;
		}
		cerr<<"Problematic pore: "<<*nn<<endl;

		for(int i=0;i<NN;i++) n[i]->tmp = n[i]->a;
		for(int i=0;i<NP;i++) p[i]->tmp = 666;
		description_note = "In control gmash: s = " + to_string(tot_steps);
		Print_network_in_debugging_style_tmp   (net_ps,*this,i);
	}

	system("ps2pdf net.ps");
	system("open net.pdf");

	exit(666);


}

/**
* This function prints the full information about podes and pores to the stderr.
* Useful for debugging.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::print_network_for_debugging (string text, string type_n, string type_p, string type_g){

	if(if_debugging_printing && tot_steps%int(s_save_data)==0){

		if (type_n == "name")                for(int i=0;i<NN;i++) n[i]->tmp = n[i]->a;
		if (type_n == "pressure")            for(int i=0;i<NN;i++) n[i]->tmp = n[i]->u;
		if (type_n == "concentration")       for(int i=0;i<NN;i++) n[i]->tmp = n[i]->cb;
		if (type_n == "acid concentration")  for(int i=0;i<NN;i++) n[i]->tmp = n[i]->cb;
		if (type_n == "concentration C")     for(int i=0;i<NN;i++) n[i]->tmp = n[i]->cc;
		if (type_n == "type") 	             for(int i=0;i<NN;i++) n[i]->tmp = n[i]->t;
		if (type_n == "x")   	             for(int i=0;i<NN;i++) n[i]->tmp = n[i]->x;

		if (type_p == "name")                for(int i=0;i<NP;i++) p[i]->tmp = p[i]->a;
		if (type_p == "flow")                for(int i=0;i<NP;i++) p[i]->tmp = p[i]->q;
		if (type_p == "length")              for(int i=0;i<NP;i++) p[i]->tmp = p[i]->l;
		if (type_p == "diameter")            for(int i=0;i<NP;i++) p[i]->tmp = p[i]->d;
		if (type_p == "concentration")       for(int i=0;i<NP;i++) p[i]->tmp = p[i]->c_in;

		if (type_g == "name")                for(int i=0;i<NG;i++) g[i]->tmp = g[i]->a;
		if (type_g == "volume A")            for(int i=0;i<NG;i++) g[i]->tmp = g[i]->Va;

		description_note = text + ": s = " + to_string(tot_steps);
		net_ps<<*this;}
}
