#include "network.h"


/**
* This function exports a file with basic information about topology of the network.
* Data describing pore - nodes connections are saved without information about grains.
* Data are saved in the following order:
* node's name <> node's position(x,y,z) <> node's type <> no of neighbours <>   list of neighbors: (node_name, pore_name)
*
* @param oput_file_name output file name
* @author Agnieszka Budek
* @date 25/09/2019
*
*/
void Network::export_topology_file (string out_file_name){

	cerr<<"Exporting network topology..."<<endl;
	ofstream_txt* os_p = NULL;
	ofstream_txt os_tmp;
	if (out_file_name =="") os_p = &net_out;
	else{
		os_tmp.open (out_file_name, ios_base::out | ios_base::trunc );
		if(os_tmp.is_open() == false) 	{cerr<<"WARNING: Problem with reading file "<<out_file_name<<endl; return;}
		os_p = &os_tmp;
	}
	ofstream_txt &os = *os_p;


	int w_tmp_n = int(log10(NN))+2; int w_tmp_p = int(log10(NP))+2;

	os << "#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl<<endl;
	os << "NN = "<<NN<<endl;
	os << "NP = "<<NP<<endl;
	os << "NG = "<<NG<<endl;
	os << "Nx = "<<N_x<<endl;
	os << "Ny = "<<N_y<<endl;
	os << "No = "<<N_wo<<endl;
	os << "Ni = "<<N_wi<<endl<<endl;
	os << "#" <<setw(8)<<"name"<<setw(36)<<" position(x,y,z)"<<setw(8)<<"type"<<setw(7)<<"b"<<"      list of neighbors: (node_name, pore_name)"<<endl;
	os << "#  ----------------------------------------------------------------------------------"<<endl;

	for (int i=0;i<NN;i++) {
		Node * nn = n[i];
		os<<setw(8)<<nn->a<<setw(12)<<nn->xy.x<<setw(12)<<nn->xy.y<<setw(12)<<nn->xy.z<<setw(8)<<int(nn->t)<<setw(8)<< nn->b;
		for(int bb=0;bb<nn->b;bb++)   os<<setw(4)<<"("<<setw(w_tmp_n)<< nn->n[bb]->a<<","<<setw(w_tmp_p)<<nn->p[bb]->a<<") ";
		os<<endl;
	}
	if (out_file_name !="") os.close();
}

/**
* This function exports a file with additional information about topology of the network connected to the grains position in the network.
* Data are saved in the following order:
* grain's name <> Va <> Ve <> no of connected nodes <> nr of connected pores <>  list of nodes: (node_name,...)    <>  list of pores: (pore_name,...)
*
* @param oput_file_name output file name
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::export_topology_file_with_grains (string out_file_name){

	cerr<<"Exporting grains topology..."<<endl;
	ofstream_txt* os_p = NULL;
	ofstream_txt os_tmp;
	if (out_file_name =="") os_p = &net_g_out;
	else{
		os_tmp.open (out_file_name, ios_base::out | ios_base::trunc );
		if(os_tmp.is_open() == false) 	{cerr<<"WARNING: Problem with reading file "<<out_file_name<<endl; return;}
		os_p = &os_tmp;
	}
	ofstream_txt &os = *os_p;


	int w_tmp_n = int(log10(NN))+2; int w_tmp_p = int(log10(NP))+2;

	os << "#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl<<endl;
	os << "NN = "<<NN<<endl;
	os << "NP = "<<NP<<endl;
	os << "NG = "<<NG<<endl;
	os << "Nx = "<<N_x<<endl;
	os << "Ny = "<<N_y<<endl;
	os << "No = "<<N_wo<<endl;
	os << "Ni = "<<N_wi<<endl<<endl;
	os << "#" <<setw(8)<<"name"<<setw(14)<<" Va"<<setw(14)<<" Ve"<<setw(7)<<"bN"<<setw(7)<<"bP"<<"   list of nodes: (node_name,...)      list of pores: (pore_name,...)"<<endl;
	os << "#  ----------------------------------------------------------------------------------"<<endl;


	for (int i=0;i<NG;i++) {
		Grain * gg = g[i];
		os<<setw(8)<<gg->a<<setw(14)<<setprecision(5)<<gg->Va<<setw(14)<<setprecision(5)<<gg->Ve<<setw(14)<<gg->Vx<<setw(14)<<gg->bN<<setw(7)<<gg->bP;
		os << "  (";
		for(int bb=0;bb<gg->bN;bb++)   os<<setw(w_tmp_n)<< gg->n[bb]->a;
		os << ")\t\t(";
		for(int bb=0;bb<gg->bP;bb++)   os<<setw(w_tmp_p)<< gg->p[bb]->a;
		os << ")"<<endl;

	}
	if (out_file_name !="") os.close();

}



/**
* This function imports basic information about topology form a file.
* Data describing pore - nodes connections are imported without information about grains.
* Data must be saved in the following order:
* node's name <> node's position(x,y,z) <> node's type <> no of neighbors <>   list of neighbors: (node_name, pore_name)
*
* @param in_file_name input file name
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::import_topology_from_file(string in_file_name){
	cerr<<"Importing topology from file "<<in_file_name<<"."<<endl;
	ifstream in_stream;

	in_stream.open (in_file_name, ios_base::in);
	if(in_stream.is_open() == false) 	{cerr<<"WARNING: Problem with reading file "<<in_file_name<<"."<<endl; return;}


	int i=0;    //counting lines
	int j=0;    //counting nodes
	int i_wi=0; //counting inlet nodes
	int i_wo=0; //counting outlet nodes
	string s;
	b_max = 0;     //maximal nr of node neighbors
	int ** n_tmp;   //temporary table to save info of node neighbors
	int ** p_tmp;   //temporary table to save info of pore-node connections

	while(getline(in_stream,s)){
		i++;

		if(s.size()<3) 	continue;
		if(s[0]=='#') 	continue;

		if(s[0]=='N'){ //reading info about NN and NP
			istringstream line(s);

			string name, e;
			double value = -1;

			line >> name >> e >> value;

			if(e!="=" || value==-1) {cerr<<"WARNING: Problem with parsing line nr "<<i<<" in file "<<in_file_name <<"."<<endl; continue;}

			if(name == "NP"){
				NP = value;
				p  = new Pore*[NP];
				for(int ii=0;ii<NP;ii++) p[ii] = new Pore(d0,l0,ii);
				cerr<<"Setting NP = "<<NP<<endl;}

			else if(name == "NN"){
				NN = value;
				n  = new Node*[NN];
				n_tmp = new int*[NN];  p_tmp = new int*[NN];
				cerr<<"Setting NN = "<<NN<<endl;}

			else if(name == "NG"){
				NG = value;
				g  = new Grain*[NG];
				cerr<<"Setting NG = "<<NG<<endl;}

			else if(name == "Ni"){
				N_wi = value;
				Q_tot = N_wi*2;
				wi = new Node* [N_wi];
				cerr<<"Setting N_wi = "<<N_wi<<endl;}

			else if(name == "No"){
				N_wo = value;
				wo = new Node* [N_wo];
				P_in   = N_wo;
				cerr<<"Setting N_wo = "<<N_wo<<endl;}

			else if(name == "Nx"){
				N_x = value;
				cerr<<"Setting N_x = "<<N_x<<endl;}

			else if(name == "Ny"){
				N_y = value;
				cerr<<"Setting N_y = "<<N_y<<endl;}


			else  {cerr<<"WARNING: Problem with parsing line nr "<<i<<" in file "<<in_file_name <<"."<<endl; continue;}
		}


		else {  //reading info about node
			if(j>=NN) {cerr<<"WARNING: To many nodes to be read in "<< in_file_name<<endl; return;}
			istringstream line(s);

			int name, t_tmp, b_tmp;
			double x_tmp, y_tmp, z_tmp;

			if(line >> name >> x_tmp >> y_tmp >> z_tmp >> t_tmp >>b_tmp){
				if (b_tmp>b_max) b_max = b_tmp;  //updating maximal nr of neighbors
				n[j] = new Node(name,b_tmp,t_tmp,Point(x_tmp,y_tmp,z_tmp));
				n_tmp[j] = new int[b_tmp]; p_tmp[j] = new int[b_tmp];
				//adding inlet/outlet nodes
				if(t_tmp == 1 ) { //setting inlet nodes
					wi[i_wi++] = n[j];
					if (i_wi>N_wi) {cerr << "To many inlet nodes i_wi = "  << i_wi << endl; exit(1);}}
				if(t_tmp == -1) { //setting outlet nodes
					wo[i_wo++] = n[j];
					if (i_wo>N_wo) {cerr << "To many outlet nodes i_wo = " << i_wo << endl; exit(1);}}
				line.clear();
			}
			else {cerr<<"WARNING: Problem with parsing beginning of line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}


			for(int bb=0;bb<b_tmp;bb++){
				char c1,c2,c3;
				int n_name, p_name;
				if(line>>c1>>n_name>>c2>>p_name>>c3){
					if(n_name>=NN || n_name<0) cerr<<"Wrong nr of node name = "<<n_name<<endl<<flush;
					if(p_name>=NP || p_name<0) cerr<<"Wrong nr of pore name = "<<p_name<<endl<<flush;
					if(c1!='(' || c2 !=',' || c3 != ')'){
						cerr<<"ERROR: Problem with reading pore neighbors in line "<<i<<" from file "<< in_file_name<< "."<<endl<<flush;
						exit(1);}
					n_tmp[j][bb] = n_name;  p_tmp[j][bb] = p_name;
					line.clear();
				}
				else{cerr<<"ERROR: Problem with parsing second part of line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}
			}
			j++;
		}
	}

	//adding info of connections
	for (int j=0;j<NN;j++){
		for(int bb=0;bb<n[j]->b;bb++){
			n[j]->n[bb] = n[n_tmp[j][bb]];
			n[j]->p[bb] = p[p_tmp[j][bb]];
			Pore* pp = p[p_tmp[j][bb]];
			if (pp->n[0] == NULL) {pp->n[0] = n[j]; pp->n[1] = n[n_tmp[j][bb]];}
		}
	}

	//free memory
	for(int ii=0;ii<NN;ii++){ delete n_tmp[ii]; delete p_tmp[ii];}
	delete n_tmp; delete p_tmp;



	cerr<<"Topology has been imported."<<endl<<endl;

	export_topology_file ("topology_tmp.out");

}



/**
* This function imports an additional information about topology of the network connected to the grains position.
* Data are saved in the following order:
* grain's name <> Va <> Ve <> no of connected nodes <> nr of connected pores <>  list of nodes: (node_name,...)    <>  list of pores: (pore_name,...)
*
* @param in_file_name input file name
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::import_grains_from_file (string in_file_name){
	cerr<<"Importing information of grains from file "<<in_file_name<<"."<<endl;
	ifstream in_stream;

	in_stream.open (in_file_name, ios_base::in);
	if(in_stream.is_open() == false) 	{
		cerr<<"WARNING: Problem with reading file "<<in_file_name<<"."<<endl;
		cerr<<"All pores have diameters set to d0 = "<<d0<<" and lengths set to l0 = "<<l0<<endl;
		return;}


	int i=0;    //counting lines
	int j=0;    //counting grains

	string s;

	while(getline(in_stream,s)){
		i++;

		if(s.size()<3) 	continue;
		if(s[0]=='#') 	continue;
		if(s[0]=='N')   continue;

		if(j>=NG) {cerr<<"WARNING: To many grains to be read in "<< in_file_name<<endl; return;}
		istringstream line(s);

		int name, bN_tmp, bP_tmp;
		double Va_tmp, Ve_tmp, Vx_tmp;

		if(line >> name >> Va_tmp >> Ve_tmp >> Vx_tmp >> bN_tmp >> bP_tmp) {
			g[j] = new Grain(name,Va_tmp,Ve_tmp,Vx_tmp,bN_tmp, bP_tmp);
			line.clear();
		}
		else {cerr<<"ERROR: Problem with parsing beginning of line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}
		//reading info about nodes
		char c;
		line >> c;
		if(c!='(') {cerr<<"ERROR: Problem 1. with parsing line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}
		for(int bb=0;bb<bN_tmp;bb++){
			int n_name;
			if(line>>n_name){
				if(n_name>=NN || n_name<0) cerr<<"Wrong nr of node name = "<<n_name<<endl<<flush;
				g[j]->n[bb] = n[n_name];
				//line.clear();
			}
			else{cerr<<"ERROR: Problem 2. with parsing second part of line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}
		}
		line >> c;
		if(c!=')') {cerr<<"ERROR: Problem 3. with parsing line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}

		//reading info about pores
		line >> c;
		if(c!='(') {cerr<<"ERROR: Problem 4. with parsing line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}
		for(int bb=0;bb<bP_tmp;bb++){
			int p_name;
			if(line>>p_name){
				if(p_name>=NP || p_name<0) cerr<<"Wrong nr of node name = "<<p_name<<endl<<flush;
				g[j]->p[bb] = p[p_name];
				Pore * pp = p[p_name];
				bool tmp = false;
				for(int s = 0;s<pp->bG;s++) if(pp->g[s]==NULL) {pp->g[s] = g[j]; tmp=true; break;}
				if(! tmp) cerr<<"Problem with adding info about pores in line "<<i<<"."<<endl;
				//line.clear();
			}
			else{cerr<<"ERROR: Problem with parsing second part of line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}
		}
		line >> c;
		if(c!=')') {cerr<<"ERROR: Problem with parsing beginning of line "<<i<< " in file "<<in_file_name<<"."<<endl<<flush; exit(1); continue;}
		j++;
	}


	if(printing_mode == "grains") save_info_about_initial_node_positions ();  //used only for printing ps/pdf files in grain style
	NN_max = NN; NP_max = NP; NG_max = NG;

	cerr<<"Information about grains has been imported."<<endl;
	export_topology_file_with_grains ("topology_g_tmp.out");
}



/**
* This function imports information about pores sizes from external file.
* @param input_file_name input file name
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::import_pore_size_from_file (string in_file_name){
	cerr<<endl<<"Reading pore sizes from file "<<in_file_name<<"."<<endl;

	ifstream in_stream;

	in_stream.open (in_file_name, ios_base::in);
	if(in_stream.is_open() == false) 	{cerr<<"ERROR: Problem with reading file "<<in_file_name<<"."<<endl; return;}

	int i = 0; //counting pores
	string s;

	while(getline(in_stream,s)){

		if(s.size()<3) 	continue;
		if(s[0]=='#') 	continue;

		stringstream line;
		line << s;
		if(i>=NP) break;  //All pores have been already read
		int name,b_tmp;
		double d_tmp, l_tmp;
		line >> name >> d_tmp >> l_tmp >>b_tmp;
		if (name>=NP || name<0) {cerr << "WARNING: Problem with pore name: name = "<<name<<endl; continue;}
		p[name]->l = l_tmp;
		p[name]->d = d_tmp;

		//new information about nr of grains
		p[name]->bG = b_tmp;
		delete [] p[name]->g;
		p[name]->g = new Grain*[b_tmp];
		for (int i=0;i<b_tmp;i++) p[name]->g[i]=NULL;

		i++;
	}
	if(i!=NP) cerr<<"WARNING: I read "<<i<<" pores. The total nr of pores is "<<NP<<"."<<endl;

	cerr<<"Pore size has been imported."<<endl;
	//cerr<<"The value of l0 has not been recalculated and is set to: "<<l0<<"."<<endl<<endl;
	//Warto by tu uaktualnic srednie l0 i d0 !!!
}

/**
* This function prints all data about network current properties.
* This data includs information about:
*
*  - pores: name, diameter, length, nr of neighboring grains, flow, Da_local, G_local
*
*  - nodes: name, type, presure, concentration of spacies B, concentration of species C
*
*  - grains: name, volube of A, volume of E
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::print_net_txt(){


	pores_out <<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	nodes_out <<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	grains_out<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;

	pores_out<<"#" <<setw(11)<<"name"<<setw(12)<<"d"<<setw(12)<<"l"<<setw(8)<<"b"<<setw(12)<<"q"<<setw(12)<<"Da_local"<<setw(12)<<"G_local"<<endl;
	pores_out<<"#  ----------------------------------------------------------------------"<<endl;

	nodes_out<<"#" <<setw(11)<<"name"<<setw(6)<<"type"<<setw(12)<<"u"<<setw(12)<<"cb"<<setw(12)<<"cc"<<endl;
	nodes_out<<"#  ------------------------------------------------"<<endl;

	grains_out<<"#" <<setw(11)<<"name"<<setw(12)<<"Va"<<setw(12)<<"Ve"<<setw(12)<<"Vx"<<endl;
	grains_out<<"#  ------------------------------------------------"<<endl;


	for(int i=0;i<NP;i++) {
		pores_out <<*(p[i]);
		//Da_local and G_local can be printed if needed
		//pores_out<<setw(12)<<p[i]->local_Da_eff(this)<<setw(12)<<p[i]->local_G(this)<<endl;
		pores_out<<endl;
	}
	for(int i=0;i<NN;i++) nodes_out  << *(n[i]);
	for(int i=0;i<NG;i++) grains_out << *(g[i]);


	pores_out  <<endl;
	nodes_out  <<endl;
	grains_out <<endl;
}

/**
* This function prints data like flow field, diameters, pressure field in table (works best for hexagonal network).
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::print_tables_txt(){

	diameters_out      <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	flow_out           <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	concentration_out  <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	concentration2_out <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	pressure_out       <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	lengths_out        <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;

	if(if_track_grains){
		VA_out		       <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
		VE_out      	   <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
		VX_out      	   <<endl<<endl<<fixed<<"#" << tot_steps<<". step of evolution: T_tot =  "<<tot_time<<endl;
	}

	if (type_of_topology == "hexagonal"){
		for(int i=0;i<NN;i++) {
			for(int b=0; b<3;b++) diameters_out  <<setprecision(7)<<setw(12)<<p[i*3+b]->d;
			for(int b=0; b<3;b++) flow_out       <<setprecision(7)<<setw(12)<<p[i*3+b]->q;
			for(int b=0; b<3;b++) lengths_out    <<setprecision(7)<<setw(12)<<p[i*3+b]->l;
			pressure_out       <<setprecision(7)<<setw(20)<<n[i]->u;
			concentration_out  <<setprecision(7)<<setw(12)<<n[i]->cb;
			concentration2_out <<setprecision(7)<<setw(12)<<n[i]->cc;
			for(int b=0; b<2;b++)VA_out  <<setprecision(7)<<setw(12)<<g[2*i+b]->Va;
			for(int b=0; b<2;b++)VE_out  <<setprecision(7)<<setw(12)<<g[2*i+b]->Ve;
			for(int b=0; b<2;b++)VX_out  <<setprecision(7)<<setw(12)<<g[2*i+b]->Vx;
			if(i%N_x ==N_x-1){
				diameters_out       <<endl;
				flow_out            <<endl;
				concentration_out   <<endl;
				concentration2_out  <<endl;
				pressure_out        <<endl;
				VA_out              <<endl;
				VE_out              <<endl;
				VX_out              <<endl;
				lengths_out         <<endl;
			}
		}
	}

	else if (type_of_topology == "diamond" || type_of_topology == "square"){
			for(int i=0;i<NN;i++) {
				for(int b=0; b<2;b++) diameters_out  <<setprecision(7)<<setw(12)<<p[i*2+b]->d;
				for(int b=0; b<2;b++) flow_out       <<setprecision(7)<<setw(12)<<p[i*2+b]->q;
				for(int b=0; b<2;b++) lengths_out    <<setprecision(7)<<setw(12)<<p[i*2+b]->l;
				pressure_out       <<setprecision(7)<<setw(20)<<n[i]->u;
				concentration_out  <<setprecision(7)<<setw(12)<<n[i]->cb;
				concentration2_out <<setprecision(7)<<setw(12)<<n[i]->cc;
				VA_out  <<setprecision(7)<<setw(12)<<g[i]->Va;
				VE_out  <<setprecision(7)<<setw(12)<<g[i]->Ve;
				VX_out  <<setprecision(7)<<setw(12)<<g[i]->Vx;
				if(i%N_x ==N_x-1){
					diameters_out       <<endl;
					flow_out            <<endl;
					concentration_out   <<endl;
					concentration2_out  <<endl;
					pressure_out        <<endl;
					VA_out              <<endl;
					VE_out              <<endl;
					VX_out              <<endl;
					lengths_out         <<endl;
				}
			}
		}



	else {   //printing nodes line by line,  mean d and q around given node are printed.
		int linia = 0;
		for(int i=0;i<NN;i++) {
			if(int(n[i]->xy.y) > linia){
				linia++;
				diameters_out      <<endl;
				flow_out           <<endl;
				concentration_out  <<endl;
				concentration2_out <<endl;
				pressure_out       <<endl;
				lengths_out        <<endl;
			}
			double d_mean = 0; double q_mean = 0; double l_mean = 0;
			for(int b=0;b<n[i]->b;b++){ d_mean+= n[i]->p[b]->d;  q_mean+= n[i]->p[b]->q; l_mean+=n[i]->p[b]->l;}
			diameters_out      <<setprecision(5)<<setw(10)<<d_mean/n[i]->b;
			flow_out           <<setprecision(5)<<setw(10)<<q_mean/n[i]->b;
			pressure_out       <<setprecision(5)<<setw(12)<<n[i]->u;
			lengths_out        <<setprecision(7)<<setw(12)<<l_mean/n[i]->b;
			if (if_streamtube_mixing)   concentration_out  <<setprecision(5)<<setw(12)<<p[i]->c_in<<"\t"<<setprecision(5)<<setw(12)<<p[NN+i]->c_in;
			else						concentration_out  <<setprecision(5)<<setw(12)<<n[i]->cb;
			concentration2_out <<setprecision(5)<<setw(12)<<n[i]->cc;
		}
	}
}

/**
* This is a main function saving all data. Is is run in each time step.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::save_all_data(bool if_save_now){


	static double Va_old = 0;

	//deciding either to save or not
	if     (s_save_data<0)                                   if_save_now = check_diss_front(print_diss_factor, pages_saved*abs(s_save_data));
	else if(s_save_data>=1 && tot_steps%int(s_save_data)==0) if_save_now = true;
	else if(s_save_data>0 && s_save_data<1){
//check the volume condition!!!
		if(tot_steps==0) Va_old = Va_tot;
		if((Va_old-Va_tot)/Va_tot>s_save_data){
			if_save_now = true;
			Va_old = Va_tot;
		}
	}


	if(if_save_now){
		cerr<<"Saving all data..."<<endl;
		description_note = "At the end of evolution step: s = " + to_string(tot_steps);
		if(if_save_ps)        net_ps<< *this;
		if(if_save_txt)       print_net_txt();
		if(if_save_table)     print_tables_txt();
		if(if_save_topology) {export_topology_file(); export_topology_file_with_grains();}
		pages_saved++;
	}


}

