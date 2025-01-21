/**
 * @file
 * @author  Agnieszka Budek
 * @date 25/09/2019
 *
 * @section DESCRIPTION
 *
 * In this file you will find functions connected to the merging.
 *
 */

#include <functional>
#include "network.h"
#include "printing.h"
#include "deque"



/**
* This is a main function doing merging.
* It checks what type of merging to execute.
* If type_of_merging is "none", then no merging is done.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::do_merging(){

	//merging
	if      (type_of_merging == "none") 					cerr<<"No merging."<<endl;
	else if (type_of_merging == "merge_empty_grains") 		merge_empty_grains();
    else if (type_of_merging == "merge_for_fracture") 		merge_for_fracture();
    else                               						cerr<<"WARNING: type_of_merging: "<<type_of_merging<<" not implemented yet."<<endl;

	//checking network connections
	if(type_of_merging != "none" and type_of_merging != "merge_for_fracture")     {
        //clear_unconeccted_pores();
        check_network_connections();}


}

/*
 * WARNING:
 * Works worse than a traditional merging.
 * Do not use this function without implementation of geometrical prediction
 * for the grain that should be in the contact with the fracture!!!
 */
void Network::merge_for_fracture() {

    cerr<<"Merging for fracture..."<<endl;

    for (int i = 0; i < NG; i++) {
        Grain *gg = g[i];
        Pore *pp = NULL;
        deque<Grain*> g_list;

        //checking if the grain belongs to the fracture
        bool is_fracture = false;
        for (int j = 0; j < gg->bP; j++)
            if (gg->p[j]->is_fracture) {
                is_fracture = true;
                pp = gg->p[j];
                break;
            }
        if (!is_fracture) continue;

        if ((gg->Va + gg->Ve + gg->Vx) / gg->V0 > merge_factor) continue;


        int number_of_g =0;     //number of grain to be added to pp
        // looking for new grains to be added to pp
        for (int j = 0; j < gg->bN; j++) for(int jj=0; jj<gg->n[j]->bG;jj++){
            if (gg->n[j]->g[jj] != gg){
                Grain *g_tmp = gg->n[j]->g[jj];

                //checking if the grain is connected to the fracture
                bool is_connected_to_the_fracture = false;
                for(int k=0;k<g_tmp->bP;k++)
                    if(g_tmp->p[k]->is_fracture)
                        is_connected_to_the_fracture=true;

                //checking if the y position if fine
                bool y_location_ok = false;
                for (int k = 0; k < g_tmp->bN; k++)
                    if( g_tmp->n[k]->xy.y>=min(pp->n[0]->xy.y,pp->n[1]->xy.y) and
                        g_tmp->n[k]->xy.y<=max(pp->n[0]->xy.y,pp->n[1]->xy.y))
                        y_location_ok = true;


                if (!is_connected_to_the_fracture and y_location_ok and g_tmp->is_lhs == gg->is_lhs){
                    g_list.push_back(g_tmp);
                    pp->add_grain(g_tmp);
                    g_tmp->add_Pore(pp);
                }
            }
        }

        for (auto g_tmp : g_list){
            g_tmp->Va+= gg->Va/g_list.size();
            g_tmp->Ve+= gg->Ve/g_list.size();
            g_tmp->Vx+= gg->Vx/g_list.size();
        }

        if(g_list.empty()) {
            cerr<<"Problem with finding grains while merging for fracture..."<<endl;
            Va_tot-=gg->Va;
            Ve_tot-=gg->Ve;
            Vx_tot-=gg->Vx;
        }
        gg->Va=0; gg->Ve=0; gg->Vx=0;

    }

}


/**
* This function merge all empty grains.
* By empty grain we mean those with Va+Ve<=0. This condition can be change  by editing Network::to_be_merge().
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_empty_grains(){

	for(int i=0;i<NN;i++) n[i]->tmp = i;
	for(int i=0;i<NP;i++) p[i]->tmp = i;
	for(int i=0;i<NG;i++) g[i]->tmp = i;

	cerr<<"Merging empty grains..."<<endl;

	if(if_debugging_printing && if_verbose) debugging_for_merging ("Before merging one grain: s = " + to_string(tot_steps));

//	for(int i=NG-1;i>=0;i--) if(g[i]->to_be_merge(this)) { //condition for merging, can be adapted later
    for(int i=0;i<NG;i++) {
        if (g[i]->to_be_merge(this) > 0) { //condition for merging, can be adapted later


            if (g[i]->a < 0) break;

            cerr<<endl << "I am merging this guy: " << *g[i] << endl;

            if (g[i]->to_be_merge(this) == 2 and g[i]->bP > 1) {
                merge_the_other_pore_grain(g[i]);
            } else {
                if (g[i]->bP <= 1) merge_one_pore_grain(g[i]);
                else merge_one_grain(g[i]);
            }

            //additional checking of network structures
            if (if_debugging_printing && if_verbose)
                debugging_for_merging("After merging one grain: s = " + to_string(tot_steps));
        }
    }

	//Clear pathological pores: (probably unnecessary)
	for(int i=0; i<NP;i++) if(p[i]->n[0]==p[i]->n[1]) {
		Pore *p_tmp = p[i];
		cerr<<"Clearing pathological pores: "<<*p_tmp<<endl;
		p_tmp->remove_info_from_attached_nodes();
		p_tmp->remove_info_from_attached_grains();
		move_pore_to_the_end(p_tmp->tmp,NP);
	}

    //clear loops in the fracture
    for(int i=0; i<NP;i++) if(p[i]->n[0]->is_fracture and p[i]->n[1]->is_fracture and !p[i]->is_fracture ) {

        cerr<<"Entering the  loop clearance for pore"<<*p[i]<<endl;
        Pore* p_tmp=NULL;   //pore belonging to the loop
        for (int b=0;b<p[i]->n[0]->b;b++)
                if(p[i]->n[0]->n[b] == p[i]->n[1] and p[i]->n[0]->p[b]!=p[i] and p[i]->n[0]->p[b]->is_fracture)
                    p_tmp=p[i]->n[0]->p[b];
        if(p_tmp==NULL) continue;
        cerr<<"The other pore is "<<*p_tmp<<endl;

        Grain* g_tmp=NULL; // grain trapped in the loop
        for (int b=0;b<p[i]->bG;b++)
            if(p[i]->g[b]->bP==2 and p[i]->g[b]->bN==2)
                if((p[i]->g[b]->p[0]==p[i] and p[i]->g[b]->p[1]==p_tmp ) or p[i]->g[b]->p[1]==p[i] and p[i]->g[b]->p[0]==p_tmp)
                    g_tmp= p[i]->g[b];
        if(g_tmp==NULL) continue;
        cerr<<"The grain to be clear is "<<*g_tmp<<endl;

        merge_one_grain(g_tmp);

        cerr<<endl;
    }


    //look for floating fracture elements and mark them as not belonging to the fracture
    for(int i=0; i<NP;i++) if( p[i]->is_fracture and !(p[i]->n[0]->is_fracture and p[i]->n[1]->is_fracture)) {
        cerr<<"This pore has problematic fracture label: "<<*p[i]<<endl;
        //p[i]->is_fracture=false;
        //exit(123);
        }

    //check the fracture percolation
    cerr<<"Checking fracture percolation"<<endl;
    for(int i=0;i<NN;i++) n[i]->x=0;
    for(int i=0;i<NP;i++) p[i]->x=0;
    function<void(Node*)> fracture_percolation = [&](Node* n_tmp) -> void {
        n_tmp->x=1;
        if(!n_tmp->is_fracture) return;
        if (n_tmp->t == -1) return;
        for(int b=0;b<n_tmp->b;b++)
            if(n_tmp->p[b]->is_fracture and n_tmp->n[b]->is_fracture and n_tmp->p[b]->x==0){
                n_tmp->p[b]->x=1;
                fracture_percolation(n_tmp->n[b]);
            }
    };
    fracture_percolation(wi[0]);
    bool f_percol=false;
    for(int i=0;i<N_wo;i++)
        if(wo[i]->x==1)f_percol=true;
    if(f_percol)
        cerr<<"Fracture percolation is ok."<<endl;
    else
        cerr<<"No fracture percolation"<<endl;


    if(if_debugging_printing && if_verbose) debugging_for_merging ("After clearing pathological pores: s = " + to_string(tot_steps));

	//Clear pathological nodes: (probably unnecessary)
	for(int i=0; i<NN;i++) if(n[i]->b==0) {
		Node *n_tmp = n[i];
		cerr<<"Clearing pathological nodes: "<<*n_tmp<<endl;
		n_tmp->remove_form_IO_list(this);
		move_node_to_the_end(i,NN);
	}

	if(if_debugging_printing && if_verbose) debugging_for_merging ("After clearing pathological nodes: s = " + to_string(tot_steps));


	//Merging pathological grains
    if(if_verbose) cerr<<"Merging pathological grains..."<<endl;
	for(int i=0;i<NG;i++) if(g[i]->is_pathological()){
        if(g[i]->a<0) break;
		if(g[i]->bP<=1) merge_one_pore_grain(g[i]);
		else            merge_one_grain(g[i]);
		if(if_debugging_printing && if_verbose) debugging_for_merging ("After merging pathological grain: s = " + to_string(tot_steps));
	}
}

/**
* This function merge all pathological grains connected only to aone pore.
* Those grains can result from merging normal grain.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_one_pore_grain (Grain *gg){


	if(if_verbose) cerr<<endl<<endl<<"Merging one-pore-grain: "<<*gg<<endl<<"g_tmp = "<<gg->tmp<<endl;
	for(int i=0;i<gg->bN;i++)	gg->n[i]->remove_Grain(gg);
	for(int i=0;i<gg->bP;i++)   gg->p[i]->remove_grain(gg);

	//updating the total volume of A species
	if (gg->bP>0){
		if(gg->p[0]->bG>0)
            for(int i=0; i<p[0]->bG;i++){
                if(g[i]!=gg){
                    p[0]->g[i]->Va+=gg->Va/p[0]->bG;
                    p[0]->g[i]->Ve+=gg->Ve/p[0]->bG;
                    }
                }
		else 				{Va_tot -= gg->Va; Ve_tot -= gg->Ve;}
    }

	else          			{Va_tot -= gg->Va; Ve_tot -= gg->Ve;}
	move_grain_to_the_end(gg->tmp,NG);

}

void Network::merge_the_other_pore_grain(Grain *gg){

    if (if_verbose) {
        cerr << endl << endl << "Merging the other grain: " << *gg << endl;
        cerr << "Pores belonging to it: " << endl;
        for (int b = 0; b < gg->bP; b++)
            cerr << *(gg->p[b]) << endl;
        cerr << endl << "Nodes belonging to it: " << endl;
        for (int b = 0; b < gg->bN; b++)
            cerr << *(gg->n[b]) << endl;
        cerr << endl;
    }

    if(gg->bP!=3) {cerr<<"ERROR: Problem in merging the other grain bP!=3"<<endl; }
    Grain* g_goal = gg;
    for(int b=0;b<gg->bP;b++)
        if(!gg->p[b]->is_fracture)
            for(int bb = 0 ;bb<gg->p[b]->bG;bb++)
                if(gg->p[b]->g[bb] != gg) g_goal = gg->p[b]->g[bb];
    if (g_goal==gg) {cerr<<"ERROR: Problem in merging the other grain"<<endl; }

    int tmp = 0;
    for(int b=0;b<g_goal->bP;b++) if(g_goal->p[b]->is_fracture) {merge_one_grain(gg); return; }
    merge_one_grain(g_goal,gg);


}

/**
* This function merge one grain.
* In a result this grain will vanish and all connected pores will projected to one of them, the same with nodes.
*
* @param gg grain to be merged
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_one_grain(Grain *gg, Grain* gg_special) {


    if (if_verbose) {
        cerr << endl << endl << "Merging grain: " << *gg << endl;
        cerr << "Pores belonging to it: " << endl;
        for (int b = 0; b < gg->bP; b++)
            cerr << *(gg->p[b]) << endl;
        cerr << endl << "Nodes belonging to it: " << endl;
        for (int b = 0; b < gg->bN; b++)
            cerr << *(gg->n[b]) << endl;
        cerr << endl;
    }

    //choose main pore
    Pore *p1 = NULL;  //master pore (the one that will remain)

    int nr_of_f =0;
    for (int b = 0; b < gg->bP; b++)
        if (gg->p[b]->is_fracture){
            p1 = gg->p[b];
            nr_of_f++;
        }

    if(nr_of_f>1){
        double q_max = 0;
        for (int b = 0; b < gg->bP; b++)
            if (gg->p[b]->is_fracture && fabs(gg->p[b]->q) >= q_max) {
                p1 = gg->p[b];
                q_max = fabs(gg->p[b]->q);
            }
    }

    if (p1 == NULL) {
        double q_max = 0;
        for (int b = 0; b < gg->bP; b++)
            if (fabs(gg->p[b]->q) >= q_max) {
                p1 = gg->p[b];
                q_max = fabs(gg->p[b]->q);
            }
    }
    if (if_verbose)cerr << "Master pore is " << *p1 << endl;
    gg->set_effective_d_and_l(p1, this);  //calculate d and l basing on S and R of the entire pore


    //temporary list of grains
    int b_tmp = 0;
    for (int b = 0; b < gg->bP; b++) { b_tmp += gg->p[b]->bG; }
    Grain **g_tmp = new Grain *[b_tmp];
    b_tmp = 0;
    for (int b = 0; b < p1->bG; b++) if (p1->g[b] != gg) { g_tmp[b_tmp++] = p1->g[b]; } //copy list of pores and nodes

    //clearing unnecessary pores
    Grain gg_copy = *gg;
    for (int b = 0; b < gg_copy.bP; b++)
        if (gg_copy.p[b] != p1) {
            Pore *p_tmp = gg_copy.p[b];
            //checking if new grain will be attached to p1 (if p_tmp nodes are not projected to the same p1 node)
            Node *n1_tmp = p1->find_closest_node(p_tmp->n[0]);
            Node *n2_tmp = p1->find_closest_node(p_tmp->n[1]);
            //adding new grains to the list of p1 grains
            if (n1_tmp != n2_tmp)
                for (int bb = 0; bb < p_tmp->bG; bb++) {
                    bool if_new = true;
                    if (p_tmp->g[bb] == gg) if_new = false;
                    for (int j = 0; j < b_tmp; j++) if (p_tmp->g[bb] == g_tmp[j]) if_new = false;
                    if (if_new) {
                        g_tmp[b_tmp++] = p_tmp->g[bb];
                        p_tmp->g[bb]->change_pores(p_tmp, p1);
                    }
                }
            //adding info about flow to master pore (for better printing only)
            p1->q = p1->q + _sign(p1->q) * abs(p_tmp->q) / 3;

            //clearing unused pore
            if (if_verbose) cerr << "Clearing pore: " << *p_tmp << endl;
            p_tmp->remove_info_from_attached_nodes();
            p_tmp->remove_info_from_attached_grains();
            move_pore_to_the_end(p_tmp->tmp, NP);
        }

    p1->bG = b_tmp;
    delete[] p1->g;
    p1->g = g_tmp;


    //merging all nodes but two (p1->n[0] and p1->n[1])
    gg_copy = *gg;
    for (int b = 0; b < gg_copy.bN; b++)
        if (gg_copy.n[b] != p1->n[0] && gg_copy.n[b] != p1->n[1]) {
            Node *n_master = p1->find_closest_node(gg_copy.n[b]);    //closer p1 node
            merge_nodes(n_master, gg_copy.n[b]);
        }

    //changing info about grains in master nodes
    p1->n[0]->remove_Grain(gg);
    p1->n[1]->remove_Grain(gg);

    int ng_tmp=0;
    for (int i = 0; i < p1->bG; i++) if(p1->g[i]->Va>0 && gg->is_lhs == p1->g[i]->is_lhs) ng_tmp++;

    if(gg_special!=NULL && gg_special!=gg) {
        gg_special->Va += gg->Va;
        gg_special->Ve += gg->Ve;
    }
    else if (ng_tmp>0) {
        for (int i = 0; i < p1->bG; i++) if (p1->g[i]->Va>0 && gg->is_lhs == p1->g[i]->is_lhs) {
            p1->g[i]->Va += gg->Va / ng_tmp;
            p1->g[i]->Ve += gg->Ve / ng_tmp;
            }
    }
    else{
        Va_tot -= gg->Va;
        Ve_tot -= gg->Ve;
    }

	move_grain_to_the_end(gg->tmp,NG);
    if(if_verbose){
        cerr<<"After moving to the end " << *gg<<endl<<endl<<endl;
    }


}


/**
* This function do special debugging for merging
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::debugging_for_merging(string text){

	check_network_connections();

	for(int i=0;i<NN;i++) n[i]->tmp = n[i]->a;
	for(int i=0;i<NP;i++) p[i]->tmp = p[i]->a;
	for(int i=0;i<NG;i++) g[i]->tmp = g[i]->Va;
	description_note = text;
	net_ps<<*this;
	for(int i=0;i<NN;i++) n[i]->tmp = i;
	for(int i=0;i<NP;i++) p[i]->tmp = i;
	for(int i=0;i<NG;i++) g[i]->tmp = i;

}

/**
* This function merge two nodes connected to the grain that will be deleted..
*
* @param n1 master node
* @param n2 node to be deleted
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::merge_nodes(Node *n1, Node *n2){

	if(if_verbose) cerr<<"Merging nodes: "<<endl<<*n1<<endl<<*n2<<endl;
	if(n1->a<0 || n2->a<0) cerr<<"ERROR in merging nodes!!!!"<<endl;


	//Updating list of neighbors
	Node  ** n_tmp = new Node  *[n1->b +n2->b];
	Pore  ** p_tmp = new Pore  *[n1->b +n2->b];
	Grain ** g_tmp = new Grain *[n1->bG+n2->bG];
	int b_tmp=0; int bg_tmp=0;
    for(int b=0;b<n1->b;b++)  if(n1->n[b]!=n2){ n_tmp[b_tmp]=n1->n[b];  p_tmp[b_tmp++]=n1->p[b];} //copy list of pores and nodes
	Node n2_copy = *n2;
    for(int b=0;b<n2_copy.b;b++)
		if(n2_copy.n[b]!=n1){
			n_tmp[b_tmp]   = n2_copy.n[b];
			p_tmp[b_tmp++] = n2_copy.p[b];
			n2_copy.n[b] -> change_node_neighbours(n2,n1,n2_copy.p[b]);
			n2_copy.p[b] -> change_pore_neighbours(n2,n1);
			}
		else{ // To not leave unattached pores
			Pore *p_del = n2_copy.p[b];
			if(if_verbose) cerr<<"Clearing pore: "<<*p_del<<endl;
			//p_del->remove_info_from_attached_nodes();  //
			p_del->remove_info_from_attached_grains();
			move_pore_to_the_end(p_del->tmp,NP);
		}

	for(int b=0;b<n1->bG;b++) g_tmp[bg_tmp++]=n1->g[b];   //copy list of grains
	for(int b=0;b<n2->bG;b++) {
		bool if_new=true;
		for(int j=0;j<bg_tmp;j++) if(g_tmp[j]==n2->g[b]) if_new = false;
		if(if_new) g_tmp[bg_tmp++]=n2->g[b];
		n2->g[b]->change_nodes(n2,n1);
	}

	//update neighbors for n1
	n1->b = b_tmp; n1->bG = bg_tmp;
	delete [] n1->n; n1->n = n_tmp;
	delete [] n1->p; n1->p = p_tmp;
	delete [] n1->g; n1->g = g_tmp;

	//change position of node n1 (if the distance is not to large - see periodic boundary conditions)
	if(n1->xy - n2->xy<N_x*2./3.){
        if(!n1->is_fracture and !n2->is_fracture) {
            double q1 = n1->total_pores_q();
            double q2 = n2->total_pores_q();
            if (q1 + q2 == 0) {
                q1 = 1;
                q2 = 1;
            }
            n1->xy = ((2 * q1 / (q1 + q2)) * n1->xy) * ((2 * q2 / (q1 + q2)) * n2->xy);
        }

        if(n2->is_fracture )
            n1->xy = n2->xy;

	}


    if(n2->is_fracture) n1->is_fracture = true;
	//change type of master node if necessary
	if(n2->t!=0){
		if(n1->t==0){
			n1->t = n2->t;
			if(n2->t==1)  for(int i=0; i<N_wi; i++) if(wi[i]==n2) {wi[i] = n1; break;}
			if(n2->t==-1) for(int i=0; i<N_wo; i++) if(wo[i]==n2) {wo[i] = n1; break;}
		}
		else  n2->remove_form_IO_list(this); //remove n2 from wi/wo list
	}


	//clear node n2
	move_node_to_the_end(n2->tmp,NN);

	if(if_verbose) cerr<<"After merging nodes: "<<endl<<*n1<<endl<<*n2<<endl<<flush;
}


