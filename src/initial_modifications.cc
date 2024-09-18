#include "network.h"
#include "printing.h"

void Network :: create_a_fracture(double factor, Node *n_1, Node *n_2) {


    cerr<<"Creating a single - layer fracture"<<endl;

    //1. Optionally finding the inlet at the center of the system:
    if (n_1 == nullptr){
        n_1 = wi[0];
        Point center = {N_x/2.,0.0};
        double dist_min = n_1->xy - center;
        for (int i=0;i<N_wi;i++){

            if(wi[i]->xy - center < dist_min){
                n_1 = wi[i];
                dist_min = wi[i]->xy - center;
            }
        }
    }

    //2. Optionally finding the outlet at the center of the system:
    if (n_2 == nullptr){
        n_2 = wo[0];
        Point center = {N_x/2.,N_y/1.};
        double dist_min = n_2->xy - center;
        for (int i=0;i<N_wo;i++){
            if(wo[i]->xy - center < dist_min){
                n_2 = wo[i];
                dist_min = wo[i]->xy - center;
            }
        }
    }

    find_shortest_path(n_1,n_2);

    for (int i=0;i<NP;i++)
        if (p[i]->tmp == 1){
            p[i]->d = p[i]->d * factor;   //the initial diameter is larger for a fracture
            p[i]->is_fracture = true;              //the fracture can be implemented in different way that the regular pore
        }
    for (int i=0;i<NN;i++)
        if (n[i]->tmp == 1)
            n[i]->is_fracture=true;


    // Moving node positions to visualize a fracture better
    //Preparing network

    int i=0;
    while (n[i]->t!=0 || n[i]->xy.x>N_x/3. || n[i]->is_fracture) i++;
    find_R_half(n[i]);
    for (int i=0;i<NG;i++){
        int tmp_lhs=0;
        for(int b=0;b<g[i]->bN;b++)
            if(g[i]->n[b]->is_LHS) tmp_lhs++;
        if(tmp_lhs>0)
            g[i]->is_lhs=true;
    }


    if(point_inlet) {
        cerr<<"Setting point inlet..."<<endl;
        //setting new outlet for point outlets
        for (int i = 0; i < N_wi; i++)  wi[i]->t =0;
        delete[] wi;
        N_wi = 1;
        wi = new Node *[N_wi];
        wi[0] = n_1;
        n_1->t = 1;
    }

    if(point_outlet) {
        cerr<<"Setting point outlet..."<<endl;
        //setting new outlet for point outlets
        for (int i = 0; i < N_wo; i++)  wo[i]->t =0;

        delete[] wo;
        N_wo = 1;
        wo = new Node *[N_wo];
        wo[0] = n_2;
        n_2->t = -1;

    }

    }

void Network::find_R_half(Node *n0) {
        n0->is_LHS=1;
    for (int i=0;i<n0->b;i++)
        if(n0->n[i]->is_LHS==0 and n0->p[i]->d!=0 and !n0->n[i]->is_fracture and n0->n[i]->xy-n0->xy<N_x/3. ){
            find_R_half(n0->n[i]);
        }



}


/**
* This function creates initial cut in the network
* - the region connected to the inlet of the system with different permeability that the rest of the system.
*
* @param cut_w width of the cut
* @param cut_l length of the cut
* @param factor factor for recycling d in the cut (by default eq 1)
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::create_an_inlet_cut(int cut_w, int cut_l, double factor){

	cerr<<"Adding an inlet cut..."<<endl;

    /// in this scenario the single-layer fracture form the middle of the system to the nearest outlet will be created;
    if(cut_w == 0) return create_a_fracture(factor);

	//if(cut_l >=N_y) {cerr<<"WARNING : Cut length is too large. No cut is made."<<endl; return;}

    // clearing info about inlet and outlet
    if(add_well || point_inlet){
        cerr<<"Changing inlets/outlets (well mode)..."<<endl;
        //clearing old inlets
        for (int i = 0; i < NN; i++) if (n[i]->t == 1) n[i]->t = 0;

        //setting new inlet
        N_wi = 0;
        delete[] wi;
        for (int i = 0; i < NN; i++)
            if (n[i]->xy.y < 3 && n[i]->xy.x >= (N_x / 2. - cut_w / 2.) &&
                n[i]->xy.x < (N_x / 2. + cut_w / 2.))  //check if 3 is ok in this equation
            {
                n[i]->t = 1;
                N_wi++;
                cerr << "Dodaje do inletow " << *n[i] << "  N_wi = " << N_wi << endl;
            }
        wi = new Node *[N_wi];
        int j = 0;
        for (int i = 0; i < NN; i++) if (n[i]->t == 1) wi[j++] = n[i];
        if (j != N_wi)
            cerr << "WARNING: Problem with setting new inlets: N_wi = " << N_wi << "   j = " << j << endl;

    }

    if (add_well) {//setting new outlet for point outlets
        //adding new outlet
        delete[] wo;
        for (int i = 0; i < NN; i++)
            if (n[i]->xy.y < 3 && n[i]->xy.x >= (-cut_w / 2.) && n[i]->xy.x < (+cut_w / 2.)) {
                n[i]->t = -1;
                N_wo++;
                cerr << "Dodaje do outletow " << *n[i] << "  N_wo = " << N_wo << endl;
            }
        wo = new Node *[N_wo];
        int j = 0;
        for (int i = 0; i < NN; i++) if (n[i]->t == -1) wo[j++] = n[i];
        if (j != N_wo)
            cerr << "WARNING: Problem with setting new outlets: N_wo = " << N_wo << "   j = " << j << endl;
    }

    if(point_outlet) {
        //setting new outlet for point outlets
        for (int i = 0; i < NN; i++) if (n[i]->t == -1) n[i]->t = 0;
        N_wo = 0;
        delete[] wo;
        for (int i = 0; i < NN; i++)
            if (n[i]->xy.y >N_y - 4 &&
                n[i]->xy.x >= (N_x / 2. - cut_w / 2.) &&
                n[i]->xy.x < (N_x / 2. + cut_w / 2.))  //check if 3 is ok in this equation
            {
                n[i]->t = -1;
                N_wo++;
                cerr << "Dodaje do outletow " << *n[i] << "  N_wi = " << N_wi << endl;
            }

        wo = new Node *[N_wo];
        int j = 0;
        for (int i = 0; i < NN; i++) if (n[i]->t == -1) wo[j++] = n[i];
        if (j != N_wo)
                cerr << "WARNING: Problem with setting new outlets: N_wo = " << N_wo << "   j = " << j << endl;

    }

//    cerr<< "Wypisuje wszystkie nody wejsciowe: "<<endl;
//    for(int i=0;i<N_wi;i++) cerr<<*wi[i]<<endl;
//
//    cerr<< "Wypisuje wszystkie nody wyjsciowe: "<<endl;
//    for(int i=0;i<N_wo;i++) cerr<<*wo[i]<<endl;


		for(int i=0;i<NN;i++) {
            if (n[i]->xy.y < cut_l && n[i]->xy.x >= (N_x / 2. - cut_w / 2.) && n[i]->xy.x < (N_x / 2. + cut_w / 2.))
                for (int j = 0; j < n[i]->b; j++)
                    if (n[i]->p[j]->d > 0){
                        n[i]->p[j]->d = d0 * factor;
                        if(!if_reactions_in_the_fracture) n[i]->p[j]->is_active = false;
                    }
        }

        if(add_well) {
            for(int i=0;i<NN;i++){
                if(n[i]->xy.y <cut_l && n[i]->xy.x >= (- cut_w/2.) && n[i]->xy.x< (+ cut_w/2.))
                    for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0) {
                            n[i]->p[j]->d = d0 * factor;
                            if(!if_reactions_in_the_fracture) n[i]->p[j]->is_active = false;
                        }
            }
        }

}


void Network::create_tilted_fracture(int cut_w, double factor){

        cerr<<"Adding a tilted fracture..."<<endl;

        //if(cut_l >=N_y) {cerr<<"WARNING : Cut length is to large. No cut is made."<<endl; return;}

        // clearing info about inlet and outlet
        if( point_inlet){
            cerr<<"Changing inlets/outlets (well mode)..."<<endl;
            //clearing old inlets
            for (int i = 0; i < NN; i++) if (n[i]->t == 1) n[i]->t = 0;

            //setting new inlet
            N_wi = 0;
            delete[] wi;
            for (int i = 0; i < NN; i++)
                if (n[i]->xy.y < 3 && n[i]->xy.x >= (N_x*2./3. - cut_w / 2.) &&
                    n[i]->xy.x < (N_x*2./3. + cut_w / 2.))  //check if 3 is ok in this equation
                {
                    n[i]->t = 1;
                    N_wi++;
                    cerr << "Dodaje do inletow " << *n[i] << "  N_wi = " << N_wi << endl;
                }
            wi = new Node *[N_wi];
            int j = 0;
            for (int i = 0; i < NN; i++) if (n[i]->t == 1) wi[j++] = n[i];
            if (j != N_wi)
                cerr << "WARNING: Problem with setting new inlets: N_wi = " << N_wi << "   j = " << j << endl;

        }


        if(point_outlet) {
            //setting new outlet for point outlets
            for (int i = 0; i < NN; i++) if (n[i]->t == -1) n[i]->t = 0;
            N_wo = 0;
            delete[] wo;
            for (int i = 0; i < NN; i++)
                if (n[i]->xy.y >N_y - 4 &&
                    n[i]->xy.x >= (N_x*2./3. - cut_w/2.) &&
                    n[i]->xy.x < (N_x*2./3. + cut_w/2.))  //check if 3 is ok in this equation
                {
                    n[i]->t = -1;
                    N_wo++;
                    cerr << "Dodaje do outletow " << *n[i] << "  N_wi = " << N_wi << endl;
                }

            wo = new Node *[N_wo];
            int j = 0;
            for (int i = 0; i < NN; i++) if (n[i]->t == -1) wo[j++] = n[i];
            if (j != N_wo)
                cerr << "WARNING: Problem with setting new outlets: N_wo = " << N_wo << "   j = " << j << endl;

        }

//    cerr<< "Wypisuje wszystkie nody wejsciowe: "<<endl;
//    for(int i=0;i<N_wi;i++) cerr<<*wi[i]<<endl;
//
//    cerr<< "Wypisuje wszystkie nody wyjsciowe: "<<endl;
//    for(int i=0;i<N_wo;i++) cerr<<*wo[i]<<endl;


        for(int i=0;i<NN;i++){
            if( n[i]->xy.x >= (4./3.*N_x/N_y*(N_y/2. - n[i]->xy.y) - cut_w/2.) &&  n[i]->xy.x< (4./3.*N_x/N_y*(N_y/2. - n[i]->xy.y)  + cut_w/2.))
                for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0){
                    n[i]->p[j]->d = d0*factor;
                    if(!if_reactions_in_the_fracture) n[i]->p[j]->is_active = false;
                    }

            if( n[i]->xy.x >= 4./3.*N_x/N_y*(n[i]->xy.y - N_y/2.) - cut_w/2. && n[i]->xy.x< (4./3.*N_x/N_y*(n[i]->xy.y - N_y/2.)) + cut_w/2.)
                for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0){
                    n[i]->p[j]->d = d0*factor;
                    if(!if_reactions_in_the_fracture) n[i]->p[j]->is_active = false;
                    }
        }

}


/**
* This function creates an initial pattern.
* It is useful when debugging.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::create_an_initial_pattern(){

//	cerr<<"Adding an initial pattern..."<<endl;
//
//	for (int i=0; i<NN; i++) if( i==7 || i == 13) p[2*i+0]->d=3*d0;
//	p[2*8+1]->d=3*d0;
//	p[2*7+1]->d=0.1*d0;

}
