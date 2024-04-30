#include "network.h"
#include "printing.h"

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

	if(cut_l >=N_y) {cerr<<"WARNING : Cut length is to large. No cut is made."<<endl; return;}

    // clearing ingo about inlet and outlet
    if(add_well || point_inlet){
        cerr<<"Changing inlets/outlets (well mode)..."<<endl;
        if(type_of_topology != "hexagonal") {
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

            if (add_well) {
                wi = new Node *[N_wi];
                int j = 0;
                for (int i = 0; i < NN; i++) if (n[i]->t == 1) wi[j++] = n[i];
                if (j != N_wi)
                    cerr << "WARNING: Problem with setting new inlets: N_wi = ." << N_wi << "   j = " << j << endl;

                //adding new outlet
                delete[] wo;
                for (int i = 0; i < NN; i++)
                    if (n[i]->xy.y < 3 && n[i]->xy.x >= (-cut_w / 2.) && n[i]->xy.x < (+cut_w / 2.)) {
                        n[i]->t = -1;
                        N_wo++;
                        cerr << "Dodaje do outletow " << *n[i] << "  N_wo = " << N_wo << endl;
                    }
                wo = new Node *[N_wo];
                j = 0;
                for (int i = 0; i < NN; i++) if (n[i]->t == -1) wo[j++] = n[i];
                if (j != N_wo)
                    cerr << "WARNING: Problem with setting new outlets: N_wo = ." << N_wo << "   j = " << j << endl;
            }
            else{
                //setting new outlet for point outlets
                N_wo = 0;
                delete[] wo;
                for (int i = 0; i < NN; i++)
                    if (n[i]->xy.y < 3 && n[i]->xy.x >= (N_x / 2. - cut_w / 2.) &&
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
                        cerr << "WARNING: Problem with setting new outlets: N_wo = ." << N_wo << "   j = " << j << endl;

            }
        }

       else  {

            cerr<<"Adding well and changing inlets/outlets for hexagonals."<<endl;
            //hexagonal case
            for (int i = 0; i < NN; i++) n[i]->t = 0;
            delete[] wi;
            N_wi = 1;
            wi = new Node *[N_wi];
            wi[0] = n[N_x / 2];
            wi[0]->t = 1;
//        cerr<<"\n Wybralam sobie wejscie:   "<< *n[N_x/2]<<endl<<endl;


            auto wo_new = new Node *[N_wo + 1];
            for (int i = 0; i < N_wo; i++) wo_new[i] = wo[i];
            wo_new[N_wo] = n[0];
            N_wo = N_wo + 1;
            delete[] wo;
            wo = wo_new;
            for (int i = 0; i < N_wo; i++) wo[i]->t = -1;
        }
    }

//
//    cerr<< "Wypisuje wszystkie nody wejsciowe: "<<endl;
//    for(int i=0;i<N_wi;i++) cerr<<*wi[i]<<endl;
//
//    cerr<< "Wypisuje wszystkie nody wyjsciowe: "<<endl;
//    for(int i=0;i<N_wo;i++) cerr<<*wo[i]<<endl;


    if(type_of_topology == "hexagonal"){
		for(int i=0; i<cut_l*N_x;i++){
			if(i%N_x >= (N_x/2. - cut_w/2.)  && i%N_x < (N_x/2. + cut_w/2.))
				for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0) {
                    n[i]->p[j]->d = d0*factor;}
        }

        if(add_well) {  //adding additional cut at x=0=Nx
            for(int i=0; i<cut_l*N_x;i++){
                if(i%N_x >= - cut_w/2.  && i%N_x < cut_w/2.)
                    for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0) n[i]->p[j]->d = d0*factor;
            }
        }
	}

	else{
		for(int i=0;i<NN;i++){
			if(n[i]->xy.y <cut_l && n[i]->xy.x >= (N_x/2. - cut_w/2.) && n[i]->xy.x< (N_x/2. + cut_w/2.))
				for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0) n[i]->p[j]->d = d0*factor;}

        if(add_well) {
            for(int i=0;i<NN;i++){
                if(n[i]->xy.y <cut_l && n[i]->xy.x >= (- cut_w/2.) && n[i]->xy.x< (+ cut_w/2.))
                    for(int j=0;j<n[i]->b;j++) if(n[i]->p[j]->d>0) n[i]->p[j]->d = d0*factor;}
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
