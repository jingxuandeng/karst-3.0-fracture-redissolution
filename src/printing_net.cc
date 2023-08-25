#include "network.h"

double z_to_print=0;
double max_distance = 4;




ofstream_ps & operator << (ofstream_ps & stream, Grain &g){
	if(g.bN!=3)  return stream;
	if(g.n[0]->xy-g.n[1]->xy<max_distance && g.n[1]->xy-g.n[2]->xy<max_distance && g.n[2]->xy-g.n[0]->xy<max_distance)

		//stream<<Trojkacik(g.n[0]->xy,g.n[1]->xy,g.n[2]->xy,g.tmp,Kolor(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX))<<endl;
		//stream<<Trojkacik(g.n[0]->xy,g.n[1]->xy,g.n[2]->xy,g.tmp,Kolor(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX))<<endl;
	//else
		stream<<Trojkacik(g.n[0]->xy,g.n[1]->xy,g.n[2]->xy,g.tmp,Kolor(0.9,0.9,0.1))<<endl;

	return stream;}



ofstream_ps & print_grain_with_scaling (ofstream_ps & stream, Grain &g, Network &S){
	if(g.bN<3) {
		if(S.if_verbose) cerr<<"Patological grain (bN = "<<g.bN<<", bP = "<<g.bP<<") has to be implemented."<<endl<<g<<endl;
		return stream; }
	if(g.Va + g.Ve <=0)  return stream;

	int b = g.bN;
	Point *PP = S.initial_xy[g.a];
	for (int i=0;i<b;i++) if(PP[i]-PP[(i+1)%b]>max_distance) return stream;

	Point *p  = new Point [b];    //initial node position
	Point *pp = new Point [b];    // rescaled position
	Point p0;   //middle of a agrain
	for (int i=0;i<b;i++) {p[i] = PP[i];  p0 = p0+p[i];}
	p0 = (1./b)*p0;
	Point p00 = (-1.)*p0;
	double factr = (g.Va + g.Ve + g.Vx)/g.calculate_maximal_volume(&S);
	if (factr<0)   factr = 0;
	if (factr>1)   factr = 1;
	if (!(factr>=0)) factr=0;

	for (int i=0;i<b;i++){
		pp[i]  = p[i] + p00;
		pp[i]  = factr*pp[i];
		pp[i]  = pp[i] + p0;
	}

	g.tmp=666; //g.Ve/(g.Ve+ g.Va); // for not printing the grain labels

	Kolor color =  Kolor(0.33,0.33,0.33);
	//fancy color for precipitation
	double factor = 2;
	if (S.if_precipitation)  color = Kolor( factor*g.Ve/(g.Va+factor*g.Ve),0,g.Va/(g.Va+factor*g.Ve));

	//if(g.bN==3)  stream<<Trojkacik(pp[0],pp[1],pp[2],g.tmp,color)<<endl<<flush;
	stream<<Wielobok (g.bN, pp,g.tmp,color)<<endl<<flush;


	delete [] p;
	delete [] pp;

	return stream;
}

ofstream_ps & operator << (ofstream_ps & stream, Node &n){
	//stream<<Kropka(n.xy,n.tmp,Kolor(1,0,0),0.1);
	if(n.xy.z != z_to_print) return stream;
	stream<<Kropa(n.xy,&n,Kolor(1,0,0),0.1);
	return stream;}

ofstream_ps & operator << (ofstream_ps & stream, Pore &p){

	Kolor kkk(0.5,0.5,0.5);

	if(p.x==1) kkk=Kolor(0.2,0.8,0.2);
	if(p.x==2) kkk=Kolor(0.8,0.2,0.2);

	bool if_debug=true;
	if(p.n[0]->xy - p.n[1]->xy < max_distance && p.d<300&& p.n[0]->xy.z == z_to_print && p.n[1]->xy.z == z_to_print){
		if (if_debug && p.x==2) stream<<Porek(p.n[0]->xy,p.n[1]->xy,p.d/10,p.tmp,kkk);
		else                    stream<<Porek(p.n[0]->xy,p.n[1]->xy,p.d/4 ,p.tmp,kkk);
		}

	return stream;}



ofstream_ps & operator << (ofstream_ps & stream, Network &S){

	if      (S.printing_mode == "debugging")   		Print_network_in_debugging_style   (stream,S);
	else if (S.printing_mode == "dissolution") 		Print_network_in_dissolution_style (stream,S);
	else if (S.printing_mode == "grains")      		Print_network_in_grain_style       (stream,S);
	else if (S.printing_mode == "grains_and_diss"){
		Print_network_in_grain_style       (stream,S);
		Print_network_in_dissolution_style (stream,S);
	}

	else
		cerr << "WARNING: Incorrect type of printing mode."<<endl;
	return stream;
}


void Print_network_in_debugging_style (ofstream_ps & stream, Network &S){

	if(S.if_save_ps==0){ cerr<<"Error during printing ps in debugging style."<<endl; return ;}


	cerr<<"Printing network for debugging ..."<<endl;

	int N=(S.N_x);
	int M=(S.N_y);
	double skala  = 400./(S.L_out*max(N,M));     //
	double x_zero = 100./skala;
	double y_zero = 750./skala;


	if(S.pages_saved==0) stream << "%!PS-Adobe-3.0" << endl<<"%%Pages:"<<S.pages_tot<<endl<<endl;
	
	stream <<"%%Page: "<<S.pages_saved+1<<" "<<S.pages_saved+1<<endl<<endl;
	
	//stream << "%%Orientation: Landscape" << endl;
	//stream << "%%DocumentMedia: a4 595 842 80 () ()" << endl;
	//stream<<"-1 -1 scale"<<endl;
	stream <<"1 setlinejoin 1 setlinecap 0.02 setlinewidth"<<endl;

	stream<<skala<<"\t"<<skala<<" scale"<<endl;
	stream<<x_zero<<"\t"<<y_zero<<" translate"<<endl;
	stream<<"1 setlinejoin"<<endl;
	stream<<"1 setlinecap" <<endl;
	//stream<<"1 -1 scale   0 -300 transform"<<endl;
	//stream<<"1 -1 scale"<<endl;

	//stream<<" 0 0 moveto 0 300 lineto stroke"<<endl;

	//title
	stream<<Kolor(0,0,0);
	stream<<"/Times-Bold findfont "<<20./skala<<" scalefont setfont"<<endl;
	stream<<"0 "<<-450./skala<<" moveto"<<endl;
	stream<<"0 0 ("<<S.description_note<<") ashow stroke"<<endl<<endl;


	//for(int i=0;i<S.NG;i++) stream<<*S.g[i];//	cerr<<"Printing grain: "<<*S.g[i]<<endl;}
	//for(int i=0;i<S.NG;i++) print_grain_with_scaling(stream,*(S.g[i]),S);
	for(int i=0;i<S.NP;i++) stream<<*S.p[i];// 	cerr<<"Printing pore: "<<*S.p[i]<<endl;}
	for(int i=0;i<S.NN;i++) stream<<*S.n[i];//  cerr<<"Printing node: "<<*S.n[i]<<endl;}

	stream << "showpage "<<endl<<flush;

}

void Print_network_in_dissolution_style (ofstream_ps & stream, Network &S){

	if(S.if_save_ps==0){ cerr<<"Error during printing ps in dissolution style."<<endl; return ;}

	cerr<<"Printing network..."<<endl;

	if(S.print_diss_factor) S.check_diss_pattern(S.print_diss_factor);

	int N=(S.N_x);
	int M=(S.N_y);
	double skala = 400./(S.L_out*max(N,M));     //
	double x_zero = 100./skala;
	double y_zero = 750./skala;


	if(S.pages_saved==0) stream << "%!PS-Adobe-3.0" << endl<<"%%Pages:"<<S.pages_tot<<endl<<endl;

	stream <<"%%Page: "<<S.pages_saved+1<<" "<<S.pages_saved+1<<endl<<endl;
	stream <<"1 setlinejoin 1 setlinecap 0.02 setlinewidth"<<endl;

	stream<<skala<<"\t"<<skala<<" scale"<<endl;
	stream<<x_zero<<"\t"<<y_zero<<" translate"<<endl;
	stream<<"1 setlinejoin"<<endl;
	stream<<"1 setlinecap" <<endl;

	//title
	stream<<Kolor(0,0,0);
	stream<<"/Times-Bold findfont "<< 20./skala<<" scalefont setfont"<<endl;
	stream<<"0 "<<-450./skala<<" moveto"<<endl;
	stream<<"0 0 ("<<S.description_note<<") ashow stroke"<<endl<<endl;


	if(S.if_precipitation)         for(int i=0;i<S.NP;i++)                  {
		Kolor kkk(0.5,0.5,0.5);
		Pore &p = *S.p[i];
		double r=0.5, g=0.5, b=0.5;
		if (p.d<S.d0 and p.d!=0){
			r = 1;
			g = ((p.d) - (S.d_min))/((S.d0) - (S.d_min));
			b = 0;
			kkk=Kolor(r,g,b);}
		if(p.x==1) kkk=Kolor(0.0,0.0,0.0);

		p.tmp=666;  // no printing names
		if(p.n[0]->xy - p.n[1]->xy < max_distance && p.d<300&& p.n[0]->xy.z == z_to_print && p.n[1]->xy.z == z_to_print){
			if (p.x == 1)           stream<<Porek(p.n[0]->xy,p.n[1]->xy,p.d ,p.tmp,kkk);
			else                    stream<<Porek(p.n[0]->xy,p.n[1]->xy,p.d ,p.tmp,kkk);
			}
	}
	else if( S.print_diss_factor)  for(int i=0;i<S.NP;i++) if(S.p[i]->x==1) {S.p[i]->tmp=666; stream<<*S.p[i];}
	else if(!S.print_diss_factor)  for(int i=0;i<S.NP;i++)                  {S.p[i]->tmp=666; stream<<*S.p[i];}

	stream << "showpage "<<endl<<flush;
}

void Print_network_in_debugging_style_tmp (ofstream_ps & stream, Network &S,int nr_tmp){

	if(S.if_save_ps==0){ cerr<<"Error during printing ps in debugging style."<<endl; return ;}


	cerr<<"Printing network for debugging ..."<<endl;

	int N=(S.N_x);
	int M=(S.N_y);
	double skala  = 400./(S.L_out*max(N,M));     //
	double x_zero = 100./skala;
	double y_zero = 750./skala;


	if(S.pages_saved==0) stream << "%!PS-Adobe-3.0" << endl<<"%%Pages:"<<S.pages_tot<<endl<<endl;

	stream <<"%%Page: "<<S.pages_saved+1<<" "<<S.pages_saved+1<<endl<<endl;

	//stream << "%%Orientation: Landscape" << endl;
	//stream << "%%DocumentMedia: a4 595 842 80 () ()" << endl;
	//stream<<"-1 -1 scale"<<endl;
	stream <<"1 setlinejoin 1 setlinecap 0.02 setlinewidth"<<endl;

	stream<<skala<<"\t"<<skala<<" scale"<<endl;
	stream<<x_zero<<"\t"<<y_zero<<" translate"<<endl;
	stream<<"1 setlinejoin"<<endl;
	stream<<"1 setlinecap" <<endl;
	//stream<<"1 -1 scale   0 -300 transform"<<endl;
	//stream<<"1 -1 scale"<<endl;

	//title
	stream<<Kolor(0,0,0);
	stream<<"/Times-Bold findfont "<<20./skala<<" scalefont setfont"<<endl;
	stream<<"0 "<<-450./skala<<" moveto"<<endl;
	stream<<"0 0 ("<<S.description_note<<") ashow stroke"<<endl<<endl;


	stream<<Kolor(1,0,1);
	stream<<" 0 0 moveto 0 -100 lineto 100 -100 lineto 100 0 lineto 0 0 lineto stroke"<<endl;
	stream<<Kolor(0,0,0);


	//for(int i=0;i<S.NG;i++) stream<<*S.g[i];//	cerr<<"Printing grain: "<<*S.g[i]<<endl;}
	Node *nn=S.n[nr_tmp];
	stream<<*nn;
	for(int bb=0;bb<nn->b;bb++) stream<<*(nn->p[bb]);
	for(int bb=0;bb<nn->b;bb++)	stream<<*(nn->n[bb]);
	stream<<*nn;

	//for(int i=0;i<S.NP;i++) stream<<*S.p[i];// 	cerr<<"Printing pore: "<<*S.p[i]<<endl;}
	//for(int i=0;i<S.NN;i++) stream<<*S.n[i];//  cerr<<"Printing node: "<<*S.n[i]<<endl;}

	stream << "showpage "<<endl<<flush;
}


void Print_network_in_grain_style (ofstream_ps & stream, Network &S){

	if(S.if_save_ps==0){ cerr<<"Error during printing ps in dissolution style."<<endl; return ;}

	cerr<<"Printing network..."<<endl;

	S.check_diss_pattern(S.print_diss_factor);

	int N=(S.N_x);
	int M=(S.N_y);
	double skala = 400./(S.L_out*max(N,M));     //
	double x_zero = 100./skala;
	double y_zero = 750./skala;


	if(S.pages_saved==0) stream << "%!PS-Adobe-3.0" << endl<<"%%Pages:"<<S.pages_tot<<endl<<endl;

	stream <<"%%Page: "<<S.pages_saved+1<<" "<<S.pages_saved+1<<endl<<endl;
	stream <<"1 setlinejoin 1 setlinecap 0.02 setlinewidth"<<endl;

	stream<<skala<<"\t"<<skala<<" scale"<<endl;
	stream<<x_zero<<"\t"<<y_zero<<" translate"<<endl;
	stream<<"1 setlinejoin"<<endl;
	stream<<"1 setlinecap" <<endl;

	//title
	stream<<Kolor(0,0,0);

	stream<<"/Times-Bold findfont "<< 20./skala<<" scalefont setfont"<<endl;
	stream<<"0 "<<-450./skala<<" moveto"<<endl;
	stream<<"0 0 ("<<S.description_note<<") ashow stroke"<<endl<<endl;


	for(int i=0;i<S.NG;i++)   print_grain_with_scaling(stream,*(S.g[i]),S);
	stream << "showpage "<<endl<<flush;

}

