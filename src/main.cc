/**
 * @author  Agnieszka Budek
 * @date 25/09/2019
 *
**/


#include "network.h"
#include "printing.h"
#include "tests.h"



int test_oscillation_in_diss_pre(string config_name){

	if (config_name == "")  config_name = "config.txt";

	Network *S = new Network(config_name);		//creation the system: reading initial parameters form the file

	double gamma = S->gamma;
	int N = 10;

	for(int i=0;i<N;i++){
		cerr<<"\n\n\n\nDissolution nr " <<i<<endl;
		S->Cb_0  = 1;
		S->Cc_0  = 0;
		S->gamma = 0;
		S->evolution(0); 							//evolution of the system

		cerr<<"\n\n\n\nPrecipitation nr " <<i<<endl;
		S->Cb_0  = 0;
		S->Cc_0  = 1;
		S->gamma = gamma;
		S->evolution(0); 							//evolution of the system
	}

	delete S;    								//closing the system

	return 0;

}


/**
 * Basic simulation run with both dissolution and precipitation (depending on the config file)
 **/

int test_dissolution(string config_name){

	if (config_name == "")  config_name = "config.txt";

	Network *S = new Network(config_name);		//creation the system: reading initial parameters form the file

	S->evolution(0); 							//evolution of the system

	delete S;    								//closing the system

	return 0;

}

int main(int argc, char** argv){

	
	bool if_show_picture = 0 ;  //if true, picture is automatically shown (working on MacOs, Linux to be checked)
	string config_name   = "";  //path to the config file
    string sim_version   = "single_run";

//inline options: normal mode, debugging mode, number of simulations, etc.

    for (int i=0;i<argc;i++) if (argv[i][0]=='-') switch (argv[i][1]) {

		case 'P':   if_show_picture = atof(argv[++i]);   fprintf(stderr,"if_show_picture=%d\n",if_show_picture); break;

		case 'C':   config_name     = string(argv[++i]); cerr<<"config_name = "<<config_name<<endl; break;

		case 'V':   sim_version     = string(argv[++i]); cerr<<"sim_version = "<<sim_version<<endl; break;

		default:	printf("Unknown argument.\n"); break;
	}


    //which version of simulateon to run
    if       (sim_version == "single_run")    test_dissolution            (config_name);
    else  if (sim_version == "oscillation")	  test_oscillation_in_diss_pre(config_name);
    else                                      cerr<<"ERROR: Improper value of sim_version!"<<endl;

//    test_triangulation();

//additional actions (working on MacOs)
	if(if_show_picture or true){
		cerr<<"Preparing pictures..."<<endl;
		system("ps2pdf net.ps");
		//system("open net.pdf");
	}
	

	return 0;

}

