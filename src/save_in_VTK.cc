/*
 * save_in_VTK.cc
 *
 *  Created on: Mar 23, 2020
 *      Author: Rishabh Sharma
 */

#include "network.h"



void Network::write_vtk_data()
{
//	cout << "Opening VTK file to write" << endl;

	ofstream network_vtk;
	network_vtk.open("network_"+to_string(tot_steps)+".vtk", ios_base::out | ios_base::trunc);

	network_vtk<< "# vtk DataFile Version 1.0" << endl;
	network_vtk<< "2D Network model"<< endl <<"ASCII"<< endl << endl <<"DATASET UNSTRUCTURED_GRID"<<endl <<"POINTS"<<" "<<NN<<" "<<"float"<<endl;


	cout <<"Writing VTK data to files" << endl;
	write_point_data("network_"+to_string(tot_steps)+".vtk");
	write_cell_data("network_"+to_string(tot_steps)+".vtk");
	write_diameter("network_"+to_string(tot_steps)+".vtk");
	write_flow_rate("network_"+to_string(tot_steps)+".vtk");
	write_concentration("network_"+to_string(tot_steps)+".vtk");
	write_concentration_c("network_"+to_string(tot_steps)+".vtk");
	network_vtk.close();

	write_grains_vtk_data();
}

/**
* This function writes a separate VTK file (grains_<step>.vtk) representing the material
* grains as polygonal cells, with the grain volume fields attached as CELL_DATA:
*   Va   - volume of mineral A
*   Va1  - volume of the less reactive mineral A1
*   Vx   - volume of the non-reacting material
*   Ve   - volume of the precipitated mineral E
*   Vtot - Va + Va1 + Vx + Ve (total solid volume left in the grain)
* Grains with fewer than 3 nodes, or grains that wrap around a periodic boundary, are
* written as a degenerate (invisible) cell so the data arrays stay aligned with NG.
*
* @author Jingxuan Deng
* @date 02/09/2026
*/
void Network::write_grains_vtk_data()
{
	string fname = "grains_"+to_string(tot_steps)+".vtk";

	ofstream f;
	f.open(fname, ios_base::out | ios_base::trunc);
	if(f.is_open() == false){ cout << "Problem in writing grain VTK file" << endl; return; }

	cout << "Writing grain volume VTK data to " << fname << endl;

	for (int i = 0; i < NN; ++i) n[i]->tmp = i;   //node -> point index

	double max_dist = (N_x > 0) ? N_x/2.0 : 1e30; //consecutive grain nodes further apart => periodic wrap

	//classify each grain: real polygon (bN>=3, no wrap) or degenerate placeholder
	bool *drawn = new bool[NG];
	long  conn_size = 0;
	for (int i = 0; i < NG; ++i){
		Grain *gg = g[i];
		bool ok = (gg->bN >= 3);
		for (int b = 0; ok && b < gg->bN; ++b)
			if ( (gg->n[b]->xy - gg->n[(b+1)%gg->bN]->xy) > max_dist ) ok = false;
		drawn[i]   = ok;
		conn_size += (ok ? gg->bN : 3) + 1;
	}

	f << "# vtk DataFile Version 1.0" << endl;
	f << "2D Network grains model" << endl << "ASCII" << endl << endl << "DATASET UNSTRUCTURED_GRID" << endl;

	f << "POINTS " << NN << " float" << endl;
	for (int j = 0; j < NN; ++j)
		f << n[j]->xy.x << " " << n[j]->xy.y << " " << n[j]->xy.z << endl;

	f << endl << "CELLS " << NG << " " << conn_size << endl;
	for (int i = 0; i < NG; ++i){
		Grain *gg = g[i];
		if (drawn[i]){
			f << gg->bN;
			for (int b = 0; b < gg->bN; ++b) f << " " << (int)gg->n[b]->tmp;
			f << endl;
		}
		else{
			int idx = (gg->bN > 0) ? (int)gg->n[0]->tmp : 0;
			f << "3 " << idx << " " << idx << " " << idx << endl;
		}
	}

	f << endl << "CELL_TYPES " << NG << endl;
	for (int i = 0; i < NG; ++i) f << 7 << endl;   //7 == VTK_POLYGON

	f << endl << "CELL_DATA " << NG << endl;

	f << "SCALARS Va float 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < NG; ++i) f << g[i]->Va << endl;

	f << "SCALARS Va1 float 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < NG; ++i) f << g[i]->Va1 << endl;

	f << "SCALARS Vx float 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < NG; ++i) f << g[i]->Vx << endl;

	f << "SCALARS Ve float 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < NG; ++i) f << g[i]->Ve << endl;

	f << "SCALARS Vtot float 1" << endl << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < NG; ++i) f << (g[i]->Va + g[i]->Va1 + g[i]->Vx + g[i]->Ve) << endl;

	delete[] drawn;
	f.close();
}

void Network::write_point_data(string file_name)
{

	ofstream tmp;
	ofstream *obj = NULL;
	tmp.open(file_name, ios::app);
	if(tmp.is_open() == false) cout << "Problem in writing VTK files" << endl;
	else obj = &tmp;

	ofstream &file = *obj;

	for (int j = 0; j < NN; ++j) file<< n[j]->xy.x << setw(10) << n[j]->xy.y << setw(10) << n[j]->xy.z << endl;
	file.close();
}

void Network::write_line_data(string file_name)
{
	ofstream tmp;
	ofstream *obj = NULL;
	tmp.open(file_name, ios::app);
	if(tmp.is_open() == false) cout << "Problem in writing VTK files" << endl;
	else obj = &tmp;
	ofstream &file = *obj;

	for (int i =0; i<NN; i++) n[i]->tmp = i;

	file << endl << endl << "LINES" << "  " <<NP << "  " << 3*NP <<endl;

	for (int j = 0; j < NP; ++j){ //Requires further modification for limit of length of pores
			if((p[j]->d != 0) || (fabs(p[j]->n[0]->xy - p[j]->n[1]->xy) < N_x/2.0))
				file<< "2" <<" "<< p[j]->n[0]->tmp << " " << p[j]->n[1]->tmp << endl;
			else
				file<< "2" <<" "<< "0" << " " << "0" << endl;}
	file.close();
}

void Network::write_cell_data(string file_name)
{
	int VTK_LINES = 3;
	ofstream tmp;
	ofstream *obj = NULL;
	tmp.open(file_name, ios::app);
	if(tmp.is_open() == false) cout << "Problem in writing VTK files" << endl;
	else obj = &tmp;
	ofstream &file = *obj;

	for (int i =0; i<NN; i++) n[i]->tmp = i;


//Writing Cell connectivity
	file << endl << endl << "CELLS" << "  " <<NP << "  " << 3*NP <<endl;

	for (int j = 0; j < NP; ++j){ //Requires further modification for limit of length of pores
			if((p[j]->d != 0) && (fabs(p[j]->n[0]->xy - p[j]->n[1]->xy) < N_x/2.0))
											{ file<< "2" <<" "<< p[j]->n[0]->tmp << " " << p[j]->n[1]->tmp << endl; }
			else  file<< "2" <<" " << "0" << " " << "0" << endl;
		}

//Writing Cell types
	file << endl << "CELL_TYPES" << " " << NP << endl;

	for (int i = 0; i < NP; ++i) file << VTK_LINES << endl;			//CELL_TYPES VTK_LINES = 3;

	file.close();
}


void Network::write_diameter(string file_name)
{

	ofstream os_tmp;
	ofstream *os_p = NULL;

	os_tmp.open(file_name, ios::app);
	if (os_tmp.is_open() == false) cout << "Problem in writing VTK files" << endl;
	else os_p = &os_tmp;

	ofstream &file = *os_p;
	file << endl << "CELL_DATA" << " " << NP <<endl <<"SCALARS Diameter float" << endl << "LOOKUP_TABLE custom_table" << endl;

	for (int i = 0; i < NP; ++i)	file<< p[i]->d << endl;

	file.close();
}


void Network::write_flow_rate(string file_name)
{

	ofstream os_tmp;
	ofstream *os_p = NULL;

	os_tmp.open(file_name, ios::app);
	if (os_tmp.is_open() == false) cout << "Problem in writing VTK files" << endl;
	else os_p = &os_tmp;

	ofstream &file = *os_p;
	file << endl <<"SCALARS Flow_Rate float" << endl << "LOOKUP_TABLE custom_table" << endl;

	for (int i = 0; i < NP; ++i)	file<< p[i]->q << endl;

	file.close();
}


void Network::write_concentration(string file_name)
{

	ofstream os_tmp;
	ofstream *os_p = NULL;

	os_tmp.open(file_name, ios::app);
	if (os_tmp.is_open() == false) cout << "Problem in writing VTK files" << endl;
	else os_p = &os_tmp;

	ofstream &file = *os_p;
	file << endl <<"SCALARS Concentration float" << endl << "LOOKUP_TABLE custom_table" << endl;

	if(if_streamtube_mixing) for (int i = 0; i < NP; ++i)	file<< p[i]->c_in << endl;
	else for (int i = 0; i < NP; ++i)	file<< p[i]->n[0]->cb << endl;
	file.close();
}

void Network::write_concentration_c(string file_name)
{

	ofstream os_tmp;
	ofstream *os_p = NULL;

	os_tmp.open(file_name, ios::app);
	if (os_tmp.is_open() == false) cout << "Problem in writing VTK files" << endl;
	else os_p = &os_tmp;

	ofstream &file = *os_p;
	file << endl <<"SCALARS Concentration_C float" << endl << "LOOKUP_TABLE custom_table" << endl;

	if(if_streamtube_mixing) for (int i = 0; i < NP; ++i)	file<< Cc_0 << endl; // this is not working in streamline routing. Need modification
	else for (int i = 0; i < NP; ++i)	file<< p[i]->n[0]->cc << endl;
	file.close();
}



