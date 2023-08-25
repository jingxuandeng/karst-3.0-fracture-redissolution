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
	network_vtk.close();
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



