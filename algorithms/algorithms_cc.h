#ifndef ALGORITHMS_cc_H
#define ALGORITHMS_cc_H 


#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>

template <class T> class Vector2;
template <class T> class Edge;
template <class T> class Triangle;

extern "C" {
	int solve_matrix(int W, int R_N, int* w_k, int* w_i, double* B, double* y);
}


void triangulation(int N_x, int N_y, \
		std::vector<Vector2< float> >  &points_out, \
		std::vector<Vector2 <float> >  &points_tmp_out, \
	    std::vector<Edge    <float> >  &edges_out, \
		std::vector<Triangle<float> >  &triangles_out, \
		bool if_regular_points, bool if_periodic_bc, double random_seed );

#endif


