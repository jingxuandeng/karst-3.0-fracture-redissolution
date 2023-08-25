#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>



//#define VISUALIZATION //only for debugging

#ifdef VISUALIZATION
#include <SFML/Graphics.hpp>
#endif


#include "delaunay-triangulation-master/vector2.h"
#include "delaunay-triangulation-master/triangle.h"
#include "delaunay-triangulation-master/delaunay.h"


using namespace std;


#ifdef VISUALIZATION
int make_visualization(const std::vector<Triangle<double> >& triangles, const std::vector<Edge<double> >& edges, const vector<Vector2<double> >& points);
#endif

float RandomFloat(float a, float b) {
    const float random = ((float) rand()) / (float) RAND_MAX;
    const float diff = b - a;
    const float r = random * diff;
    return a + r;
}

bool f_weights(const Vector2<double> &v1, const Vector2<double> &v2){

	if       ((int) v1.y < (int) v2.y) return true;
	else if  ((int) v1.y > (int) v2.y) return false;
	else                               return v1.x <  v2.x;
}

void triangulation(int N_x, int N_y, \
		std::vector<Vector2 <double> >  &points_out, \
		std::vector<Vector2 <double> >  &points_tmp_out, \
	    std::vector<Edge    <double> >  &edges_out, \
		std::vector<Triangle<double> >  &triangles_out, \
		bool if_regular_points, bool if_periodic_bc, double random_seed){


	double eps =1.;  //parameter for regular network threshold
	bool if_visualization = true;
	bool if_verbose_mode  = false;


	int NN = N_x*N_y; //total number of nodes

	std::vector<Vector2<double> > points;
	//select nodes

	if (random_seed==-1) srand (time(NULL));
	else                 srand (random_seed);
	cerr<<"Random seed: "<<random_seed<<endl;

	if (if_regular_points){

		cerr << "Generating " << NN<< " regular points." << endl;
		for(int j=0; j<N_y;j++) for(int i=0; i<N_x; i++)  {
			points.push_back(Vector2<double>(i + RandomFloat(-eps, eps)+0.5*(j%2), j + RandomFloat(-eps, eps),i+j*N_x));
		}
	}
	else{
		cerr << "Generating " << NN<< " random points." << endl;
		for(int i=0; i<NN;i++) {
			points.push_back(Vector2<double>(RandomFloat(0, N_x), RandomFloat(0, N_y),i));
		}
		cerr<<"Sorting points (for more readable printing)..."<<endl;
		std::stable_sort (points.begin(), points.end(), f_weights);
		int i=0;
		for(auto &p : points) p.a = i++;
	}
	if (if_verbose_mode){
		cerr<<endl<<"After generation points."<<endl;
		for(const auto &p : points)  cerr << p << std::endl;
	}

	std::vector<Vector2<double> > points_tmp;
	points_tmp = points;
	int eps_bc = 20;   //vicinity of mirror nodes to get periodic boundary conditions
	if(if_periodic_bc){
		//doubling the point set to get boundary conditions
		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x+N_x,p.y    ,p.a  +NN));
		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x-N_x,p.y    ,p.a-2*NN));
		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x,    p.y+N_y,p.a+3*NN));
		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x,    p.y-N_y,p.a-4*NN));

		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x-N_x,p.y-N_y,p.a-5*NN));
		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x+N_x,p.y+N_y,p.a+6*NN));
		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x-N_x,p.y+N_y,p.a+7*NN));
		for(const auto &p : points_tmp)		points.push_back(Vector2<double>(p.x+N_x,p.y-N_y,p.a-8*NN));

	    for (auto it = points.begin(); it != points.end(); ) {
	        Vector2<double> p = *it;
	        if(p.x>N_x+eps_bc || p.x<-eps_bc || p.y>N_y+eps_bc || p.y<-eps_bc) it = points.erase(it);
	        else ++it;
	    }

	}
	if(if_verbose_mode){
		cerr<<endl<<"After adding periodic bc:"<<endl;
		for(const auto &p : points)  cerr << p << std::endl;
	}

	cerr<<"Triangulation..."<<endl;
	Delaunay<double> triangulation;
	std::vector<Triangle<double> > triangles = triangulation.triangulate(points);
	std::cout << triangles.size() << " triangles generated\n";
	std::vector<Edge<double> > edges = triangulation.getEdges();


	if(if_verbose_mode){
		std::cout << " ========= ";

		std::cout << "\nPoints : " << points.size() << std::endl;
		for(const auto &p : points)
			std::cout << p << std::endl;

		std::cout << "\nTriangles : " << triangles.size() << std::endl;
		for(const auto &t : triangles)
			std::cout << t << std::endl;

		std::cout << "\nEdges : " << edges.size() << std::endl;
		for(const auto &e : edges)
			std::cout << e << std::endl;
	}

	#ifdef VISUALIZATION
	make_visualization(triangles,edges,points_tmp);
	#endif
	edges_out       = edges;
	triangles_out   = triangles;
	points_out      = points;
	points_tmp_out  = points_tmp;


	return;
}




#ifdef VISUALIZATION
int make_visualization(const std::vector<Triangle<double> >& triangles, const std::vector<Edge<double> >& edges, const std::vector<Vector2<double> >& points  ){

	double max_x=0, min_x=0, max_y=0, min_y=0;
	for(const auto &t : triangles) for (const auto &p : {t.p1,t.p2,t.p3}){
		if (p.x>max_x) max_x=p.x;
		if (p.x<min_x) min_x=p.x;
		if (p.y>max_y) max_y=p.y;
		if (p.y<min_y) min_y=p.y;
}
	//min_x-=1; min_y-=1; max_y+=1; max_y+=1;
	int delta_x = max_x - min_x ;
	int delta_y = max_y - min_y ;

	// SFML window
	sf::RenderWindow window(sf::VideoMode((delta_x)*10, (delta_y)*10), "Delaunay triangulation");

	// Transform each points of each vector as a rectangle
	std::vector<sf::RectangleShape*> squares;

	for(const auto p : points) {
		sf::RectangleShape *c1 = new sf::RectangleShape(sf::Vector2f(4, 4));
		c1->setPosition((p.x-min_x)*10, (p.y-min_y)*10);
		squares.push_back(c1);
	}

	// Make the lines
	std::vector<std::array<sf::Vertex, 2> > lines;
	for(const auto &e : edges) {
		lines.push_back({{
			sf::Vertex(sf::Vector2f((e.p1.x-min_x)*10 + 2, (e.p1.y - min_y)*10 + 2)),
			sf::Vertex(sf::Vector2f((e.p2.x-min_x)*10 + 2, (e.p2.y - min_y)*10 + 2))
		}});
	}

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear();

		// Draw the squares
		for(const auto &s : squares) {
			window.draw(*s);
		}

		// Draw the lines
		for(const auto &l : lines) {
			window.draw(l.data(), 2, sf::Lines);
		}

		window.display();
	}

	return 0;
}
#endif

template std::ostream& operator<< <double>(std::ostream&, Edge<double> const&);
