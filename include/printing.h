#ifndef DRUKOWANIE_H
#define DRUKOWANIE_H 


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include "constants.h"

// #include "network.h"
class Node;
class Pore;
class Grain;
class Network;


using namespace std;


class Kolor {
	public:
		double r,g,b;
	public:
		Kolor(double rr=0,double gg=0, double bb=0){r=rr; g=gg; b=bb;}	
		
		friend ostream & operator<< (ostream &os, Kolor k);
};


class Point{
	
	public:	
		double x, y, z;
	
	public:
		Point( double xx = 0., double yy = 0., double zz = 0. ){
    		x = xx; 
    		y = yy;
    		z = zz;
  		}
		
		friend ostream      & operator<< (ostream      &os, Point p);
		friend ofstream_txt & operator<< (ofstream_txt &os, Point p);
};


class Linia{
	
	protected:	
		Point a, b;
		double w;
	
	public:
		Linia( Point aa = Point(), Point bb = Point(), double ww = 0.1){
    		a = aa; 
    		b = bb;
			w = ww;
  		}
		
		friend ostream & operator<< (ostream &os, Linia l);
};


class Porek : public Linia{
	
	private:	
		float podpis;
		Kolor k;
	
	public:
		Porek( Point aa = Point(), Point bb = Point(), double ww = 0.1, float pp = 0, Kolor kk=Kolor() ){
    		a = aa; 
    		b = bb;
			w = ww;
			podpis=pp; k=kk;
  		}
		
		friend ostream & operator<< (ostream &os, Porek l);
};

class Kropka
{
	protected:
		Point a;
		Kolor k;
		float podpis;
		double r;
	public:
		Kropka(Point aa = Point(), float i=666, Kolor kk = Kolor(), double rr=1)
		{
			a = aa;
			k = kk;
			podpis=i;
			r=rr;
		}
		
		friend ostream & operator<< (ostream &os, Kropka kr);
};
	
class Kropa : public Kropka
{
	public:
		Node *nod;
	
	public:
		Kropa(Point aa = Point(), Node * nnod=NULL, Kolor kk = Kolor(), double rr=1)
		{
			a = aa;
			k = kk;
			nod=nnod;
			r=rr;
			podpis = -8;
		}
		
		friend ostream & operator<< (ostream &os, Kropa kr);
};


class Trojkacik{
	
	public:	
		Point n1,n2,n3;
	
		float podpis;
		Kolor k;
		
	
	public:
		Trojkacik( Point nn1 = Point(), Point nn2 = Point(), Point nn3 = Point(), float ppodpis = 0, Kolor kk=Kolor() ){
    		n1 = nn1; 
    		n2 = nn2; 
    		n3 = nn3;

			podpis=ppodpis; k=kk;
  		}
		
		friend ostream & operator<< (ostream &os, Trojkacik tr);
};


class Wielobok{

	public:
		int    b;     ///<  nr of vertices
		Point* p;     ///<  list of vertices

		float podpis;
		Kolor k;


	public:
		Wielobok (int bb, Point *pp, float ppodpis = 0, Kolor kk=Kolor() ){

			b = bb;
			p = pp;
			podpis = ppodpis; k=kk;
  		}

		friend ostream & operator<< (ostream &os, Wielobok w);
};



ofstream_ps & operator << (ofstream_ps & stream, Node    &n);
ofstream_ps & operator << (ofstream_ps & stream, Pore    &p);
ofstream_ps & operator << (ofstream_ps & stream, Grain   &g);
ofstream_ps & operator << (ofstream_ps & stream, Network &S);
ofstream_ps & print_grain_with_scaling (ofstream_ps & stream, Grain &g, Network &S);


///< double operator - (Point  &p1, Point &p2); ///< return distance between two points
double operator - (Point   p1, Point  p2); ///< return distance between two points
///< Point  operator + (Point  &p1, Point &p2); ///< return sum of two points
Point  operator + (Point   p1, Point  p2); ///< return sum of two points
Point  operator * (Point  &p1, Point &p2); ///< return middle of two points
Point  operator * (Point   p1, Point  p2); ///< return middle of two points
Point  operator * (double a,  Point &p );  ///< multiply point by double


void Print_network_in_debugging_style     (ofstream_ps &stream, Network &S);
void Print_network_in_dissolution_style   (ofstream_ps &stream, Network &S);
void Print_network_in_debugging_style_tmp (ofstream_ps &stream, Network &S,int nr_tmp);
void Print_network_in_grain_style         (ofstream_ps &stream, Network &S);

#endif


