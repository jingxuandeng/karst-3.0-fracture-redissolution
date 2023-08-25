#include "printing.h"
#include "node.h"

bool if_print_labels = true;

ostream & operator<< (ostream &os, Kolor k) {
	os<<k.r<<" "<<k.g<<" "<<k.b<<" setrgbcolor"<<endl;
	return os;}

ostream & operator<< (ostream &os, Linia l){

	os<<l.w<<" setlinewidth"<<endl;
	os<<l.a<<"moveto "<<l.b<<"lineto stroke"<<endl;

	return os;}


ostream & operator<< (ostream &os, Porek l){
	
	os<<setprecision(4)<<l.k;
	os<<l.w<<" setlinewidth"<<endl;
	os<<l.a<<" moveto "<<l.b<<"lineto stroke"<<endl;
	if(if_print_labels) os<<"/Times-Bold findfont "<<(l.b-l.a)/7.<<" scalefont setfont "<<Kolor(0,0,1)<<endl;
	os<<Point((l.b.x+l.a.x)/2, (l.b.y+l.a.y)/2)<<"moveto"<<endl;
	if(if_print_labels){
		if (l.podpis != 666) os <<"0 0 ("<<setprecision(3)<< l.podpis<<") ashow stroke"<<endl;
		else                 os << "stroke"<<endl;}
	else os << "stroke"<<endl;

	return os;}


ostream & operator<< (ostream &os, Point p) {
	if(p.z==0) os<<setprecision(4)<<setw(12)<<p.x<<setw(12)<<-p.y<<" ";
	else       os<<setw(12)<<p.x<<setw(12)<<p.y<<setw(12)<<p.z<<" ";
	return os;}

ofstream_txt & operator<< (ofstream_txt &os, Point p) {
	os<<"("<<setw(12)<<setprecision(4)<<p.x<<setw(12)<<p.y<<setw(12)<<p.z<<")";
	return os;}



//double operator - (Point &p1, Point &p2) {return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));} //return distance between two points
double operator - (Point  p1, Point  p2) {return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));} //return distance between two points
Point  operator * (Point &p1, Point &p2) {return Point((p1.x+p2.x)/2,(p1.y+p2.y)/2,(p1.z+p2.z)/2);}         //return middle of two points
Point  operator * (Point p1,  Point p2)  {return Point((p1.x+p2.x)/2,(p1.y+p2.y)/2,(p1.z+p2.z)/2);}         //return middle of two points
//Point  operator + (Point &p1, Point &p2) {return Point((p1.x+p2.x),(p1.y+p2.y),(p1.z+p2.z));}               //return sum of two points
Point  operator + (Point  p1, Point  p2) {return Point((p1.x+p2.x),(p1.y+p2.y),(p1.z+p2.z));}               //return sum of two points
Point  operator * (double a,  Point &p ) {return Point(a*p.x,a*p.y,a*p.z);}

ostream & operator<< (ostream &os, Kropka kr) {

	if(kr.podpis>=0) 	os<<kr.k;
	else				os<<Kolor(1-kr.k.r,1-kr.k.g,1-kr.k.b);
	os<<kr.a<<kr.r<<" 0 360 arc fill closepath"<<endl;
	os<<Kolor(0,0,1);
	os<<"/Times-Bold findfont "<<kr.r<<" scalefont setfont"<<endl;
	os<<Point(kr.a.x-kr.r/3,kr.a.y-kr.r/3)<<"moveto"<<endl;
	if(if_print_labels){
		if (kr.podpis != 666) os<<"0 0 ("<<setprecision(3)<<kr.podpis<<") ashow stroke"<<endl;}
	else os << "stroke"<<endl;

	return os;}

ostream & operator<< (ostream &os, Kropa kr) {
	
	if     (kr.nod->t==0) 	kr.k=Kolor(0.7,0.7,0.7);
	else if(kr.nod->t==-1) 	kr.k=Kolor(0,0,0);
	else if(kr.nod->t== 1) 	kr.k=Kolor(0,0,1);
	if     (kr.nod->t==0 && kr.nod->cb==1) kr.k=Kolor(0,1,0); //for stream-tube mixing

	if(kr.nod->x>=1)   kr.k=Kolor(0.8,0.2,0.2);
	if(kr.nod->tmp>=0) kr.k=Kolor(0.2,0.2,0.8);


	os<<kr.k;

	os<<kr.a<<kr.r<<" 0 360 arc fill closepath"<<endl;
	os<<Kolor(1,1,1);
	os<<"/Times-Bold findfont "<<kr.r<<" scalefont setfont"<<endl;
	os<<Point(kr.a.x-kr.r/3,kr.a.y-kr.r/3)<<"moveto"<<endl;
	if(if_print_labels){
		if (kr.nod->tmp != 666) os<<"0 0 ("<<setprecision(5)<<kr.nod->tmp<<") ashow stroke"<<endl;}
	else os << "stroke"<<endl;
	
	return os;}


ostream & operator<< (ostream &os, Trojkacik tr) {
	
	os<<tr.k;
	os<<tr.n1<<"moveto "<<tr.n2<<"lineto "<<tr.n3<<"lineto closepath fill stroke"<<endl;
	os<<Kolor(0,0,0);//os<<Kolor(0.5,0.5,0.5);
	//os<<"/Times-Bold findfont "<<(tr.n1-tr.n2)/5.<<" scalefont setfont"<<endl;
	os<<"/Times-Bold findfont "<<0.1<<" scalefont setfont"<<endl;
	os<<Point((tr.n1.x+tr.n2.x+tr.n3.x)/3,(tr.n1.y+tr.n2.y+tr.n3.y)/3)<<"moveto"<<endl;
	if(if_print_labels){
		if (tr.podpis != 666) os<<"0 0 ("<< setprecision(3) <<tr.podpis<<") ashow stroke"<<endl;}
	else os << "stroke"<<endl;
	
	return os;}

ostream & operator<< (ostream &os, Wielobok w) {

	os<<w.k;
	os<<w.p[0]<<"moveto ";
	for (int i=1;i<w.b;i++)
	os<<w.p[i]<<"lineto ";
	os<<"closepath fill stroke"<<endl;
	os<<Kolor(0,0,0);//os<<Kolor(0.5,0.5,0.5);
	//os<<"/Times-Bold findfont "<<(tr.n1-tr.n2)/5.<<" scalefont setfont"<<endl;
	//os<<"/Times-Bold findfont "<<0.1<<" scalefont setfont"<<endl; FIXME:edit 28.11.2022
	Point sr = Point();
	for (int i=0;i<w.b;i++) sr = sr +w.p[i];
	os<<Point(sr.x/w.b,sr.y/w.b)<<"moveto"<<endl;
	if(if_print_labels){
		if (w.podpis != 666) os<<"0 0 ("<< setprecision(3) <<w.podpis<<") ashow stroke"<<endl;}
	else os << "stroke"<<endl;

	return os;}


