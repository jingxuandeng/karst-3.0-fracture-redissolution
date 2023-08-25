#ifndef H_TRIANGULATION
#define H_TRIANGULATION

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <limits>

//numeric.h

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp=2);

template<class T>
T half(T x){}

template <>
float half(float x);

template <>
double half(double x);


//vector2.h
template <typename T>
class Vector2
{
	public:
		//
		// Constructors
		//

		Vector2():x(0), y(0), a(0){}

		Vector2(T _x, T _y): x(_x), y(_y), a(0){}
		Vector2(T _x, T _y, int _a): x(_x), y(_y), a(_a){}

		Vector2(const Vector2 &v): x(v.x), y(v.y), a(v.a){}

		//
		// Operations
		//
		T dist2(const Vector2 &v) const
		{
			T dx = x - v.x;
			T dy = y - v.y;
			return dx * dx + dy * dy;
		}

		T dist(const Vector2 &v) const
		{
			return sqrt(dist2(v));
		}

		T norm2() const
		{
			return x * x + y * y;
		}

		T x;
		T y;
		int a;  //point name

};

template <>
float Vector2<float>::dist(const Vector2<float> &v) const ;

template <>
double Vector2<double>::dist(const Vector2<double> &v) const ;

template<typename T>
std::ostream &operator << (std::ostream &str, Vector2<T> const &point);

template<typename T>
bool operator == (const Vector2<T>& v1, const Vector2<T>& v2);

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(const Vector2<T>& v1, const Vector2<T>& v2, int ulp=2);


//edge.h


template <class T>
class Edge
{
	public:
		using VertexType = Vector2<T>;

		Edge(const VertexType &p1, const VertexType &p2) : p1(p1), p2(p2), isBad(false) {};
		Edge(const Edge &e) : p1(e.p1), p2(e.p2), isBad(false) {};

		VertexType p1;
		VertexType p2;

		bool isBad;
};

template <class T>
inline std::ostream &operator << (std::ostream &str, Edge<T> const &e);

template <class T>
inline bool operator == (const Edge<T> & e1, const Edge<T> & e2);
template <class T>
inline bool almost_equal (const Edge<T> & e1, const Edge<T> & e2);


//triangle.h

template <class T>
class Triangle
{
	public:
		using EdgeType = Edge<T>;
		using VertexType = Vector2<T>;

		Triangle(const VertexType &_p1, const VertexType &_p2, const VertexType &_p3)
		:	p1(_p1), p2(_p2), p3(_p3),
			e1(_p1, _p2), e2(_p2, _p3), e3(_p3, _p1), isBad(false)
		{}

		bool containsVertex(const VertexType &v) const
		{
			// return p1 == v || p2 == v || p3 == v;
			return almost_equal(p1, v) || almost_equal(p2, v) || almost_equal(p3, v);
		}

		bool circumCircleContains(const VertexType &v) const
		{
			const T ab = p1.norm2();
			const T cd = p2.norm2();
			const T ef = p3.norm2();

			const T circum_x = (ab * (p3.y - p2.y) + cd * (p1.y - p3.y) + ef * (p2.y - p1.y)) / (p1.x * (p3.y - p2.y) + p2.x * (p1.y - p3.y) + p3.x * (p2.y - p1.y));
			const T circum_y = (ab * (p3.x - p2.x) + cd * (p1.x - p3.x) + ef * (p2.x - p1.x)) / (p1.y * (p3.x - p2.x) + p2.y * (p1.x - p3.x) + p3.y * (p2.x - p1.x));

			const VertexType circum(half(circum_x), half(circum_y));
			const T circum_radius = p1.dist2(circum);
			const T dist = v.dist2(circum);
			return dist <= circum_radius;
		}

		VertexType p1;
		VertexType p2;
		VertexType p3;
		EdgeType e1;
		EdgeType e2;
		EdgeType e3;
		bool isBad;
};

template <class T>
inline std::ostream &operator << (std::ostream &str, const Triangle<T> & t);

template<typename T>
inline bool operator == (const Vector2<T>& v1, const Vector2<T>& v2)
{
	//return (v1.x == v2.x) && (v1.y == v2.y);
	return (v1.a==v2.a);
}

template <class T>
inline bool operator == (const Triangle<T> &t1, const Triangle<T> &t2)
{
	return	(t1.p1 == t2.p1 || t1.p1 == t2.p2 || t1.p1 == t2.p3) &&
			(t1.p2 == t2.p1 || t1.p2 == t2.p2 || t1.p2 == t2.p3) &&
			(t1.p3 == t2.p1 || t1.p3 == t2.p2 || t1.p3 == t2.p3);
}

template <class T>
inline bool operator == (const Edge<T> & e1, const Edge<T> & e2)
{
	return 	(e1.p1 == e2.p1 && e1.p2 == e2.p2) ||
			(e1.p1 == e2.p2 && e1.p2 == e2.p1);
}


template <class T>
inline bool almost_equal(const Triangle<T> &t1, const Triangle<T> &t2);

void triangulation(int N_x, int N_y, \
		std::vector<Vector2 <double> >  &points_out, \
		std::vector<Vector2 <double> >  &points_tmp_out, \
	    std::vector<Edge    <double> >  &edges_out, \
		std::vector<Triangle<double> >  &triangles_out, \
		bool if_regular_points, bool if_periodic_bc, double random_seed=-1);

#endif
