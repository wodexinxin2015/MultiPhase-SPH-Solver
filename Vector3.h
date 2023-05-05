/************************************************************************************
Multi-Phase and Parallelized SPH program
--parallelized and optimized for SMP (Symmetric multiprocessing ) architecture.
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#ifndef __VECTOR3_H_INCLUDED__
#define __VECTOR3_H_INCLUDED__

#include <math.h>
class Vector3
{
public:
	double x, y, z;

	Vector3()
	{
		x = y = z = 0.0;
	}

	Vector3(const Vector3 &a) : x(a.x), y(a.y), z(a.z) {}
	Vector3(double nx, double ny, double nz) : x(nx), y(ny), z(nz) {}

	//zero vector
	void Zero()
	{
		x = y = z = 0.0;
	}

	//normalization of vectot

	void normalize()
	{
		double magsq = x * x + y * y + z * z; //magnitue
		if (magsq > 0.0)					  //
		{
			double temp = 1.0 / sqrt(magsq);
			x *= temp;
			y *= temp;
			z *= temp;
		}
	}

	//vector -
	Vector3 operator-() const { return Vector3(-x, -y, -z); }

	Vector3 operator*(double a /*vector * a*/) const
	{
		return Vector3(x * a, y * a, z * a);
	}
	Vector3 operator/(double a) const
	{
		double temp = 1.0 / a;
		return Vector3(x * temp, y * temp, z * temp);
	}

	Vector3 operator*=(double a)
	{
		x *= a;
		y *= a;
		z *= a;
		return *this;
	}

	Vector3 operator/=(double a)
	{
		double temp = 1.0 / a;
		x *= temp;
		y *= temp;
		z *= temp;
		return *this;
	}

	Vector3 operator+(const Vector3 &a) const
	{
		return Vector3(x + a.x, y + a.y, z + a.z);
	}

	Vector3 operator+=(const Vector3 &a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	}

	Vector3 operator-=(const Vector3 &a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;
		return *this;
	}

	Vector3 operator-(const Vector3 &a) const
	{
		return Vector3(x - a.x, y - a.y, z - a.z);
	}

	//dot production

	double operator*(const Vector3 &a) const
	{
		return x * a.x + y * a.y + z * a.z;
	}
};

//magnitude of vector
inline double vectorMag(const Vector3 &a)
{
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

//vector * k
inline Vector3 operator*(double k, const Vector3 &v)
{
	return Vector3(k * v.x, k * v.y, k * v.z);
}

//distance between two vector
inline double distance(const Vector3 &a, const Vector3 &b)
{
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	double dz = a.z - b.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}

//cross production
inline Vector3 crossProduct(const Vector3 &a, const Vector3 &b)
{
	return Vector3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}
#endif