#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
#include <iomanip>
#include "Math270A_DoubleArray2D.h" 

using namespace std; // Using the "standard" (std) standard template library components
#ifndef _GridFun2D_
#define _GridFun2D_

#ifdef _OPENMP
#include <omp.h>
#endif

class GridFun2D
{
public:

	GridFun2D()
	{
		initialize();
	};

	GridFun2D(const GridFun2D& G)
	{
		initialize(G);
	};

	GridFun2D(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax)
	{
		initialize(xPanel, xMin, xMax, yPanel, yMin, yMax);
	};

	void initialize()
	{
		xMin = 0;
		xMax = 0;
		xPanel = 1;
		hx = 0;
		yMin = 0;
		yMax = 0;
		yPanel = 1;
		hy = 0;
	};

	void initialize(const GridFun2D& G)
	{
		hx = G.hx;
		xMin = G.xMin;
		xMax = G.xMax;
		xPanel = G.xPanel;
		hy = G.hy;
		yMin = G.yMin;
		yMax = G.yMax;
		yPanel = G.yPanel;
		values = G.values;

	};

	void initialize(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax)
	{
		this->xPanel = xPanel;
		this->xMin = xMin;
		this->xMax = xMax;
		this->yPanel = yPanel;
		this->yMin = yMin;
		this->yMax = yMax;
		this->hy = (yMax - yMin) / (yPanel);
		this->hx = (xMax - xMin) / (xPanel);
		values.initialize(xPanel + 1, yPanel + 1);

	};

	void operator=(const GridFun2D& G)
	{
		initialize(G);
	};
#ifndef _OPENMP
	void operator+=(const GridFun2D& G)
	{
		double* dValues = values.getDataPointer();
		double* gValues = G.values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] += gValues[i];
		}
	};
#else

	void operator+=(const GridFun2D& G)
	{
		double* dValues = values.getDataPointer();
		double* gValues = G.values.getDataPointer();
		long i;
#pragma omp parallel for default(shared) private(i) schedule(static)
		for (i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] += gValues[i];
		}
	};

#endif


	void operator-=(const GridFun2D& G)
	{
		double* dValues = values.getDataPointer();
		double* gValues = G.values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] -= gValues[i];
		}
	};

	void operator*=(const double alpha)
	{
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] *= alpha;
		}
		return;
	};

	void operator/=(const double alpha)
	{
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] /= alpha;
		}
		return;
	};

	void setToValue(double d)
	{
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			dValues[i] = d;
		}
		return;
	};

	double normInf()
	{
		//fix
		double norm = 0.0;
		double* dValues = values.getDataPointer();
		for (long i = 0; i < values.getDataSize(); i++)
		{
			if (fabs(dValues[i]) > norm)
				norm = fabs(dValues[i]);
		}
		return norm;
	};

	friend ostream& operator<<(ostream& outStream, const GridFun2D& V)
	{
		//fix
		double* dValues = V.values.getDataPointer();
		for (long i = 0; i < V.values.getDataSize(); i += 10)
		{
			outStream << setw(5) << dValues[i] << " ";
			outStream << endl;
		}
		return outStream;
	};

	void extractXslice(long yIndex, vector<double>& u1Dx) const
	{
		u1Dx.resize(xPanel + 1);
		for (long i = 0; i <= xPanel; i++)
		{
			u1Dx[i] = this->values(i, yIndex);
		}
	};


	void extractYslice(long xIndex, vector<double>& u1Dy) const
	{
		u1Dy.resize(yPanel + 1);
		for (long j = 0; j <= yPanel; j++)
		{
			u1Dy[j] = this->values(xIndex, j);
		}
	};

	void insertXslice(long yIndex, const vector<double>& u1Dx)
	{
		for (long i = 0; i <= xPanel; i++)
		{
			this->values(i, yIndex) = u1Dx[i];
		}
	};

	void insertYslice(long xIndex, const vector<double>& u1Dy)
	{
		for (long j = 0; j <= yPanel; j++)
		{
			this->values(xIndex, j) = u1Dy[j];
		}
	};



	double hx;
	double xMin;
	double xMax;
	long xPanel;
	double hy;
	double yMin;
	double yMax;
	long yPanel;
	Math270A::DoubleArray2D values;
};

#endif