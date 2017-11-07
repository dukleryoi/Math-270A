#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
#include <iomanip>

#include "GridFun2D.h"  // 1D grid function class
#include "TriSolver.h"  // 1D Crank-Nicholson relaxation operator class
#include "Math270A_DoubleArray2D.h" 

using namespace std; // Using the "standard" (std) standard template library components
#ifndef _RelaxOp2D_
#define _RelaxOp2D_






class RelaxOp2D
{
public:
	void initialize(double dt, double alphaX, double alphaY, const GridFun2D& f)
	{
		this->alphaX = alphaX;
		this->alphaY = alphaY;
		this->dt = dt;
		this->f = f;
		long M = f.xPanel;
		long N = f.yPanel;
		vector <double> loDiag(M);
		vector <double> diag(M+1);
		vector <double> upDiag(M);
		vector <double> loDiagY(N);
		vector <double> diagY(N+ 1);
		vector <double> upDiagY(N);
		double hx = f.hx;
		double hy = f.hy;

		//TrisolverX
		long i;

		// First equation is the identity

		i = 0;
		diag[0] = 1.0 ;
		upDiag[0] = 0.0;

		// Interior equation coefficients associated with a standard
		// second order finite difference approximation

		for (long i = 1; i < M; i++)
		{
			loDiag[i - 1] = -0.5*dt*alphaX / (hx*hx);
			upDiag[i] = -0.5*dt*alphaX / (hx*hx);
			diag[i] = dt*alphaX / (hx*hx) + 1;
		}


		// Last equation is the identity

		i = M;
		diag[M] = 1.0 ;
		loDiag[M - 1] = 0.0;

		triSolverX.initialize(M + 1, loDiag, diag, upDiag);



		//TrisolverY
		
		

		// First equation is the identity

		i = 0;
		diagY[0] = 1.0;
		upDiagY[0] = 0.0;

		// Interior equation coefficients associated with a standard
		// second order finite difference approximation

		for (long i = 1; i < N; i++)
		{
			loDiagY[i - 1] = -0.5*dt*alphaY / (hy*hy);
			upDiagY[i] = -0.5*dt*alphaY / (hy*hy);
			diagY[i] = dt*alphaY / (hy*hy) + 1;
		}


		// Last equation is the identity

		i = N;
		diagY[N] = 1.0;
		loDiagY[N - 1] = 0.0;
		triSolverY.initialize(N + 1, loDiagY, diagY, upDiagY);
	};

	void apply(const GridFun2D& uIn, GridFun2D& uOut)
	{
		//step 1 multiply by operator and subtract fdt
		GridFun2D ustar = uIn;
		vector<double> f1D_X(uIn.xPanel + 1);
		vector<double> u1D_X(uIn.xPanel + 1);

		vector<double> f1D_Y(uIn.yPanel + 1);
		vector<double> u1D_Y(uIn.yPanel + 1);

		double hx = uIn.hx; //space mesh
		double hy = uIn.hy;
		
			for(long j = 1; j < uIn.yPanel; j++)
			{
				// Interior equation coefficients associated with a standard
				// second order finite difference approximation
				for (long i = 1; i < uIn.xPanel; i++)
				{
					ustar.values(i, j) = uIn.values(i, j)*(1 - dt*(alphaX/(hx*hx) + 2*alphaY/(hy*hy))) + (alphaX/(hx*hx))*dt*0.5*(uIn.values(i + 1, j) + uIn.values(i - 1, j)) + (alphaY/(hy*hy))*dt*(uIn.values(i, j + 1) + uIn.values(i, j - 1))-dt*f.values(i,j);
				}
			}

		//step 2 TriSolverX
			for (long j= 1; j < uIn.yPanel; j++)
			{
				ustar.extractXslice(j, f1D_X);
				triSolverX.apply(f1D_X, u1D_X);
				ustar.insertXslice(j, u1D_X);
			}

		
		//step 3 get second intermideate
			for (long j = 1; j < uIn.yPanel; j++)
			{
				// Interior equation coefficients associated with a standard
				// second order finite difference approximation
				for (long i = 1; i < uIn.xPanel; i++)
				{
					ustar.values(i, j) -= 0.5*dt*(alphaY/(hy*hy))*(uIn.values(i, j + 1) + uIn.values(i, j - 1) - 2 * uIn.values(i, j));
				}
			}
				
				//step 4 TriSolverY
			for (long i = 1; i < uIn.xPanel; i++)
			{
				ustar.extractYslice(i, f1D_Y);
				triSolverY.apply(f1D_Y, u1D_Y);
				ustar.insertYslice(i, u1D_Y);
			}
			uOut = ustar;
		
		//cout << ukp1;
	};




	TriSolver triSolverX;
	TriSolver triSolverY;
	double alphaX;
	double alphaY;
	double dt;
	GridFun2D f;
};


#endif