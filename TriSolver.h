#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
using namespace std; // Using the "standard" (std) standard 

#ifndef _TriSolver_
#define _TriSolver_


class TriSolver
{
public:
	TriSolver()
	{
		initialize();
	};


	TriSolver(const TriSolver& T) // called when you declare an instance with an existing instance
	{
		initialize(T);
	};

	TriSolver(long systemSize, vector<double>& loDiag, vector<double>& diag,
		vector<double>& upDiag)
	{
		initialize(systemSize, loDiag, diag, upDiag);
	}

	virtual ~TriSolver() {};
	// Do nothing initializer

	void initialize(long systemSize, vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
	{
		this->systemSize = systemSize;
		this->loDiag = loDiag;
		this->diag = diag;
		this->upDiag = upDiag;
	}

	void initialize()
	{
		systemSize = 0;
		loDiag.clear();
		diag.clear();
		upDiag.clear();
	}


	void initialize(const TriSolver& T)
	{
		systemSize = T.systemSize;
		loDiag = T.loDiag;
		diag = T.diag;
		upDiag = T.upDiag;
	};

	// Thomas algorithm as in https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

	void apply(vector<double>& f, vector<double>& u)
	{

		u.resize(systemSize);
		vector<double> cp(systemSize);
		if (diag[0] == 0)
		{
			cout << "0 Error" << endl;
			return;
		}
		cp[0] = upDiag[0] / (diag[0]);
		for (long i = 1; i < systemSize - 1; i++)
		{
			double denom = diag[i] - loDiag[i - 1] * cp[i - 1];
			if (denom == 0)
			{
				cout << "0 Error" << endl;
				return;
			}

			cp[i] = upDiag[i] / (denom);

		}

		vector<double> dp(systemSize);
		dp[0] = f[0] / (diag[0]);
		for (long i = 1; i < systemSize; i++)
		{
			double denom = diag[i] - loDiag[i - 1] * cp[i - 1];
			if (denom == 0)
			{
				cout << "0 Error" << endl;
				return;
			}

			dp[i] = (f[i] - loDiag[i - 1] * dp[i - 1]) / (denom);
		}

		u[systemSize - 1] = dp[systemSize - 1];

		for (long i = systemSize - 2; i > -1; i--)
		{
			u[i] = dp[i] - cp[i] * u[i + 1];
		}
		return;
	}

	friend ostream& operator<<(ostream& outStream, const TriSolver& V)
	{
		for (long i = 0; i < V.diag.size(); i++)
		{
			outStream << setw(5) << V.diag[i] << " ";
			outStream << endl;
		}
		return outStream;
	};


	long        systemSize;
	vector<double>  loDiag;
	vector<double>    diag;
	vector<double>  upDiag;
};
#endif // !_TriSolver_
