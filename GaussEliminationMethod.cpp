
/*
Program: Gauss Elimination Method
All array indexes are assumed to start from 1
*/

#include"GaussEliminationMethod.h"

using namespace std;

bool gaussEliminationMethod(int n, double** AB, double* X)
{
	int i, j, k;
	double m, s;
	double eps = 0.000001;
	// elimination of coefficients

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// calculating unknowns

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}