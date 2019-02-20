/*
File name: matrix_multiplication.cpp
Date: 17 Jan 2017
Author: Andrew Roberts
*/

#include <iostream>
#include "gsl/gsl_matrix.h"
#include <iomanip>

using namespace std;

int main()
{
	int i,j;
	gsl_matrix * A = gsl_matrix_alloc(6,6);
	gsl_matrix * B = gsl_matrix_alloc(6,6);

	for(i=0; i<=5; ++i)
	{
		for(j=0; j<=5; ++j)
		{
			gsl_matrix_set(A,i,j,(i+1)*(j+1));
			gsl_matrix_set(B,i,j,(i+1)*(j+1));
		}
	}

	gsl_matrix_mul_elements(A,B);

    cout << "matrix A: \n";

    for(int i=0; i<=5; i++)
    {
        for(int j=0; j<=5; j++)
        {
            cout << setw(5) << gsl_matrix_get(A,i,j);
        }
        cout << "\n";
    }

	return 0;

}
