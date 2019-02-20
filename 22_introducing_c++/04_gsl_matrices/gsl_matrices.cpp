#include <iostream>
#include <iomanip>
#include <gsl/gsl_matrix.h>

using namespace std;

int main()
{

    gsl_matrix * a = gsl_matrix_alloc(2,2);
    gsl_matrix * b = gsl_matrix_alloc(2,2);

    gsl_matrix_set_zero(a);
    gsl_matrix_set_zero(b);

    for(int i=0; i<=1; i++)
    {
        for(int j=0; j<=1; j++)
        {
            gsl_matrix_set(a,i,j,1.0);
            gsl_matrix_set(b,i,j,2.0);
        }
    }

    gsl_matrix_add(a,b);

    cout << "matrix a: \n";

    for(int i=0; i<=1; i++)
    {
        for(int j=0; j<=1; j++)
        {
            cout << setw(3) << gsl_matrix_get(a,i,j);
        }
        cout << "\n";
    }

    cout << "matrix b: \n";

    for(int i=0; i<=1; i++)
    {
        for(int j=0; j<=1; j++)
        {
            cout << setw(3) << gsl_matrix_get(b,i,j);
        }
        cout << "\n";
    }

    gsl_matrix_free(a);
    gsl_matrix_free(b);

    return 0;
}
