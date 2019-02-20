/*
File name: quadratic_equation.cpp
Date: 17 Jan 2017
Author: Andrew Roberts
*/

#include <iostream>
#include <cmath>

using namespace std;

int main()
{

	double root1, root2, a, b, c, root;

	cout << "Enter the coefficients a, b and c \n";
	cin >> a >> b >> c;

	root = sqrt(b*b - 4.0*a*c);
	root1 = 0.5*(root - b)/a;
	root2 = -0.5*(root + b)/a;

	cout << "The roots are " << root1 << " and " << root2 << "\n";

	return 0;
}
