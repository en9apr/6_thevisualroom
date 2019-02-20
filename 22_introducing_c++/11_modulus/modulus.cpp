/*
File name: modulus.cpp
Date: 22 Jan 2017
Author: Andrew Roberts
*/

#include <iostream>

using namespace std;

int main()
{

	int i,j,m,n,a,b;
	cout << "Enter integer: ";
	cin >> i;
	cout << "Enter integer: ";
	cin >> j;
	m=i/j;
	n=i%j;
	a=m*j;
	b=a+n;
	cout << "i= " << i << " j= " << j << "\n\n"
		 << "i/j= " << m << " i%j= " << n << "\n\n"
		 << "The following result should be " << i << ":\n"
		 << "(i/j)*j+i%j= " << b << "\n";
	return 0;

}
