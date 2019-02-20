/*
File name: overflow_underflow.cpp
Date: 22 Jan 2017
Author: Andrew Roberts
*/

#include <iostream>
#include <climits>

using namespace std;

int main()
{
	int i;
	//i = INT_MAX +10;
	i = INT_MAX;
	cout << "i= " << i << "\n";
	//i = INT_MIN -10;
	i = INT_MIN;
	cout << "i= " << i << "\n";

	return 0;

}
