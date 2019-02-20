/*
File name: data_types.cpp
Date: 9 Dec 2016
Author: Andrew Roberts
*/

// allows cin, cout:
#include <iostream>

// standard namespaces:
using namespace std;

// main program:
int main()
{
	// parameters:
	char letter = 'a';
	short age = 10;
	int count = 575;
	long numStars = 985467528;
	float pi = 3.14;
	double price = 89.65;
	string season = "summer";

	// output:
	cout << "Letter: " << letter << endl;
	cout << "Age: " << age << endl;
	cout << "Count: " << count << endl;
	cout << "Number of stars in the sky: " << numStars << endl;
	cout << "Pi: " << pi << endl;
	cout << "Price: " << price << endl;
	cout << "Season: " << season << endl;

	// return integer:
	return 0;
}
