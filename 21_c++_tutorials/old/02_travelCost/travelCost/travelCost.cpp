/*
travelCost.cpp
8 December 2016
Andrew Roberts
*/

// to allow cin, cout, fixed:
#include <iostream>

// to allow setprecision:
#include <iomanip>

// standard namespaces:
using namespace std;

// main function:
int main()
{
	// variables:
	double mpg, distance, gallons, pricePerGallon, totalCost;

	// input mpg:
	cout << "Enter mpg" << endl;
	cin >> mpg;

	// input distance:
	cout << "Enter distance in miles: " << endl;
	cin >> distance;

	// input price per gallon:
	cout << "Enter price for one gallon of gas: " << endl;
	cin >> pricePerGallon;

	// calculate gallons needed:
	gallons = distance/mpg;

	// calculate total cost:
	totalCost = gallons * pricePerGallon;

	// output total cost to screen:
	// fixed point notation, with a precision of 2, e.g. 1000.00 not 1e3:
	cout << "Total trip cost: $" << fixed << setprecision(2) << totalCost << endl;

	return 0;
}
