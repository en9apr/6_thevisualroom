/*
File name: dollars_and_cents.cpp
Date: 9 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <cmath>

using namespace std;

int main()
{

	// variables:
	int dollars, cents;
	double price, convert;

	// implicit conversion:
	char dollar_sign = 36;

	// input price:
	cout << "Enter price in dollars and cents: " << endl;
	cin >> price;

	// round down price to dollars:
	dollars = static_cast<int>(floor(price));

	// convert price to cents and round cents up:
	convert = price * 100 - static_cast<double>(dollars * 100);
	cents = static_cast<int>(ceil(convert));

	// output to screen:
	cout << "Total dollars: " << dollar_sign << dollars << " cents: " << cents << endl;

	return 0;

}
