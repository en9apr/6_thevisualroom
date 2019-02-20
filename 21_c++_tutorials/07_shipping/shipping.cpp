/*
File name: shipping.cpp
Date: 11 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
	double price, shipping;

	cout<<"enter the total price of the order: "<<endl;
	cin>>price;

	if(price > 75)
		shipping = 0;
	else if(price > 50)
		shipping = 5;
	else if(price > 25)
		shipping = 10;
	else if(price > 0)
		shipping = 15;

	cout<<"Total price including shipping is: "<<fixed<<setprecision(2)<<price + shipping<<endl;

	return 0;

}
