/*
File name: loan.cpp
Date: 11 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>

using namespace std;

int main()
{
	int months = 0;
	double balance = 5000;
	double payment = 150;

	while(balance > 0)
	{
		months=months+1; 	//or months++;
		balance = balance - payment;	//or balance -= payment
	}
	cout<<"It will take you "<<months<<" months to pay off your loan"<<endl;

}
