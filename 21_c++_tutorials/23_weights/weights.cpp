/*
File name: heights.cpp
Date: 21 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <string>
using namespace std;

void getWeights(double[], int);
double findAverage(double[], int);
double findHeaviest(double[], int);

int main()
{
	double weights[10];
	getWeights(weights, 10);
	cout<<"The average weight of all students is: "
			<<findAverage(weights, 10)<<" pounds"<<endl;
	cout<<"The heaviest student is: "
				<<findHeaviest(weights, 10)<<" pounds"<<endl;

	return 0;
}
void getWeights(double h[], int size)
{
	for(int i = 0; i<size; i++)
	{
		cout<<"Enter student weight in pounds: \n";
		cin>>h[i];
	}
}
double findAverage(double h[], int size)
{
	int count = 0;
	double total = 0;
	for(int i = 0; i<size; i++)
		{
			count++;
			total += h[i];
		}
	return total/count;
}
double findHeaviest(double h[], int size)
{
	double heaviest = h[0];
	for(int i = 1; i<size; i++)
	{
		if(h[i] > heaviest)
		{
			heaviest = h[i];
		}
	}
	return heaviest;
}

