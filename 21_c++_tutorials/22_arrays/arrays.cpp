/*
File name: arrays.cpp
Date: 21 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

int main()
{
	//array declarations
	string seasons[4];
	int numbers[10];
	for(int i = 0;i<10;i++)
	{
		numbers[i] = rand();
	}

	//declare and initialize the array
	string weekdays[] = {"Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"};

	//change the value of an array
	seasons[0] = "summer";
	seasons[1] = "fall";
	seasons[2] = "winter";
	seasons[3] = "spring";

	cout<<seasons<<endl;
	cout<<seasons[0]<<endl;
	cout<<seasons[1]<<endl;
	cout<<seasons[2]<<endl;
	cout<<seasons[3]<<endl;

	cout<<numbers[0]<<endl;
	return 0;
}
