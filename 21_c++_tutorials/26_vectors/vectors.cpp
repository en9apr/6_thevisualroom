/*
File name: vectors.cpp
Date: 24 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <vector>
#include <cstdlib>
using namespace std;

int main()
{
	//create three vectors
	//carVins contains the vin of each car on a showroom floor
	//(note: carVins starts as an empty vector)
	vector<int> carVins;

	//num contains 10 integer values
	vector<int> num(10);

	//prices contains 5 prices that are initially zero
	vector<double> prices(5,0.0);

	int numCars;
	int vin;

	cout<<"How many cars are in the showroom?\n";
	cin>>numCars;
	for(int i = 0; i<numCars; i++)
	{
		cout<<"Enter the vin for car "<<i+1<<" :";
		cin>>vin;
		//add the new vin to the end of the vector
		carVins.push_back(vin);
	}

	//print the vins:
	cout<<"\nVin numbers for cars in showroom:\n";
	for(int i = 0; i<carVins.size();i++)
	{
		cout<<carVins[i]<<endl;
	}

	//add 10 random integers to the num vector
	int randomNum;
	for(int i = 0; i < num.size(); i++)
	{
		num[i] = rand();
	}

	//print the random numbers:
	cout<<"\n10 random integers: \n";
	for(int i = 0; i<num.size();i++)
	{
		cout<<num[i]<<endl;
	}

	//Other vector functions include:
	//.at returns the value located at element 0
	int value = num.at(0);

	//.capacity() returns the max number of elements
	//that may be stored
	//without adding additional memory
	int value2 = num.capacity();

	//.pop_back removes the last element from a vector
	num.pop_back();

	//.clear() is used to clear a vector of all elements
	num.clear();

	//.empty() returns true if vector is empty
	if(num.empty()==true)
		cout<<"vector num is empty\n";

	//.swap(vector2) swaps the contents

}
