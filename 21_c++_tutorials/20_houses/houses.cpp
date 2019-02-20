/*
File name: houses.cpp
Date: 17 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include "House.h"
#include <string>
using namespace std;

void printHouse(House);

int main()
{
	//Create an instance of the House class
	//using the default constructor
	House house1;

	//Create a second instance
	//provide values for color, bed, bath, and sqft
	House house2("Tan", 3, 2, 2500);

	//print house information
	printHouse(house1);
	printHouse(house2);
	house2.setColor("White");
	house2.setNumBath(3);
	house2.setNumBed(4);
	house2.setSqft(3500);
	printHouse(house2);

	return 0;
}
void printHouse(House house)
{
	//print the house information using
	//the dot operator
	cout<<"House Color: "<<house.getColor()
			<<"\nNumber Bathrooms: "
			<<house.getNumBath()
			<<"\nNumberBedrooms: "
			<<house.getNumBed()<<endl;
	cout<<"Total Square feet: "
			<<house.getSqft()<<endl;
	cout<<endl<<endl;

	//since it is a void function,
	//I do not need a return statement
}
