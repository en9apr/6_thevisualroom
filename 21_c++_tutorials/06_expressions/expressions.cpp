/*
File name: expressions.cpp
Date: 11 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>

using namespace std;

int main()
{
	// variables:
	int length, width, area, perimeter;
	int radius, circleArea, circumference;
	string firstName, lastName, fullName;

	// area and perimeter rectangle:
	cout<<"Enter the length and width of a rectangle"<<endl;
	cin>>length;
	cin>>width;
	area = length*width;
	perimeter = 2*length + 2*width;
	cout<<"The area of the rectangle is: "<<area<<endl;
	cout<<"The perimeter of the rectangle is: "<<perimeter<<endl;

	// area and circumference of circle:
	cout<<"Enter the radius of a circle:"<<endl;
	cin>>radius;
	circleArea = 3.14 * radius * radius;
	circumference = 2 * 3.14 * radius;
	cout<<"The area of the circle is: "<<circleArea<<endl;
	cout<<"The circumference of the circle is: "<<circumference<<endl;

	// add two strings:
	cout<<"Enter the first and last name"<<endl;
	cin>>firstName;
	cin>>lastName;
	fullName = firstName + " " +lastName;
	cout<<"The full name is: "<<fullName<<endl;
	return 0;

}
