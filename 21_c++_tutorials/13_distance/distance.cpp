/*
File name: distance.cpp
Date: 16 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <cmath>
using namespace std;

double getDistance(int x1, int y1, int x2, int y2)
{
	double distance;
	distance = sqrt(pow((x2-x1),2) + pow((y2 - y1), 2));
	return distance;
}

int main()
{
	int x1, x2, y1, y2;

	cout<<"Please enter the first point x1, then y1"<<endl;
	cin>>x1;
	cin>>y1;

	cout<<"Now enter the second point"<<endl;
	cin>>x2;
	cin>>y2;

	cout<<"The distance between these two points is: "<<
			getDistance(x1,y1,x2,y2)<<endl;

	return 0;
}
