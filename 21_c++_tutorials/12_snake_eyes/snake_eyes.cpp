/*
File name: snake_eyes.cpp
Date: 14 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
using namespace std;
int main()
{
	srand(time(0));
	int count = 0;
	for (int i = 0;i<100;i++)
	{
		int die1 = rand()%6 + 1;
		int die2 = rand()%6 + 1;
		cout<<setw(6)<<die1<<setw(6)<<die2;
		if(die1 == 1 && die2 == 1)
		{
			cout<<"  Snake Eyes";
			count ++;
		}
		cout<<endl;
	}
	cout<<"There were "<<count<<" snake eyes"<<endl;
}
