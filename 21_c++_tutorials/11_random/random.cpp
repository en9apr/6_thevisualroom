/*
File name: random.cpp
Date: 14 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;
int main()
{
	srand(time(0));
	for(int i = 0;i<10;i++)
	{
		cout<<rand()<<endl;
	}
	return 0;
}
