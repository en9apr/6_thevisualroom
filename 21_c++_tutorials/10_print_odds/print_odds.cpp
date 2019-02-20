/*
File name: print_odds.cpp
Date: 11 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>

using namespace std;

int main()
{
	for(int num = 1; num <= 100; num++)
	{
		if(num %2 != 0)
			cout<<num<<endl;
	}
	return 0;

}
