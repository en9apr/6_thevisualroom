/*
File name: bool_function.cpp
Date: 16 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

/*
 * display a message
 */
int number;
bool guess(int);

int main()
{
	int usernum;
	srand(time(0));
	number = rand() % 50 + 1;
	cout<<"Guess a number between 1 and 50:"<<endl;
	cin>>usernum;
	while(guess(usernum) != true)
	{
		cout<<"Try again:"<<endl;
		cin>>usernum;
	}
	cout<<"You guessed it!!! "<<endl;
	return 0;
}
bool guess(int num)
{
	if(num == number)
		return true;
	if(num > number)
		cout<<"Too high"<<endl;
	else
		cout<<"Too low"<<endl;
	return false;
}

