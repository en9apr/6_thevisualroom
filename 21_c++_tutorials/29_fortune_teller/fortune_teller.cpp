/*
File name: fortune_teller.cpp
Date: 13 Jan 2017
Author: Andrew Roberts
*/

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream> // reading, writing and creating files
#include <ctime>
using namespace std;

string getRandomReply(string [], int);

int main()
{
	ifstream inputfile;
	//declare an input file
	inputfile.open("replies.txt", ios::in);

	char answer[30];
	string answers[20];
	int pos = 0;

	//read from the file until end of file (eof)
	while(!inputfile.eof())
	{
		inputfile.getline(answer, 30);
		answers[pos] = answer;
		pos ++;
	}

	cout<<"Think of a question for the fortune teller, "
			"\npress enter for the answer "<<endl;
	cin.ignore();
	cout<<getRandomReply(answers,20)<<endl;

}
string getRandomReply(string replies[],int size)
{
	srand(time(0));
	int randomNum = rand()%20;
	return replies[randomNum];
}
