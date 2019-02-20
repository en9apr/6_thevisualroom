/*
File name: challenge_file_read.cpp
Date: 14 Jan 2017
Author: Andrew Roberts
*/

#include <iostream>
#include <fstream> // reading, writing and creating files
#include <iomanip> // setprecision
#include <vector>
#include <stdlib.h>
#include <numeric>
using namespace std;

int add_columns(vector<int>, vector<int>);
double average_column(vector<int>);
int error_code=666;

int main()
{
	ifstream inputfile;
	try
	{
		inputfile.open("numbers.txt", ios::in);
		if(inputfile.fail())
		{
			throw error_code;
		}
	}
	catch(int error_number)
	{
		if(error_number==error_code)
		{
			cout << "File reading error." << endl;
			exit(1);
		}
	}

	int number_1;
	int number_2;
	vector<int> column_1;
	vector<int> column_2;

	while(!inputfile.eof())
	{
		inputfile >> number_1;
		column_1.push_back(number_1);
		inputfile >> number_2;
		column_2.push_back(number_2);
	}

	int total = add_columns(column_1, column_2);
	cout << "Total: " << total << endl;

	double average_1 = average_column(column_1);
	cout << "Average column 1: " << fixed << setprecision(2) << average_1 << endl;

	double average_2 = average_column(column_2);
	cout << "Average column 2: " << fixed << setprecision(2) << average_2 << endl;
}

int add_columns(vector<int> column_1, vector<int> column_2)
{
	double total_1 = accumulate(column_1.begin(), column_1.end(), 0.0);
	double total_2 = accumulate(column_2.begin(), column_2.end(), 0.0);
	return total_1+total_2;
}

double average_column(vector<int> column_1)
{
	double average = accumulate(column_1.begin(), column_1.end(), 0.0)/column_1.size();
	return average;
}

