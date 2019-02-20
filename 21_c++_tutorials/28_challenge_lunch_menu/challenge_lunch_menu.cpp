/*
File name: challenge_lunch_menu.cpp
Date: 24 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <string>
//#include <algorithm>
#include <vector>
// to allow setprecision:
#include <iomanip>
using namespace std;

//int findNumber(string array[], int size, string target);

int search(string[], int, string);
void print_total_cost(int[], double[], int[], double);
void get_choices(int, string[], int[], int, string[]);
void get_quantities(int, string[], int[]);

int main()
{
	string menu[] = {"hotdog", "hamburger", "french fries", "salad", "exit"};
	double prices[] = {1.20, 3.34, 0.99, 1.11, 0.00};
	string choices[] = {"exit", "exit", "exit", "exit"};
	int quantities[] = {0, 0, 0, 0};
	int positions[] = {4, 4, 4, 4};
	int arraysize = 5;
	double cost = 0.0;

	for(int i=0; i<4; i++)
	{
		get_choices(i, choices, positions, arraysize, menu);

		while(positions[i] == -1)
		{
			cout << "Invalid choice." << endl;
			get_choices(i, choices, positions, arraysize, menu);
		}

		if (choices[i].compare("exit") == 0)
		{
			break;
		}

		get_quantities(i, choices, quantities);
	}

	print_total_cost(quantities, prices, positions, cost);

	return 0;

}

void get_choices(int i, string choices[], int positions[], int arraysize, string menu[])
{
	cout << "Please enter choice number " << i+1 << " (hotdog, hamburger, french fries, salad or exit)." << endl;
	getline(cin, choices[i]);
	positions[i] = search(menu, arraysize, choices[i]);
}

void get_quantities(int i, string choices[], int quantities[])
{
	cout << "Please enter quantity of " << choices[i] << endl;
	cin >> quantities[i];
	cin.ignore();
}

void print_total_cost(int quantities[], double prices[], int positions[], double cost)
{
	cout << "This is the total cost:" << endl;
	for(int i=0; i<4; i++)
	{
		cost += quantities[i]*prices[positions[i]];
	}
	cout << "$" << fixed << setprecision(2) << cost << endl;
}

int search(string numbers[], int size, string value)
{

	bool exit = false;
	int pos = 0;
	while(!exit)
	{
		if(numbers[pos] == value)
		{
			exit = true;
		}
		else if(pos + 1 == size)
		{
			pos = -1;
			exit = true;
		}
		else
		{
			pos = pos + 1;
		}
	}
	return pos;
}

