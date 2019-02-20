/*
File name: challenge_change.cpp
Date: 16 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// function prototypes
void get_items(string&);
void get_price_and_paid(int&, int&, int&, int&);
void set_change(int&, int, int, int, int);
void set_coins(int&, int&, int&, int&, int&, int&);
void print_coins(int&, int&, int&, int&, int&, int&);

int main()
{
	int price_dollars, price_cents, paid_dollars, paid_cents;
	int dollars, quarters, dimes, nickels, pennies;
	string items;
	int change = 0; // initialises change to zero

	cout << "-----------------------" << endl;
	cout << "Welcome to the store!!!" << endl;
	cout << "-----------------------" << endl;

	get_items(items);

	while (items == "y")
	{
		get_price_and_paid(price_dollars, price_cents, paid_dollars, paid_cents);
		set_change(change, price_dollars, price_cents, paid_dollars, paid_cents);

		while (change < 0)
		{
			cout << "I'm sorry, the price is wrong, or the customer has not paid the full amount!" << endl;
			get_price_and_paid(price_dollars, price_cents, paid_dollars, paid_cents);
			set_change(change, price_dollars, price_cents, paid_dollars, paid_cents);
		}

		set_coins(change, dollars, quarters, dimes, nickels, pennies);
		print_coins(change, dollars, quarters, dimes, nickels, pennies);
		get_items(items);
	}

	cout << "-------------------------------" << endl;
	cout << "Thank you, please come again!!!" << endl;
	cout << "-------------------------------" << endl;
	return 0;
}

void get_items(string& items)
{
	cout << "Do you have any items to pay for? (y or n):" << endl;
	cin >> items;

	while (items != "y" && items != "n")
	{
		cout << "I'm sorry, please enter y or n if you have items:" << endl;
		cin >> items;
	}
}

void get_price_and_paid(int& price_dollars, int& price_cents, int& paid_dollars, int& paid_cents)
{
	cout << "Please enter the price of the item (in dollars): " << endl;
	cin >> price_dollars;
	cout << "Please enter the price of the item (in cents): " << endl;
	cin >> price_cents;
	cout << "Please enter the amount the customer paid (in dollars):" << endl;
	cin >> paid_dollars;
	cout << "Please enter the amount the customer paid (in cents):" << endl;
	cin >> paid_cents;
}

void set_change(int& change, int price_dollars, int price_cents, int paid_dollars, int paid_cents)
{
	change = 100*(paid_dollars-price_dollars) + (paid_cents-price_cents);
}

void set_coins(int& change, int& dollars, int& quarters, int& dimes, int& nickels, int& pennies)
{
	cout << "The change in pennies is: " << change << endl;
	dollars = change / 100;
	quarters = (change % 100) / 25;
	dimes = (change % 100 % 25) / 10;
	nickels = (change % 100 % 25 % 10) / 5;
	pennies = (change % 100 % 25 % 10 % 5) / 1;
	change = change - (dollars * 100) - (quarters * 25) - (dimes * 10) - (nickels * 5) - (pennies * 1);
}

void print_coins(int& change, int& dollars, int& quarters, int& dimes, int& nickels, int& pennies)
{
	cout << "The number of dollars in the change is: " << dollars << endl;
	cout << "The number of quarters in the change is: " << quarters << endl;
	cout << "The number of dimes in the change is: " << dimes << endl;
	cout << "The number of nickels in the change is: " << nickels << endl;
	cout << "The number of pennies in the change is: " << pennies << endl;
	cout << "The change in pennies is now: " << change << endl;
}
