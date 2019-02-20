/*
File name: Movie.cpp
Date: 18 Dec 2016
Author: Andrew Roberts
*/

#include "Movie.h"
#include <iostream>

using namespace std;

// default constructor:
Movie::Movie()
{
	name = "Star Wars";
	year = 1977;
	rating = "PG";
}

// constructor:
Movie::Movie(string name, int year, string rating)
{
	this->name = name;
	this->year = year;
	this->rating = rating;
}

// destructor:
Movie::~Movie()
{

}

// accessor methods

string Movie::get_name()
{
	return name;
}

int Movie::get_year()
{
	return year;
}

string Movie::get_rating()
{
	return rating;
}

// mutator methods

void Movie::set_name(string name)
{
	this->name = name;
}

void Movie::set_year(int year)
{
	while(year < 1950 || year > 2016)
	{
		cout << "Sorry, please enter movie year between 1950 and 2016: ";
		cin >> year;
		cin.ignore(); // needed so the getline function does not use the return as input
	}
	this->year = year;
}

void Movie::set_rating(string rating)
{
	this->rating = rating;
}
