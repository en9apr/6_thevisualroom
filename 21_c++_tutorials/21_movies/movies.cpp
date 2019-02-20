/*
File name: movies.cpp
Date: 18 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>
#include "Movie.h"
#include <string>
#include <cctype>
using namespace std;

void print_movie(Movie);
void set_movie(Movie&);

int main()
{
	cout << "Please enter your three favourite movies: " << endl;

	Movie movie_one;
	set_movie(movie_one);
	Movie movie_two;
	set_movie(movie_two);
	Movie movie_three;
	set_movie(movie_three);

	print_movie(movie_one);
	print_movie(movie_two);
	print_movie(movie_three);

	return 0;
}

void print_movie(Movie movie)
{
	cout << "Movie name: " << movie.get_name() << endl;
	cout << "Movie year: " << movie.get_year() << endl;
	cout << "Movie rating: " << movie.get_rating() << endl;
	cout << "" << endl;
}

void set_movie(Movie& new_movie)
{
	string new_name;
	int new_year;
	string new_rating;

	cout << "Enter movie name: ";
	getline(cin, new_name);
	new_movie.set_name(new_name);

	cout << "Enter movie year: ";
	cin >> new_year;
	cin.ignore(); // needed so the getline function does not use the return as input
	new_movie.set_year(new_year);

	cout << "Enter movie rating: ";
	getline(cin, new_rating);
	new_movie.set_rating(new_rating);

	cout << "" << endl;
}




