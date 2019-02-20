/*
File name: Movie.h
Date: 18 Dec 2016
Author: Andrew Roberts
*/

#ifndef MOVIE_H_
#define MOVIE_H_

#include <string>

using namespace std;

class Movie
{
	private:
		string name;
		int year;
		string rating;

	public:
		// constructors
		Movie();
		Movie(string, int, string);

		// destructor
		~Movie();

		string get_name();
		int get_year();
		string get_rating();

		void set_name(string);
		void set_year(int);
		void set_rating(string);
};

#endif /* MOVIE_H_ */
