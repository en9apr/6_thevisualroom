/*
File name: mailing_addresses.cpp
Date: 13 Jan 2017
Author: Andrew Roberts
*/

#include <iostream>
#include <fstream>
using namespace std;

//create a new datatype for holding addresses
struct Address
{
	string fullName;
	char houseNum[5];
	char streetName[30];
	char city[30];
	char state[3];
	int zipcode;
};

void writeToFile(Address m[], int size);

int main()
{
	int numees;
	string fullname;
	string fname, lname;
	cout<<"How many employees do you have? \n";
	cin>>numees;
	Address mailing[numees];

	for(int i = 0;i < numees;i++)
	{
		cout<<"Please enter first name: \n";
		cin>>fname;
		cout<<"Enter last name: \n";
		cin>>lname;
		fullname = fname + " " + lname;
		mailing[i].fullName = fullname;
		cout<<"Enter house number: \n";
		cin>>mailing[i].houseNum;
		cin.ignore();
		cout<<"Enter street name: \n";
		cin.getline(mailing[i].streetName, 30);
		cout<<"Enter city: \n";
		cin>>mailing[i].city;
		cout<<"Enter state: \n";
		cin>>mailing[i].state;
		cout<<"Enter zipcode: \n";
		cin>>mailing[i].zipcode;
	}

	writeToFile(mailing,numees);

}

void writeToFile(Address m[], int size)
{
	ofstream output;
	output.open("mailingaddress.txt", ios::out);
	for(int i = 0; i < size; i++)
	{
		output<<m[i].fullName<<endl;
		output<<m[i].houseNum<<" "<<m[i].streetName<<endl;
		output<<m[i].city<<", "<<m[i].state<<" " <<m[i].zipcode<<endl;
		output<<endl;

	}
}
