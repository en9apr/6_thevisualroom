/*
File name: employees.cpp
Date: 23 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>

using namespace std;

//define the new data structure
struct Address
{
	char houseNum[5];
	char streetName[30];
	char city[30];
	char state[3];
	int zipcode;
};
struct DateOfBirth
{
	int month;
	int day;
	int year;

};
struct Employee
{
	char firstName[30];
	char lastName[30];
	int eeNumber;
	double hourlyWage;
	Address homeAddress;
	DateOfBirth dob;
};
//printEes prototype
void printEes(Employee[],int);

int main()
{
	int numEmployees;
	cout<<"How many employees do you have? \n";
	cin>>numEmployees;
	//create an array of type Employee
	Employee employees[numEmployees];

	//loop for creating new employee data structure variables
	for(int i = 0;i < numEmployees;i++)
	{
		cout<<"Please enter first name: \n";
		cin>>employees[i].firstName;
		cout<<"Enter last name: \n";
		cin>>employees[i].lastName;
		cout<<"Enter employee id: \n";
		cin>>employees[i].eeNumber;
		cout<<"Enter hourly wage: \n";
		cin>>employees[i].hourlyWage;
		cout<<"Enter house number: \n";
		cin>>employees[i].homeAddress.houseNum;
		cin.ignore();
		cout<<"Enter street name: \n";
		cin.getline(employees[i].homeAddress.streetName, 30);
		cout<<"Enter city: \n";
		cin>>employees[i].homeAddress.city;
		cout<<"Enter state: \n";
		cin>>employees[i].homeAddress.state;
		cout<<"Enter zipcode: \n";
		cin>>employees[i].homeAddress.zipcode;
		cout<<"Enter birth month, day and year:  \n";
		cin>>employees[i].dob.month;
		cin>>employees[i].dob.day;
		cin>>employees[i].dob.year;
	}

	printEes(employees,numEmployees);
}

void printEes(Employee employees[], int size)
{
	for(int i = 0; i < size; i++)
	{
		cout<<employees[i].firstName <<" "
				<<employees[i].lastName<<endl;
		cout<<employees[i].homeAddress.houseNum<<" "
				<<employees[i].homeAddress.streetName<<endl;
		cout<<employees[i].homeAddress.city<<", "
				<<employees[i].homeAddress.state<<" "
				<<employees[i].homeAddress.zipcode<<endl;
		cout<<"Hourly wage: "<<employees[i].hourlyWage;
		cout<<endl;
		cout<<"Date of Birth: ";
		switch(employees[i].dob.month)
		{
		case 1:
			cout<<"January ";
			break;
		case 2:
			cout<<"February ";
			break;
		case 3:
			cout<<"March ";
			break;
		case 4:
			cout<<"April ";
			break;
		case 5:
			cout<<"May ";
			break;
		case 6:
			cout<<"June ";
			break;
		case 7:
			cout<<"July ";
			break;
		case 8:
			cout<<"August ";
			break;
		case 9:
			cout<<"September ";
			break;
		case 10:
			cout<<"October ";
			break;
		case 11:
			cout<<"November ";
			break;
		case 12:
			cout<<"December ";
			break;
		default:
			cout<<"Invalid Month";
			break;
		};
		cout<<employees[i].dob.day<<", ";
		cout<<employees[i].dob.year<<endl;
	}
}
