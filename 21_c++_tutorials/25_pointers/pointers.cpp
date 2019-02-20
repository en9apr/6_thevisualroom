/*
File name: pointers.cpp
Date: 24 Dec 2016
Author: Andrew Roberts
*/

#include <iostream>

using namespace std;

int main()
{
/*
	int JohnsAcct = 12345;
	int* acctPtr = &JohnsAcct;

	cout<<"Johns acct number: "<<JohnsAcct<<" is at: "<<acctPtr<<endl;
	cout<<"acctPtr points to this value: *acctPtr:"<<*acctPtr<<endl;
*/
	//another example
	cout<<endl<<"Next example, declares num = 95:"<<endl;
	int num = 95;
	cout<<"Address to find num value: "<<&num<<endl;
	cout<<"size of num value based on datatype: "<<sizeof(num)<<endl;
	cout<<"value at this location: "<<num<<endl;

	cout<<endl<<"Next example, add numPtr:"<<endl;
	//lets add a pointer to hold the address of num
	int* numPtr = &num;

	//lets change the value of num using the pointer
	cout<<"add 5 to num using the pointer and not the variable:"<<endl;
	*numPtr += 5;  //num value is now 100
	cout<<"Num is now: "<<num<<endl;
	cout<<"Num location is still the same: "<<&num<<endl;
	cout<<"Location of numPtr is: "<<numPtr<<endl;
	cout<<"Value of numPtr is: "<<*numPtr<<endl;

	cout<<endl<<"Next example, using characters"<<endl;
	char letter = 'A';
	cout<<"size of letter value based on datatype: "<<
			sizeof(letter)<<endl;
	cout<<"value at this location: "<<letter<<endl;
	//lets add a pointer to hold the address of letter
	char* letterPtr = &letter;

	//lets change the value of letter using the pointer
	cout<<"\nadd 25 to letter:"<<endl;
	*letterPtr += 25;  //letter value is now Z
	cout<<"Letter is now: "<<letter<<endl;

	return 0;
}



