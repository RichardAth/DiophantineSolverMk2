// Diophantinechecker.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"

/* get a number from the input stream */
long long getnumber(std::string msg) {
	long long result;

	std::cout << msg;
	while (1) {
		std::cin >> result;     // as result is an integer only numeric input is allowed
		if (std::cin.good()) {
			break;			// valid number
		}
		else {
			// not a valid number
			std::cout << "Invalid Input! Please input a numerical value.  ";
			std::cin.clear();               // clear error flags
			std::cin.ignore(100000, '\n');  // discard invalid input up to newline/enter
		}
	}
	return result;
}

int main(int argc, char* argv[]) {
	long long int a, b, c, d, e, f, x, y, value;
	long long int p, q, r, s, k, l, x1, y1;
	char yn = '\0';
	bool test = false;

	std::cout << "Check Diophantine equations of the form: Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0" << "\n";
	a = getnumber("Enter value for A ");
	b = getnumber("Enter value for B ");
	c = getnumber("Enter value for C ");
	d = getnumber("Enter value for D ");
	e = getnumber("Enter value for E ");
	f = getnumber("Enter value for F ");
	x = getnumber("Enter value for x ");
	y = getnumber("Enter value for y ");
	value = a*x*x + b*x*y + c*y*y + d *x + e*y + f;

	if (value == 0)
		std::cout << "x, y values are valid\n";
	else {
		std::cout << "** value not zero; it is " << value << "\n";
		system("PAUSE");
		return 0;
	}

	while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
		std::cout << "Enter P, Q, K, R, S, L (enter 0 for unused values)? (Y/N): ";
		std::cin >> yn;
		test = (toupper(yn) == 'Y');
	}

	if (test) {
		p = getnumber("Enter value for P ");
		q = getnumber("Enter value for Q ");
		k = getnumber("Enter value for K ");
		r = getnumber("Enter value for R ");
		s = getnumber("Enter value for S ");
		l = getnumber("Enter value for L ");

	/*  X(n+1) = P X(n) + Q Y(n)        - if k and l are zero
        Y(n+1) = R X(n) + S Y(n)  
		
		it follows that Y(n-1) = (X(n)/P - Y(n)/R) / (Q/P - S/R)
		and             X(n-1) = (X(n)/Q - Y(n)/S) / (P/Q - R/S)
		provided that these are integers and don't cause division by zero */
		x1 = p*x + q*y + k;
		y1 = r*x + s*y + l;
		value = a*x1*x1 + b*x1*y1 + c*y1*y1 + d*x1 + e*y1 + f;
		if (value == 0)
			std::cout << " P, Q, K, R, S, L values are valid" << " x1=" << x1 << " y1=" << y1 << "\n";
		else {
			std::cout << "** value not zero; it is " << value << " x1=" << x1 << " y1=" << y1 << "\n";
		}
		system("PAUSE");
	}
    return 0;
}

