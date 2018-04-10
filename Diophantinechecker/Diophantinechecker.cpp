// Diophantinechecker.cpp : Defines the entry point for the console application.
// This is a little program that can take the output from diophantine solver as input 
// and check that the values are correct.

#include <boost/multiprecision/gmp.hpp>
#include <io.h>
#include "stdafx.h"
using namespace boost::multiprecision;
bool IpFromFile = false;

/* get a number from the input stream */
mpz_int getbignumber(std::string msg) {
	mpz_int result;
	if (!IpFromFile)
		std::cout << msg;
	while (1) {
		std::cin >> result;     // as result is an integer only numeric input is allowed
		if (std::cin.good()) {
			break;			// valid number
		}
		else {
			// not a valid number
			std::cout << "Invalid Input! Please input a numerical value.  ";
			Beep(1200, 500);
			std::cin.clear();               // clear error flags
			std::cin.ignore(100000, '\n');  // discard invalid input up to newline/enter
		}
	}
	return result;
}
long long getnumber(std::string msg) {
	long long result;
	if (!IpFromFile)
		std::cout << msg;
	while (1) {
		std::cin >> result;     // as result is an integer only numeric input is allowed
		if (std::cin.good()) {
			break;			// valid number
		}
		else {
			// not a valid number
			std::cout << "Invalid Input! Please input a numerical value.  ";
			Beep(1200, 500);
			std::cin.clear();               // clear error flags
			std::cin.ignore(100000, '\n');  // discard invalid input up to newline/enter
		}
	}
	return result;
}

/* test boost multi-precision integers */
void fact() {
	using namespace boost::multiprecision;

	mpz_int v = 1;

	// Do some arithmetic:
	for (unsigned i = 1; i <= 1000; ++i)
		v *= i;

	std::cout << "1000!= " << v << std::endl; // prints 1000!
}
/* try to calculate previous term */
bool PrevValue(mpz_int x, mpz_int y, long long p, long long q, long long k, 
	long long r, long long s, long long l,mpz_int* xprev, mpz_int*yprev) {
/*    - if k and l are zero
      Y(n-1) = (X(n)/P - Y(n)/R) / (Q/P - S/R)
             = (x(n)*R - Y(n)*P) / (Q*R - S*P)  - multiply num & den by PR
and   X(n-1) = (X(n)/Q - Y(n)/S) / (P/Q - R/S)
             = (X(n)*S - Y(n)*Q) / (P*S - R*Q)  - multiply num & den by QS 
whith K, l non-zero we have:
	Y(n-1) = ((x(n)-K)*R - (Y(n)-L)*P) / (Q*R - S*P)
	X(n-1) = ((X(n)-K)*S - (Y(n)-L)*Q) / (S*P - Q*R) 
	
	it seems the denominator is always 1 or -1 (why?)
	*/

	mpz_int t1, t2;
	mpz_int num, denom;

	t1 = ((x-k)*r - (y-l)* p);
	t2=  (q*r - s*p);

	if (t2 == 0)
		return false;
	if (t2 != 1 && t2 != -1)
		std::cout << "q*r-s*p = " << t2 << "\n";
	if (t1%t2 != 0) {
		//std::cout << "yprev not integer\n";
		return false;
	}
	*yprev = t1 / t2;
	//std::cout << "yprev =" << t1 / t2 << "\n";
	t1 = (x-k)*s - (y-l)*q;
	t2 = p*s - r*q;

	if (t2 == 0)
		return false;
	if (t1%t2 != 0) {
		//std::cout << "xprev not integer\n";
		return false;
	}
	*xprev = t1 / t2;
	//std::cout << "xprev =" << t1 / t2 << "\n";
	return true;
}

int main(int argc, char* argv[]) {
	//fact();
	long long int a, b, c, d, e, f, listcount = 0;
	mpz_int value, x, y, x1, y1, xlist[100], ylist[100];
	long long int p, q, r, s, k, l;
	char yn = '\0';
	bool test = false;
	if (argc > 1) {
		IpFromFile = true;
		std::cout << "assume input is from a file\n";
	}
	else std::cout << "assume input is from console window\n";

	std::cout << "Check Diophantine equations of the form: Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0" << "\n";
	a = getnumber("Enter value for A ");
	b = getnumber("Enter value for B ");
	c = getnumber("Enter value for C ");
	d = getnumber("Enter value for D ");
	e = getnumber("Enter value for E ");
	f = getnumber("Enter value for F ");

	/* get x,y values and check them. */
	while (true) {
		x = getbignumber("Enter value for x ");
		y = getbignumber("Enter value for y ");
		value = a*x*x + b*x*y + c*y*y + d *x + e*y + f;

		if (value == 0) {
			std::cout << "x=" << x << "  y=" << y << " are valid\n";
			xlist[listcount] = x;
			ylist[listcount] = y;
			listcount++;
		}
		else {
			std::cout << "x=" << x << "  y=" << y << " are ** NOT VALID **\n";
		}
		yn = '\0';
		while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
			if (!IpFromFile)
				std::cout << "Enter new values for x,y) (Y/N): ";
			std::cin >> yn;
			if (toupper(yn) != 'Y' && toupper(yn) != 'N') {
				Beep(1200, 500);
				std::cin.clear();               // clear error flags
				std::cin.ignore(100000, '\n');  // discard invalid input up to newline/enter
			}
			test = (toupper(yn) == 'Y');
		}
		if (!test) 
			break;  // no more x,y pairs to enter
	}

	if (listcount == 0)
		return EXIT_FAILURE;

	for (int i = 0; i < listcount; i++) {
		std::cout << i << "  x=" << xlist[i] << "  \ty=" << ylist[i] << "\n";
	}
	
	/* get P, Q, K, R, S, L values and check them */
	while (true) {
		yn = '\0';
		while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
			if (!IpFromFile)
				std::cout << "Enter P, Q, K, R, S, L (enter 0 for unused values)? (Y/N): ";
			std::cin >> yn;
			if (toupper(yn) != 'Y' && toupper(yn) != 'N') {
				Beep(1200, 500);
				std::cin.clear();               // clear error flags
				std::cin.ignore(100000, '\n');  // discard invalid input up to newline/enter
			}
			test = (toupper(yn) == 'Y');
		}

		if (test) {
			p = getnumber("Enter value for P ");
			q = getnumber("Enter value for Q ");
			k = getnumber("Enter value for K ");
			r = getnumber("Enter value for R ");
			s = getnumber("Enter value for S ");
			l = getnumber("Enter value for L ");
			if (IpFromFile)
				std::cout << "P,Q,K,R,S,L are: " << p << " " << q << " "
				<< k << " " << r << " " << s << " " << l << "\n";

			for (int i = 0; i < listcount; i++) {
				x1 = p*xlist[i] + q*ylist[i] + k;
				y1 = r*xlist[i] + s*ylist[i] + l;
				value = a*x1*x1 + b*x1*y1 + c*y1*y1 + d*x1 + e*y1 + f;
				if (value == 0) {
					if (!IpFromFile)
						std::cout << " P, Q, K, R, S, L values are valid";
					std::cout << " x=" << xlist[i] << " y=" << ylist[i] << "\n";
					std::cout << " x1=" << x1 << " y1=" << y1 << "\n";
				}
				else {
					std::cout << "** value not zero; it is " << value;
					std::cout << " x=" << xlist[i] << " x1=" << x1;
					std::cout << " y=" << ylist[i] << " y1=" << y1 << "\n";
				}
				if (PrevValue(xlist[i], ylist[i], p, q, k, r, s, l, &x1, &y1))
				{
					value = a*x1*x1 + b*x1*y1 + c*y1*y1 + d*x1 + e*y1 + f;
					if (value == 0) {
						std::cout  << " xprev=" << x1<< " yprev=" << y1 << "\n";
					}
				}
			}
		}
		else break;
	}

	system("PAUSE");   // press any key to continue
    return 0;
}

