﻿#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>  

#include "stdafx.h"
#include <windows.h> 
// Diophantine Quadratic Equation Solver
// ax^2 + bxy + cy^2 + dx + ey + f = 0 (unknowns x,y integer numbers)
//
// Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
// (Last updated December 15th, 2003)
//
// No part of this code can be used for commercial purposes without
// the written consent from the author. Otherwise it can be used freely.
// 
// converted fron Java to C++ by R Atherton, Octobler 2017
// because of problems with java not being supported by web browsers.

//using namespace std;

/* following used to detect memory leakage */
_CrtMemState memState1;          // initial memory state
_CrtMemState memState2;          // current memory state
_CrtMemState memStatDiff;        // difference between initial and current memory state

std::string info;
std::string txt;
//const std::string sq =  "²";        // would be nice to use, but causes problems 

const std::string sq = "^2";
const std::string msg = "There are no solutions !!!\n";
std::string number;
bool also = false, ExchXY = false, teach = false;
bool allSolsFound;
const std::string divgcd = "Dividing the equation by the GCD we obtain: \n";

long long g_A, g_B, g_C, g_D, g_E, g_F;
long long g_Xi, g_Xl, g_Yi, g_Yl;
long long g_CY1, g_CY0;
long long g_A2, g_B2;
long long g_Disc;

mpz_t Dp_NUM, Dp_DEN;

unsigned long long SqrtDisc;

/* large numbers; shown by Bi_ prefix */
mpz_t Bi_H1, Bi_H2, Bi_K1, Bi_K2, Bi_L1, Bi_L2;

int NbrSols = 0, NbrCo;
int digitsInGroup = 6;     // used in ShowLargeNumber

/* typedef for vector-type arrays */
typedef  std::vector <void *> ArrayLongs;

ArrayLongs sortedSolsXv ;
ArrayLongs sortedSolsYv ;

void listLargeSolutions();

/* print a string variable. This is a relic from the original java program */
void w(std::string texto) {
	//	info = info.append(texto);
	std::cout << texto;
}

/* calculate integer square root */
unsigned long long llSqrt(unsigned long long num) {
	unsigned long long num1 = 0;
	unsigned long long num2 = (long long)65536 * (long long)32768;  /* 2^31 */
	while (num2 != 0) {
		if ((num1 + num2)*(num1 + num2) <= num) {
			num1 += num2;
		}
		num2 /= 2;
	}
	return num1;
}

/* Out = Nbr */
void LongToDoublePrecLong(long long Nbr, mpz_t Dp_Out) {
	mpz_set_si(Dp_Out, Nbr);
}

/* calculate Greatest Common Divisor */
long long gcd(long long M, long long N) {
	long long P = M;
	long long Q = N;
	while (P != 0) {
		long long R = Q%P;
		Q = P;
		P = R;
	}
	return abs(Q);
}

/* calculate Greatest Common Divisor for Bigints*/
void gcd(const mpz_t Dp_Nbr1, const mpz_t Dp_Nbr2, mpz_t Dp_Gcd) {
	mpz_gcd(Dp_Gcd, Dp_Nbr1, Dp_Nbr2);
}

/* convert double to a string */
//std::string doubToStr(double num) {
//	char numbr[21];
//	sprintf_s(numbr, sizeof(numbr), "%g", num);
//	number = numbr;   // copy to C++ string 
//	return number;
//}

void showAlso() {
	if (also) {
		printf("and also: ");
	}
	else {
		also = true;
	}
}

/* print values of X and Y */
void ShowXY(long long X, long long Y) {
	showAlso();
	if (!ExchXY) {
		//w("x = " + numToStr(X));
		printf("x = %lld  y = %lld \n", X, Y);
		//w("\ny = " + numToStr(Y) + "\n");
	}
	else {
		//w("x = " + numToStr(Y));
		printf("x = %lld  y = %lld \n", Y, X);
		//w("\ny = " + numToStr(X) + "\n");
	}
}

/* calculate and print values of X and Y. */
void ShowX1Y1(long long X1, long long Y1, long long A, long long B, long long D,
	long long D1, long long E1) {
	long long X, Y;

	if ((Y1 - E1) % D1 == 0) {
		Y = (Y1 - E1) / D1;
		if ((X1 - B*Y - D) % (2 * A) == 0) {
			X = (X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);          // display solution
		}
		if (X1 != 0 && (-X1 - B*Y - D) % (2 * A) == 0) {
			X = (-X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
	}

	if (Y1 != 0 && (-Y1 - E1) % D1 == 0) {
		Y = (-Y1 - E1) / D1;
		if ((X1 - B*Y - D) % (2 * A) == 0) {
			X = (X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
		if (X1 != 0 && (-X1 - B*Y - D) % (2 * A) == 0) {
			X = (-X1 - B*Y - D) / (2 * A);
			ShowXY(X, Y);           // display solution
		}
	}
}

/* output value of num. Uses global variable txt.
iff t is odd, use '+' for +ve numbers.
if t =2 result is left in txt and not output.
if num is zero do nothing */
int Show(long long num, std::string str, int t) {
	if (t == 2) {
		txt = "";    // clear contents of txt
	}
	if (num != 0) {
		std::string str1 = "";
		if ((t & 1) != 0 && num>0) {
			str1 += " +";
		}
		if (num<0) {
			str1 += " -";
		}
		if (abs(num) != 1) {
			str1 += " " + numToStr(abs(num));
		}
		str1 += str;
		if ((t & 2) != 0) {
			txt += str1;
		}
		else {
			//printf("%s", str1.c_str());
			std::cout << str1;
		}
		return (t | 1);  // return type bit for next time
	}
	return t;
}

void Show1(long long num, int t) {
	char u = Show(num, "", t);
	if ((u & 1) == 0 || abs(num) == 1) {
		if ((u & 2) != 0) {
			txt += numToStr(abs(num));
		}
		else {
			//w("" + numToStr(abs(num)));
			std::cout << abs(num);
		}
	}
}

/* print values of D, E and F */
void ShowLin(long long D, long long E, long long F, std::string x, std::string y) {
	int t = Show(D, x, 0);
	t = Show(E, y, t);
	Show1(F, t);
}

/* print "ax^2 + Bxy + Cy^2 +Dx + Ey + F" */
void ShowEq(long long A, long long B, long long C, long long D, long long E, long long F,
	std::string x, std::string y) {
	int t = Show(A, x + sq, 0);
	t = Show(B, x + y, t);
	t = Show(C, y + sq, t);
	t = Show(D, x, t);
	t = Show(E, y, t);
	Show1(F, t);
}

/* print " so there are no integer solutions" and set also to 'true' */
void NoSol() {
	printf(" so there are no integer solutions.\n");
	also = true;
}

void NoGcd(long long F) {
	if (teach) {
		std::cout << "This gcd is not a divisor of the constant term (" << F << ")," << "\n";
	}
	NoSol();
}

/* Nbr = -Nbr*/
void ChangeSign(mpz_t Bi_Nbr) {
	mpz_neg(Bi_Nbr, Bi_Nbr);
}

/* for some reason these functions don´t exist already int he GMP/MPIR library*/
void mpz_add_si(mpz_t rop, const mpz_t op1, mpir_si op2) {
	mpir_ui ui;
	if (op2 >= 0) {
		mpz_add_ui(rop, op1, op2);
	}
	else {
		// careful: overflow when negating INT_MIN
		ui = abs(op2);  // change sign of op2
		mpz_sub_ui(rop, op1, ui);
	}
}
void mpz_sub_si(mpz_t rop, const mpz_t op1, mpir_si op2) {
	mpir_ui ui;
	if (op2 >= 0) {
		mpz_sub_ui(rop, op1, op2);
	}
	else {
		// careful: overflow when negating INT_MIN
		ui = abs(op2);  // change sign of op2
		mpz_add_ui(rop, op1, ui);
	}
}


long long floordiv(long long num, long long den) {
	if ((num<0 && den>0 || num>0 && den<0) && num%den != 0) {
		return num / den - 1;
	}
	return num / den;
}

// commented out because it's not actually used
//long long ceildiv(long long num, long long den) {
//	if ((num>0 && den>0 || num<0 && den<0) && num%den != 0) {
//		return num / den + 1;
//	}
//	return num / den;
//}

/* q = n/d. return quotient as a bigint. Uses floor division.
i.e. rounds q down towards -infinity. The only difference from
< DivLargeNumber> is that this function does not return a value but
< DivLargeNumber> returns the remainder as a long long. */
void DivideDoublePrecLong(const mpz_t n, const mpz_t d, mpz_t q) {
	mpz_fdiv_q(q, n, d);          // note use of floor division
}

/* q = n/d    return value is remainder. Uses floor division 
i.e. rounds q down towards -infinity, and r will have the same sign as d.
the 3 different types of division only give different results when
the remainder is non-zero and n or d (or both) are negative. 
NB. The choice of floor or truncation type division is sometimes
critical, and the correct choice is not obvious. */
long long DivLargeNumber(const mpz_t n, long long d, mpz_t q) {
	long long remainder;
	mpz_t mpr, mpd, nsave;
	
	if (d == 0) {
		fprintf(stderr, "** divide by zero error\n");
		abort();
	}

	/* because no mpz_fdiv_q_si exists in the library we need to convert the divisor to a bigint */
	mpz_init_set_si(mpd, d);       // mpd = d
	mpz_init_set(nsave, n);        // nsave = n
	mpz_init(mpr);
	mpz_fdiv_qr(q, mpr, n, mpd);  // q = n/d, mpr=remainder
	remainder = mpz_get_si(mpr);  // magnitude of remainder is < d, so it can't overflow
	if (remainder != 0) {
		//gmp_printf("**temp DivLargeNumber %Zd/%lld = %Zd  rem=%lld \n", nsave, d, q, remainder);
	}
	mpz_clears(mpr, mpd, nsave, NULL);
	return remainder;
}

/* q = n/d    return value is remainder. Uses 'truncate' division 
i.e. rounds q towards zero, and remainder will have the same sign as n.
the 3 different types of division only give different results when
the remainder is non-zero and n or d (or both) are negative
NB. The choice of floor or truncation type division is sometimes
critical, and the correct choice is not obvious. */
long long tDivLargeNumber(const mpz_t n, long long d, mpz_t q) {
	long long remainder;
	mpz_t mpr, mpd, nsave;

	if (d == 0) {
		fprintf(stderr, "** divide by zero error\n");
		abort();
	}

	/* because no mpz_tdiv_q_si exists in the library we need to convert the divisor to a bigint */
	mpz_init_set_si(mpd, d);       // mpd = d
	mpz_init_set(nsave, n);        // nsave = n
	mpz_init(mpr);
	mpz_tdiv_qr(q, mpr, n, mpd);  // q = n/d, mpr=remainder
	remainder = mpz_get_si(mpr);  // magnitude of remainder is < d, so it can't overflow
	if (remainder != 0) {
		//gmp_printf("**temp tDivLargeNumber %Zd/%lld = %Zd  rem=%lld \n", nsave, d, q, remainder);
	}
	mpz_clears(mpr, mpd, nsave, NULL);
	return remainder;
}

/* return quotient (= n/d) as a normal integer. uses floor division 
i.e. rounds q down towards -infinity, and r will have the same sign as d.*/
long long DivDoublePrec(const mpz_t n, const mpz_t d) {
	mpz_t q, r;
	long long llquot;

	mpz_inits(q, r, NULL);
	mpz_fdiv_qr(q, r, n, d);             // note use of floor division
	assert(mpz_fits_slong_p(q) != 0);   // is quotient too big for normal integer?
	llquot = mpz_get_si(q);          // convert to normal integer
	/*if (!IsZero(r)) {
		gmp_printf("**temp DivDoublePrec: %Zd = %Zd/%Zd  rem=%Zd \n", q, n, d, r);
	}*/
	mpz_clears(q, r, NULL);		           // avoid memory leakage
	return llquot;
}

/*  Dest = Cprev*Bi_Prev + CAct*Bi_Act */
void MultAddLargeNumbers(long long CPrev, const mpz_t Bi_Prev,
	long long CAct, const mpz_t Bi_Act, mpz_t Bi_Dest) {
	mpz_t temp1, temp2;
	mpz_inits(temp1, temp2, NULL);    
	mpz_mul_si(temp1, Bi_Prev, CPrev);  // temp1 = Cprev*Bi_Prev
	mpz_mul_si(temp2, Bi_Act, CAct);    // temp2 =  CAct*Bi_Act
	mpz_add(Bi_Dest, temp1, temp2);     // Dest = Cprev*Bi_Prev + CAct*Bi_Act
	mpz_clears(temp1, temp2, NULL);
}

/* calculate K and L. values are returned in Bi_H1, Bi_K1, Bi_L1 and Bi_L2
It returns true if K and L are valid 
Uses input global variables g_A, g_B, g_C, g_D, g_E, Bi_H2, Bi_K2
called from ShowSols */
bool CalculateKandL() {
	/*std::cout << "**temp CalculateKandL " 
		<< " Bi_H2=" << numToString(Bi_H2) << " Bi_K2=" << numToString(Bi_K2);
	std::cout << " A=" << g_A << " B=" << g_B << " C=" << g_C << " D=" << g_D
		<< " E=" << g_E << "\n";*/

	MultAddLargeNumbers(2, Bi_H2, g_B, Bi_K2, Bi_L1);  /* 2r + Bs */
	AddLarge(Bi_L1, -2, Bi_H1);  /* Kd = 2r + Bs - 2 */
	AddLarge(Bi_H2, -1, Bi_K1);  /* r - 1 */
	MultAddLargeNumbers(-g_B, Bi_K1, -2 * g_A*g_C, Bi_K2, Bi_K1);  /* Ke */
	MultAddLargeNumbers(g_C*g_D, Bi_H1, g_E, Bi_K1, Bi_L1);  /* K(4AC - BB) */
	if (tDivLargeNumber(Bi_L1, 4 * g_A*g_C - g_B*g_B, Bi_L1) != 0) {  /* K */
		//std::cout <<"**temp CalculateKandL: K not integer; return false\n";
		return false;               /* K not integer */
	}
	MultAddLargeNumbers(g_D, Bi_K1, g_A*g_E, Bi_H1, Bi_L2);
	/*std::cout << "**temp CalculateKandL " << " Bi_L2=" << numToString(Bi_L2) 
		<< " Divisor = " <<   4 * g_A*g_C - g_B*g_B << "\n";*/
	if (tDivLargeNumber(Bi_L2, 4 * g_A*g_C - g_B*g_B, Bi_L2) != 0) {
		/*std::cout << "**temp CalculateKandL: L not integer; return false\n";*/
		return false;               /* L not integer */
	}
	MultAddLargeNumbers(1, Bi_L2, g_D, Bi_K2, Bi_L2);    /* L */
	//std::cout << "**temp CalculateKandL: Bi_L2=" << numToString(Bi_L2) << "\n";
	return true;
}

/* print large number.  */
void ShowLargeNumber(const mpz_t Bi_Nbr){
	std::string nbrOutput = "";
	char* buffer = NULL;
	size_t msglen, index=0;

	// convert to null-terminated ascii string, base 10, with sign if -ve
	buffer = mpz_get_str(NULL, 10, Bi_Nbr);
	msglen = strnlen(buffer, 50000);  // arbitrary limit of 50000
	if (buffer[0] == '-') {  // if number is -ve
		nbrOutput = "-";
		index = 1;
	}
	for (; index < msglen; index++) {
		if ((msglen - index) % digitsInGroup == 0) {
			// divide digits into groups, if there are more than 6 digits
			if ((index > 0  && buffer[0] != '-')|| (index > 1))
			nbrOutput += " ";  
		}
		nbrOutput += buffer[index];
	}
	if (buffer[0] == '-')   // if number is -ve
		msglen--;
	free(buffer);		// avoid memory leakage
	std::cout << nbrOutput;
	if (msglen > 6)
		std::cout << " (" << msglen << " digits)";
}

/* Dest = Nbr*Coef */
void MultLargeNumber(long long Coef, const mpz_t Bi_Nbr, mpz_t Bi_Dest) {
	mpz_mul_si(Bi_Dest, Bi_Nbr, Coef);
}

/* Dest = Nbr1 *Nbr2 */
void Mult2LargeNumbers(const mpz_t Bi_Nbr1, const mpz_t Bi_Nbr2, mpz_t Bi_Dest) {
	mpz_mul(Bi_Dest, Bi_Nbr1, Bi_Nbr2);
}

/* called from ShowRecursionRoot. Uses global variables:
Bi_H1, Bi_H2, Bi_K2, Bi_L1, Bi_L2, A, C 
type = hyperbolic_homog (homogenous) or hyperbolic_gen (general hyperbolic)
return 1 if K or L not integers, otherwise zero
If K and L are integers, print values for P, Q, K, R, S, L*/
int ShowSols(equation_class type) {

	//std::cout << "**temp ShowSols  Bi_H2=" << numToString(Bi_H2)
	//	<< " Bi_K2=" << numToString(Bi_K2) << " g_C=" << g_C << "\n";
	if (type == hyperbolic_gen) {
		if (!CalculateKandL()) {     /* if K or L not integers */
			return 1;                      /* bye */
		}
		/* values for K and L now in Bi_K1,  Bi_L1, Bi_L2*/
	}
	if (teach) {
		printf("m = ");
		ShowLargeNumber(Bi_H2);
		printf("\nn = ");
		ShowLargeNumber(Bi_K2);
		printf("\nUsing the formulas: P = m\n");
		printf("Q = -Cn \n");
		if (type == hyperbolic_gen) {
			printf("K =CD(P+S-2) + E(B-Bm-2ACn)4AC-B^2\n");
		}
		printf("R = An \nS = m + Bn\n");
		if (type == hyperbolic_gen) {
			printf("L =D(B-Bm-2ACn) + AE(P+S-2)4AC - B^2 + Dn\n");
		}
		printf("we obtain:\n");
	}
	printf("P = ");
	ShowLargeNumber(Bi_H2);
	printf("\nQ = ");
	MultLargeNumber(-g_C, Bi_K2, Bi_H1);  // H1 = -C*K2
	ShowLargeNumber(Bi_H1);
	if (type == hyperbolic_gen) {
		printf("\nK = ");
		ShowLargeNumber(Bi_L1);
	}
	printf("\nR = ");
	MultLargeNumber(g_A, Bi_K2, Bi_H1);   // H1 = A*K2
	ShowLargeNumber(Bi_H1);
	printf("\nS = ");
	MultAddLargeNumbers(1, Bi_H2, g_B, Bi_K2, Bi_H1);
	ShowLargeNumber(Bi_H1);
	if (type == hyperbolic_gen) {
		printf("\nL = ");
		ShowLargeNumber(Bi_L2);
	}
	putchar('\n');
	return 0;
}

/* Sum = Nbr1 +Nbr2 */
void AddLarge(const mpz_t Bi_Nbr1, const mpz_t Bi_Nbr2, mpz_t Bi_Sum) {
	mpz_add(Bi_Sum, Bi_Nbr1, Bi_Nbr2);
}

/* dest = Src + Nbr */
void AddLarge(const mpz_t Bi_Src, long long Nbr, mpz_t Bi_Dest) {
	mpz_add_si(Bi_Dest, Bi_Src, Nbr);
}

/* called from ShowRecursion
type = hyperbolic_homog (homogenous) or hyperbolic_gen (general hyperbolic) 
Uses global variables Bi_H1, Bi_H2, Bi_K2, Bi_L1, Bi_L2, A, B, C, D, E */
void ShowRecursionRoot(equation_class type) {
	char t;   // 1 if K and L not integers, otherwise 0
	t = ShowSols(type);
	if (type == hyperbolic_gen) {
		ChangeSign(Bi_H2);
		ChangeSign(Bi_K2);
		if (t == 0) putchar('\n');    /* separate 2 sets of values for  P, Q, K, R, S, L 
						  If there is only 1 set we just get an unneeded blank line.*/
		//if (t == 1 && ShowSols(type) == 1) {
		/* changed this so that ShowSols is called again even if t is 1.
		ShowSols has an important side-effect; it prints values for P, Q, K, R, S, L
		if it can find them. Sometimes each call produces a different but valid
		set of values */
		if (ShowSols(type) == 1 && t ==1) {
			/* if we get to here we haven't yet got any values for P, Q, K, R, S, L.
			The 3rd way below tends to produce very large numbers. */
			if (teach) {
				printf("m = ");
				ShowLargeNumber(Bi_H2);
				printf("\nn = ");
				ShowLargeNumber(Bi_K2);
				std::cout << "\nUsing the formulas: P = m" << sq << " - ACn" << sq <<
					" \nQ = -Cn(2m + Bn) \n K = -n(Em + CDn)\n" <<
					"R = An(2m + Bn)\nS = m" << sq << " + 2Bmn + (B" + sq << " - AC)n" << sq <<
					"\nL = n(Dm + (BD-AE)n) \nwe obtain:\n";
			}
			Mult2LargeNumbers(Bi_H2, Bi_H2, Bi_H1);   /* m^2 */
			Mult2LargeNumbers(Bi_H2, Bi_K2, Bi_K1);   /* mn */
			Mult2LargeNumbers(Bi_K2, Bi_K2, Bi_L1);   /* n^2 */
			MultAddLargeNumbers(1, Bi_H1, -g_A*g_C, Bi_L1, Bi_L2);
			printf("P = ");
			ShowLargeNumber(Bi_L2);
			MultAddLargeNumbers(-2 * g_C, Bi_K1, -g_B*g_C, Bi_L1, Bi_L2);
			printf("\nQ = ");
			ShowLargeNumber(Bi_L2);
			MultAddLargeNumbers(-g_E, Bi_K1, -g_C*g_D, Bi_L1, Bi_L2);
			printf("\nK = ");
			ShowLargeNumber(Bi_L2);
			MultAddLargeNumbers(2 * g_A, Bi_K1, g_A*g_B, Bi_L1, Bi_L2);
			printf("\nR = ");
			ShowLargeNumber(Bi_L2);
			MultAddLargeNumbers(g_B*g_B - g_A*g_C, Bi_L1, 2 * g_B, Bi_K1, Bi_L2);
			AddLarge(Bi_L2, Bi_H1, Bi_L2);
			printf("\nS = ");
			ShowLargeNumber(Bi_L2);
			MultAddLargeNumbers(g_D, Bi_K1, g_B*g_D - g_A*g_E, Bi_L1, Bi_L2);
			printf("\nL = ");
			ShowLargeNumber(Bi_L2);
			putchar('\n');
		}
	}
}

/* called from solveEquation. Returns true if equation has no solutions*/
bool CheckMod(long long R, long long S, long long X2, long long X1, long long X0) {
	long long Y2 = gcd(R, S);
	long Y1 = 2;
	int indH = 1;
	long long D1 = abs(Y2);
	long long factors[64];
	if (teach) {
		if (Y2 != 1) {
			printf("We try to check the equation modulo the prime divisors of ");
			if (R != 0 && S != 0) {
				//w("gcd(" + numToStr(abs(R)) + "," + numToStr(abs(S)) + ") = " + numToStr(Y2) + ".");
				printf("gcd(%lld, %lld) = %lld.", R, abs(S), Y2);
			}
			else {
				//w(numToStr(Y2) + ".\n");
				printf("%lld. \n", Y2);
			}
		}
	}
	/* store factors */
	while (D1 >= Y1*Y1) {
		int T = 0;
		while (D1%Y1 == 0) {
			T++;
			if (T == 1) {
				factors[indH++] = Y1;  // only save factor once
			}
			D1 /= Y1;      // remove factor
		}
		Y1++;
		if (Y1>3) {
			Y1++;  // skip even numbers other than 2
		}
	}
	if (D1>1) {
		factors[indH++] = D1;   
	}

	/* now have list of factors. indH is the number of factors*/
	for (int T = 1; T<indH; T++) {
		if (factors[T]>1) {
			long long Z = ((X1*X1 - 4 * X0*X2) % factors[T]) + factors[T];
			long long N = (factors[T] - 1) / 2;
			long long Y = 1;
			while (N != 0) {
				if (N % 2 != 0) {
					Y = (Y*Z) % factors[T];
				}
				N /= 2;
				Z = (Z*Z) % factors[T];
			}
			if (Y>1) {
				if (teach) {
					printf("There are no solutions modulo %lld,", factors[T]);
					NoSol();
				}
				return true;   // equation has no solutions
			}
		}
	}
	return false;  // equation appears to have solution(s)
}

/* reverse sign of Nbr */
void ChSign(mpz_t Dp_Nbr) {
	mpz_neg(Dp_Nbr, Dp_Nbr);
}

/* return true iff N1 == N2 */
bool AreEqual(const mpz_t Dp_N1, const mpz_t  Dp_N2) {
	int rv = mpz_cmp(Dp_N1, Dp_N2);
	return (rv == 0);
}

/* return true iff N1 == 0 */
bool IsZero(const mpz_t Dp_N1) {
	int rv = mpz_sgn(Dp_N1);
	return (rv == 0);
}

bool IsNeg(const mpz_t Dp_N1) {
	int rv = mpz_sgn(Dp_N1);
	return (rv < 0);
}

/* return true iff A == 1 */
bool IsOne(const mpz_t Dp_A) {
	int rv = mpz_cmp_si(Dp_A, 1);
	return (rv == 0);
}

/* Diff = Nbr1- Nbr2 */
void SubtDoublePrecLong(const mpz_t Dp_Nbr1, const mpz_t Dp_Nbr2, mpz_t Dp_Diff) {
	mpz_sub(Dp_Diff, Dp_Nbr1, Dp_Nbr2);
}

/* convert Nbr to string */
std::string numToStr(const mpz_t Dp_Nbr) {
	std::string nbrOutput;
	char* buffer = NULL;
	// convert to null-terminated ascii string, base 10, with sign if -ve
	buffer = mpz_get_str(NULL, 10, Dp_Nbr);   
	nbrOutput = buffer;  // copy from C-style string to STL-type string
	free(buffer);		// avoid memory leakage
	return nbrOutput;
}

/* convert integer to a string */
std::string numToStr(long long num) {
	char numbr[21];
	sprintf_s(numbr, sizeof(numbr), "%lld", num);
	number = numbr;  // copy to C++ string 
	return number;
}

/* print ax^2 +bxy +cy2  (with apppropriate signs)*/
void ShowBigEq(mpz_t Dp_A, mpz_t Dp_B, mpz_t Dp_C, std::string x, std::string y) {
	gmp_printf("%Zd%s%s ", Dp_A, x.c_str(), sq.c_str());
	if (!IsZero(Dp_B)) {
		gmp_printf("%+Zd%s%s ", Dp_B, x.c_str(), y.c_str());
	}
	if (!IsZero(Dp_C)) {
		gmp_printf("%+Zd", Dp_C);
	}
	if (y.length() == 0)
		std::cout << " ";
	else
		std::cout << " " << y << sq << " ";
}

/* convert biginteger to normal. Checks for overflow */
long long DoublePrecToLong(const mpz_t x) {
	long long rv;
	assert(mpz_fits_slong_p(x) != 0);   // is x too big for normal integer?
	rv = mpz_get_si(x); // convert to normal integer
	return rv;
}

/* prepare for solving by continued fractions.
return values in global variables Disc, SqrtDisc, Dp_NUM, Dp_DEN */
void GetRoot(mpz_t Dp_A, mpz_t Dp_B, mpz_t Dp_C) {
	long long A, B, C, M, P, T, Z;
	mpz_t Dp_P, Dp_M, Dp_Z, Dp_G, Dp_K, Dp_zz;
	bool NUMis0;

	mpz_inits(Dp_P, Dp_M, Dp_Z, Dp_G, Dp_K, Dp_zz, NULL);

	A = DoublePrecToLong(Dp_A);            // assume that Dp_A will fit int 64 bits
	B = DoublePrecToLong(Dp_B);
	C = DoublePrecToLong(Dp_C);
	g_Disc = B*B - 4 * A*C;                        // Discriminant = B^2 -4AC

	mpz_set(Dp_NUM, Dp_B);
	ChSign(Dp_NUM);                // NUM = -Dp_B
	NUMis0 = IsZero(Dp_NUM);					/* check whether NUM == 0 */
	AddLarge(Dp_A, Dp_A, Dp_DEN);       // DEN = Dp_A*2

	if (teach) {
		printf("We have to find the continued fraction expansion of the roots of \n");
		ShowBigEq(Dp_A, Dp_B, Dp_C, "t", "");
		printf("= 0, that is, ");
		printf("(Sqrt(%lld) ", g_Disc);
		if (!NUMis0) {  // print NUM if it's not zero
						//w("(Sqrt(" + numToStr(Disc) + ") ");
						//printf("(Sqrt(%lld) ", Disc);
			gmp_printf("%+Zd) ", Dp_NUM);  // print num with sign + or -
		}

		if (IsNeg(Dp_DEN )) {
			std::cout << " / (" << numToStr(Dp_DEN) << ")";  // print /(-den) if negative
		}
		else {
			std::cout << " / " << numToStr(Dp_DEN);          // print /den if positive
		}
		putchar('\n');
	}

	gcd(Dp_NUM, Dp_DEN, Dp_G);               // Dp_G = gcd(NUM,DEN)
	Mult2LargeNumbers(Dp_G, Dp_G, Dp_Z);        // Dp_Z = Dp_G^2 = gcd(NUM,DEN)^2
	LongToDoublePrecLong(g_Disc, Dp_G);		     // Dp_G = Disc = B^2-4AC
	gcd(Dp_Z, Dp_G, Dp_G);	                    // Dp_G = gcd(NUM^2, DEN^2, Disc)
	A = DoublePrecToLong(Dp_G);                 // convert gcd. assume that gcd fit into 64 bits!! 
	B = 1;
	T = 3;

	// remove any repeated factors from A (=Dp_G), put sqrt(product of repeated factors) in B
	while (A % 4 == 0) {
		A /= 4;   // remove factor 2^2
		B *= 2;
	}
	while (A >= T*T) {
		while (A % (T*T) == 0) {
			A /= T*T;   // remove factor T^2
			B *= T;
		}
		T += 2;
	}
	/* B is the product of the factors removed from A (which is the GCD of NUM^2, DEN^2 , Disc)*/
	g_Disc /= B*B;
	DivLargeNumber(Dp_NUM, B, Dp_NUM);   // NUM /= B
	DivLargeNumber(Dp_DEN, B, Dp_DEN);   // DEN /= B

	if (teach && B != 1) {
		bool DENis1 = IsOne(Dp_DEN);  /* check whether DEN == 1*/
		printf("Simplifying, ");
		if (!DENis1 && !NUMis0) {
			printf("(");
		}
		//w("Sqrt(" + numToStr(Disc) + ")");
		printf("Sqrt(%lld)", g_Disc);
		if (!NUMis0) {
			gmp_printf("%+Zd", Dp_NUM);   // print NUM with sign + or - 
		}
		if (!DENis1) {
			if (NUMis0) {
				printf(" / ");
			}
			else {
				printf(") / ");
			}
		}
		if (IsNeg(Dp_DEN)) {
			std::cout << "(" << numToStr(Dp_DEN) << ")";  // print (-DEN)
		}
		else {
			if (!DENis1) {
				std::cout << numToStr(Dp_DEN);            // print DEN
			}
		}
		putchar('\n');
	}
	g_Disc *= B*B;
	mpz_set(Dp_NUM, Dp_B);                     // NUM = B
	ChSign(Dp_NUM);                // NUM = -B
	AddLarge(Dp_A, Dp_A, Dp_DEN);       // DEN = A*2
	SqrtDisc = llSqrt(g_Disc);
	/* temporary */
	//std::cout << "**temp** getroot: DP_A=" << numToString(Dp_A);
	//std::cout << " DP_B=" << numToString(Dp_B);
	//std::cout << " DP_C=" << numToString(Dp_C);
	//std::cout << " DP_NUM=" << numToString(Dp_NUM);
	//std::cout << " DP_DEN=" << numToString(Dp_NUM);
	//std::cout << " SqrtDisc =" << SqrtDisc << "  Disc =" << g_Disc << "\n";
	/* end temporary */

	if (teach) {
		printf("The continued fraction expansion is: \n");
		mpz_set(Dp_P, Dp_DEN);   // P = DEN

		/* if DEN >= 0 then K = SqrtDisc else K = SqrtDisc+1 */
		LongToDoublePrecLong(SqrtDisc + (IsNeg(Dp_DEN) ? 1 : 0), Dp_K);
		AddLarge(Dp_K, Dp_NUM, Dp_K);  // K = K+NUM
		Z = DivDoublePrec(Dp_K, Dp_DEN);      // Z = K/DEN
		LongToDoublePrecLong(Z, Dp_M);            // M = K/DEN
		Mult2LargeNumbers(Dp_M, Dp_DEN, Dp_K);   // K = M+DEN
		SubtDoublePrecLong(Dp_K, Dp_NUM, Dp_M);   // M = K-NUM
		 //w("" + numToStr(Z));
		printf("%lld", Z);
		std::string sep = "+ //";
		int cont = -1;
		B = P = C = M = -1;

		while (cont<0 || B != P || C != M) {
			std::cout << sep;
			sep = ", ";
			mpz_add_si(Dp_zz, Dp_M, SqrtDisc);
			if (cont<0 &&
				/*  is P > 0 and P < SqrtDisc+M ?*/
				(mpz_cmp(Dp_P, Dp_zz) < 0) &&
				(mpz_sgn(Dp_P) > 0) &&
				

				/*  is M > 0 and M <= SqrtDisc ?*/
				(mpz_cmp_si(Dp_M, SqrtDisc) <= 0) &&
				(mpz_sgn(Dp_M) > 0) )
			{

				putchar('\n');
				B = P = DoublePrecToLong(Dp_P);
				C = M = DoublePrecToLong(Dp_M);
				cont = 0;
			}
			if (cont >= 0) {
				P = (g_Disc - M*M) / P;  /* both numerator and denominator are positive */
				Z = (SqrtDisc + M) / P;
				M = Z*P - M;
				cont++;
			}
			else {
				LongToDoublePrecLong(g_Disc, Dp_Z);      // Z = Disc
				Mult2LargeNumbers(Dp_M, Dp_M, Dp_G);  // G = M^2
				SubtDoublePrecLong(Dp_Z, Dp_G, Dp_G);  // G = Disc-M^2
				DivideDoublePrecLong(Dp_G, Dp_P, Dp_K);  // K=G/P
				mpz_set(Dp_P, Dp_K);                   // P=K
				LongToDoublePrecLong(SqrtDisc + (IsNeg(Dp_P) ? 1 : 0), Dp_Z);
				AddLarge(Dp_Z, Dp_M, Dp_K);
				Z = DivDoublePrec(Dp_K, Dp_P);
				LongToDoublePrecLong(Z, Dp_G);
				Mult2LargeNumbers(Dp_G, Dp_P, Dp_Z);
				SubtDoublePrecLong(Dp_Z, Dp_M, Dp_M);
			}
			//w("" + numToStr(Z));
			printf("%lld", Z);
		}
		printf("//\nwhere the periodic part is on 2nd line");
		if (cont>1) {
			//w(" (the period has " + numToStr(cont) + " coefficients)");
			printf(" (the period has %d coefficients)", cont);
		}
		putchar('\n');
	}
	mpz_clears(Dp_P, Dp_M, Dp_Z, Dp_G, Dp_K, Dp_zz, NULL);
	assert(_CrtCheckMemory());
}

/* return gcd (P,Q,R), unless gcd is not a factor of S.
In this case return 0 */
long long DivideGcd(long long P, long long Q, long long R, long long S, 
	const std::string x1, const std::string y) {
	int t;
	long long N = gcd(P, gcd(Q, R));
	if (N != 1) {
		if (teach) {
			printf("We must get the gcd of all terms except the constant: \n");
			//w("gcd(" + numToStr(P));
			printf("gcg(%lld", P);
			if (Q != 0) {
				//w(", " + numToStr(Q));
				printf(", %lld", Q);
			}
			if (R != 0) {
				//w(", " + numToStr(R));
				printf(", %lld", R);
			}
			//w(") = " + numToStr(N) + "\n");
			printf(") = %lld \n", N);
		}
		if (S%N != 0) {
			NoGcd(S);
			return 0;
		}
		if (teach) {
			std::cout << divgcd;    // show new values after division by gcd
			P /= N;
			Q /= N;
			R /= N;
			S /= N;
			ShowEq(P, 0, 0, Q, R, S, x1, y);
			printf(" = 0 \n");
		}
	}
	if (teach) {
		if (abs(R) != 1) {
			printf("This means that ");
			ShowEq(P, 0, 0, Q, 0, S, x1, y);
			//w(" should be a multiple of " + numToStr(abs(R)) + "\n");
			printf(" should be a multiple of %lld \n", abs(R));
			std::cout << "To determine this, we should try all values of " << x1 << " from 0 to " <<
				(abs(R) - 1) << " to check if the condition holds.\n";
			t = 0;
			for (long long u = 0; u<abs(R); u++) {
				if ((P*u*u + Q*u + S) % R == 0) {
					if (t != 0) {
						std::cout << ", " << u;
					}
					else {
						std::cout << "The values of " << x1 << " (mod " <<
							abs(R) << ") are: " << u;
						t = 1;
					}
				}
			}
			if (t == 0) {
				std::cout << "The modular equation is not satisfied by any " << x1;
				NoSol();
				return 0;
			}
			putchar('\n');
		}
		else {
			t = 1;
		}
	}
	return N;
}

/* compare 2 Bi numbers*/
int Compare(const mpz_t Bi_array, const mpz_t Bi_K1) {
	return mpz_cmp(Bi_array, Bi_K1);
}

int Comparev(const void * Bi_array, const mpz_t Bi_K1) {
	mpz_t temp;
	memcpy(temp, Bi_array, sizeof(mpz_t));
	return mpz_cmp(temp, Bi_K1);
}

/* Insert new number into sorted solutions vector                             */
/* performing a binary search                                                 */
/* called from ShowLargeXY                                                    */
/* Note: an mpz_t is in fact a structure which is not a type that is supported 
by the std::vector STL so the mpz_t structures and the values they contain are
are copied from the stack to  blocks in the heap and the vectors contain pointers 
to the copies on the heap. */
void InsertNewSolution(const mpz_t Bi_H1, mpz_t Bi_K1) {
	ptrdiff_t indexVector= 0, compare;
	size_t sizeVector = 0, increment;
	/* temporary  */
	//gmp_printf("\n**temp InsertNewSolution X = %Zd  Y = %Zd\n", Bi_H1, Bi_K1);  
	/* end temporary */

	sizeVector = sortedSolsYv.size();
	if (sizeVector > 0) {  // is there already something in the list?
		increment = 1;
		while (increment * 2 <= sizeVector) {
			increment *= 2;
		}
		while (increment > 0) {  /* Perform binary search */
			if (indexVector + increment <= sizeVector) {
				compare = Comparev(sortedSolsXv.at(indexVector + increment - 1), Bi_H1);
				if (compare == 0) {
					compare = Comparev(sortedSolsYv.at(indexVector + increment - 1), Bi_K1);
				}
				if (compare == 0) {
					/* temporary  */
					/*printf("** duplicate solution discarded: X = ");
					ShowLargeNumber(Bi_H1);
					printf("  Y = ");
					ShowLargeNumber(Bi_K1);
					putchar('\n');*/
					/* end temporary */
					return;  // if solution already in list don't add it again
				}
				if (compare < 0) {
					indexVector += increment;
				}
			}
			increment /= 2;
		}
	}

	mpz_t H1, K1;             // copy the values. The structure itself
	mpz_init_set(H1, Bi_H1);  // is copied to a stack variable. the value is copied to
	mpz_init_set(K1, Bi_K1);  // a new block from the heap
							  

	void *Hcopy, *Kcopy;
	Hcopy = malloc(sizeof(mpz_t));
	assert(Hcopy != NULL);
	Kcopy = malloc(sizeof(mpz_t));
	assert(Kcopy != NULL);
	memcpy(Hcopy, H1, sizeof(mpz_t));  // copy the mpz_t structures to a heap area
	memcpy(Kcopy, K1, sizeof(mpz_t));

	NbrSols++;
	//cout << "**temp InsertNewSolution NbrSols =" << NbrSols << "\n";
	ArrayLongs::iterator itX = sortedSolsXv.begin() + indexVector;
	ArrayLongs::iterator itY = sortedSolsYv.begin() + indexVector;
	sortedSolsXv.insert(itX, Hcopy);
	sortedSolsYv.insert(itY, Kcopy);
	//listLargeSolutions();    // temporary
	assert(_CrtCheckMemory());   // check for heap corruption
}

/* convert num to digits, in brackets if -ve*/
std::string par(long long num) {
	if (num<0) {
		return "(" + numToStr(num) + ")";
	}
	return "" + numToStr(num);
}

/* convert num to digits, but return empty string if num == 1*/
std::string par1(long long num) {
	return (num == 1) ? "" : par(num);
}

/* this function either calls InsertNewSolution to save the solution
or outputs the the solution: "x=<val>  y=<val>"                   
In some cases the -ve value of x and/or y is also stored/printed */
void ShowLargeXY(std::string x, std::string y, mpz_t Bi_X, mpz_t Bi_Y, 
	bool sol, const std::string eqX, const std::string eqY) {
	/* if sol == true use both +ve and -ve of solution */

	if ((x == "X") && (y == "Y")) {
		
		/*std::cout << "**temp ShowLargeXY H1="; 	ShowLargeNumber(Bi_H1);
		std::cout << " K1=";  ShowLargeNumber(Bi_K1);
		std::cout << "  sol=" << sol << "\n";*/

		if (!teach && !allSolsFound) {
			also = true;
			/* store solutions in sorted solutions vector and exit*/
			if (sol && IsZero(Bi_Y)) {
				InsertNewSolution(Bi_X, Bi_Y);
				ChangeSign(Bi_X);
				InsertNewSolution(Bi_X, Bi_Y);
				ChangeSign(Bi_X);
			}
			else {
				if (sol && IsNeg(Bi_Y)) {
					ChangeSign(Bi_X);
					ChangeSign(Bi_Y);
					InsertNewSolution(Bi_X, Bi_Y);
					ChangeSign(Bi_X);
					ChangeSign(Bi_Y);
				}
				else {
					InsertNewSolution(Bi_X, Bi_Y);
				}
			}
			if (!allSolsFound)
				std::cout  << NbrSols << " solutions \n";
			return;
		}
	}

	if (teach && y == "Y") {
		putchar('\n');
	}
	if (y == "Y") {
		showAlso();
	}
	std::cout << x << " = ";
	if (teach) {
		std::cout << eqX;
	}
	ShowLargeNumber(Bi_X);
	std::cout << "    " << y << " = ";
	if (teach) {
		std::cout << eqY;
	}
	ShowLargeNumber(Bi_Y);
	putchar('\n');

	if (y == "Y" && sol) {
		/* show -ve of solution as well */
		ChangeSign(Bi_X);
		ChangeSign(Bi_Y);
		std::cout << "\n  and also:  " << x << "0 = ";
		if (teach) {
			w(eqX == "" ? "" : "-" + eqX);
		}
		ShowLargeNumber(Bi_X);
		std::cout << "\n" << y << "0 = ";
		if (teach) {
			w(eqY == "" ? "" : "-" + eqY);
		}
		ShowLargeNumber(Bi_Y);
		putchar('\n');
		ChangeSign(Bi_X);
		ChangeSign(Bi_Y);
	}
	if (teach && y == "Y") {
		//w("</B>");
	}
}

/* show all solutions stored in sorted solutions vector, stored by
InsertNewSolution function.
Note that solutions are deleted after they are printed, making the vector
ready for reuse, but this function can only be called once. */

void ShowAllLargeSolutions() {
	size_t i;
	mpz_t xtemp, ytemp;    // note these are NOT initialised by mpz_set
	allSolsFound = true;
	also = false;

	for (i = 0; i < sortedSolsYv.size(); i++) {
		memcpy(xtemp, sortedSolsXv.at(i), sizeof(mpz_t));
		memcpy(ytemp, sortedSolsYv.at(i), sizeof(mpz_t));
		std::cout << "X= ";     ShowLargeNumber(xtemp);
		std::cout << "  \tY= "; ShowLargeNumber(ytemp);  // tab makes output a bit more tidy
		std::cout << "\n";

		mpz_clear(xtemp);  // avoid memory leakage
		mpz_clear(ytemp);
		free(sortedSolsXv[i]);
		free(sortedSolsYv[i]);
	}
	sortedSolsXv.clear();		// clear vectors to initial state.
	sortedSolsXv.shrink_to_fit();   // avoid spurious memory leakage report
	sortedSolsYv.clear();
	sortedSolsYv.shrink_to_fit();
	assert(_CrtCheckMemory());     // check for heap corruption
}


/* quick print of sortedSols - for testing */
void listLargeSolutions() {
	long long  i;
	long long sizeVector = sortedSolsYv.size();
	mpz_t xtemp, ytemp;
	for (i = 0; i < sizeVector; i++) {
		memcpy(xtemp, sortedSolsXv.at(i), sizeof(mpz_t));
		memcpy(ytemp, sortedSolsYv.at(i), sizeof(mpz_t));
		gmp_printf("X= %Zd, Y= %Zd \n", xtemp, ytemp);
	}
}

/* uses global vars Xi, Yi, X1, Y1
print x = xi +X1*t,  y = yi +Y1*t*/
void PrintLinear(std::string va) {
	if (va == "t") {         // actually, va always = t
		showAlso();
	}

	if (ExchXY) {
		long long T = g_Xi; g_Xi = g_Yi; g_Yi = T;  /* swap Xi and Yi */
		T = g_Xl; g_Xl = g_Yl; g_Yl = T;            /* swap X1 and Y1 */
	}

	printf("x = ");
	if (g_Xi == 0 && g_Xl == 0) {
		printf("0");
	}
	if (g_Xi != 0) {
		//w("" + numToStr(Xi));
		printf("%lld", g_Xi);
	}
	if (g_Xl<0) {
		printf(" - ");
	}
	if (g_Xl>0 && g_Xi != 0) {
		printf(" + ");
	}
	if (g_Xl != 0) {
		if (abs(g_Xl) == 1) {
			std::cout << " " << va;
		}
		else {
			std::cout << abs(g_Xl) << " " << va;
		}
	}

	printf("   y = ");
	if (g_Yi == 0 && g_Yl == 0) {
		printf("0");
	}
	if (g_Yi != 0) {
		//w("" + numToStr(Yi));
		printf("%lld", g_Yi);
	}
	if (g_Yl<0) {
		printf(" - ");
	}
	if (g_Yl>0 && g_Yi != 0) {
		printf(" + ");
	}
	if (g_Yl != 0) {
		if (abs(g_Yl) == 1) {
			std::cout << " " << va;
		}
		else {
			std::cout << abs(g_Yl) << " " << va;
		}
	}
	putchar('\n');
	return;
}

/* called from solveEquation, SolveParabolic & SolveDiscIsSq.
return 0 if solution found, 1 if no solutions exist, 2 if there are an infinite 
number of solutions. 
The equation to be solved is of the form Dx + Ey + F =0 */
int Linear(long long D, long long E, long long F) {
	long long Tx;
	int t;
	if (teach) {
		printf("This is a linear equation ");
		ShowLin(D, E, F, "x", "y");
		printf(" = 0 \n");
	}
	if (D == 0) {
		if (E == 0) {
			if (F != 0) {
				return 1;             // No solutions
			}
			else {
				printf("x, y: any integer\n");  /* 0X + 0Y +0 = 0 */
				return 2;             // Infinite number of solutions
			}
		}
		if (F%E != 0) {
			return 1;               // No solutions
		}
		else {
			g_Xi = 0;
			g_Xl = 1;
			g_Yi = -F / E;
			g_Yl = 0;
			PrintLinear("t");
			return 0;               // Solution found
		}
	}
	if (E == 0) {
		if (F%D != 0) {
			return 1;               // No solutions
		}
		else {
			g_Xi = -F / D;
			g_Xl = 0;
			g_Yi = 0;
			g_Yl = 1;
			PrintLinear("t");
			return 0;               // Solution found
		}
	}

	long long Q = gcd(D, E);
	if (Q != 1 && Q != -1) {
		if (teach) {
			std::cout << "To solve it, we first find the gcd of the linear coefficients, that is: gcd(" <<
				D << ", " << E << ") = " << Q << ".\n";
		}
		if (F%Q != 0) {
			NoGcd(F);
			return 1;               // No solutions
		}
		D = D / Q; E = E / Q; F = F / Q;
	}
	if (teach) {
		std::cout << divgcd;    // show new values after division by gcd
		ShowLin(D, E, F, "x", "y");
		printf(" = 0 \n");
		printf("Now we must apply the Generalized Euclidean algorithm: \n");
	}
	long long U1 = 1;
	long long U2 = 0;
	long long U3 = D;
	long long V1 = 0;
	long long V2 = 1;
	long long V3 = E;
	t = 1;
	while (V3 != 0) {
		if (teach) {
			std::cout << "Step " << t << ": " << par(U1) << " * " << par(D) << " + " << par(U2) << " * "
				<< par(E) << " = " << U3 << "\n";
		}
		long long q = floordiv(U3, V3);
		long long T1 = U1 - q*V1;
		long long T2 = U2 - q*V2;
		long long T3 = U3 - q*V3;
		U1 = V1; U2 = V2; U3 = V3;
		V1 = T1; V2 = T2; V3 = T3;
		t++;
	}
	g_Xi = -U1*F / U3; g_Xl = E; g_Yi = -U2*F / U3; g_Yl = -D;
	if (teach) {
		std::cout << "Step " << t << ": " << par(U1) << " * " << par(D) << " + " << par(U2) << " * " <<
			par(E) << " = " << U3 + "\n";
		std::cout << "\nMultiplying the last equation by " << par(-F / U3) << " we obtain:\n";
		std::cout << par(g_Xi) << " * " << par(D) << " + " << par(g_Yi) << " * " << par(E) << " = " << (-F) << "\n";
		std::cout << "Adding and subtracting " << par(D) << " * " << par(E) << " t' we obtain:\n";
		std::cout << "(" << g_Xi << " + " << par(E) << " t') * " << par(D) << " + " << g_Yi << " - " << par(D)
			<< " t') * " << par(E) << " = " << (-F) << "\n";
		printf("So, the solution is given by the set:\n");
		PrintLinear("t'");
	}
	V1 = D*D + E*E;
	Tx = floordiv((D*g_Yi - E*g_Xi) + V1 / 2, V1);
	if (teach) {
		std::cout << "By making the substitution t = " << Tx << " + t' we finally obtain:\n";
	}
	g_Xi += E*Tx;
	g_Yi += -D*Tx;
	if (g_Xl<0 && g_Yl<0) {
		g_Xl = -g_Xl;
		g_Yl = -g_Yl;
	}
	PrintLinear("t");
	return 0;
}

/* uses global variables A, B, C, D, E, F.
returns true if there are no solutions. */
bool Mod(long long mod) {
	for (long long x = 0; x<mod; x++) {
		long long z = g_A*x*x + g_D*x + g_F;
		long long t = g_B*x + g_E;
		for (long long y = 0; y<mod; y++) {
			if ((z + y*(t + g_C*y)) % mod == 0) {
				if (teach) {
					std::cout << "solution found using x= " << x << ", y=" << y;
					std::cout << " mod = " << mod << "\n";
				}
				return false;  // there is a solution
			}
		}
	}

	/* no solution found */
	if (teach) {
		std::cout << "No solutions found using mod " << mod << ",";
		NoSol();
	}
	return true;
}

/* called from SolveElliptical */
void ShowElipSol(long long A, long long B, long long D, long long u, std::string x, std::string x1,
	std::string y, long long w2) {
	if (teach) {
		putchar('\n');
		ShowLin(2 * A, B, D, x, y);
		std::cout << " = " << w2 << "\n";
		ShowLin(2 * A, B, D, x, par(u));
		std::cout << " = " << w2 << "\n";
		std::cout << par1(2 * A) << x << " = " << (w2 - D - B*u) << "\n";
		if ((w2 - D - B*u) % (2 * A) == 0) {
			std::cout << x << " = " << ((w2 - D - B*u) / (2 * A));
		}
		else {
			std::cout << "There is no integer solution for " << x << " \n";
		}
	}
	if ((w2 - D - B*u) % (2 * A) == 0) {
		if (teach) {
			putchar('\n');
		}
		ShowXY((w2 - D - B*u) / (2 * A), u);        // display solution
		putchar('\n');
		also = true;
	}
}

/* called from solveEquation
type = hyperbolic_homog (homogenous) or hyperbolic_gen (general hyperbolic)
uses global variables g_A, g_B, g_C, Bi_H2*/
void ShowRecursion(equation_class type) {
	const std::string t1 = " integer solution of the equation \nm" + sq + " + bmn + acn" + sq + " = ";
	mpz_t Dp_A, Dp_B, Dp_C;

	mpz_inits(Dp_A, Dp_B, Dp_C, NULL);

	/*std::cout << "**temp ShowRecursion A=" << g_A << " B=" << g_B << " C=" << g_C << " Bi_H2=";
	ShowLargeNumber(Bi_H2);
	std::cout << " type=" << type << "\n";*/

	if (type == hyperbolic_gen) {
		std::cout << "X(n+1) = P X(n) + Q Y(n) + K \n";		
		std::cout << "Y(n+1) = R X(n) + S Y(n) + L \n";
	}
	else {
		std::cout << "X(n+1) = P X(n) + Q Y(n) \n";
		std::cout << "Y(n+1) = R X(n) + S Y(n) \n";
	}

	if (teach) {
		std::cout << "In order to find the values of P, Q, R, S we have to find first an" << t1;
		ShowEq(1, g_B, g_A*g_C, 0, 0, 0, "m", "n");
		printf(" = 1. \n");
	}
	
	LongToDoublePrecLong(g_A, Dp_A);       // Dp_A = A
	LongToDoublePrecLong(g_C, Dp_B);       // Dp_B = C
	Mult2LargeNumbers(Dp_A, Dp_B, Dp_C);   // Dp_C = A*C
	LongToDoublePrecLong(1, Dp_A);         // Dp_A = 1
	LongToDoublePrecLong(g_B, Dp_B);       // Dp_B = B
	GetRoot(Dp_A, Dp_B, Dp_C);          /* return values in Disc, SqrtDisc, Dp_NUM, Dp_DEN */
	ContFrac(Dp_A, 2, 1, 0, 0, 1, g_A);

	if (teach) {
		std::cout << "An" << t1;
		ShowEq(1, g_B, g_A*g_C, 0, 0, 0, "m", "n");
		printf(" = 1 is: \n");
	}

	ShowRecursionRoot(type);

	/* big problem here!! it prints rubbish values!! */
	//if (g_B != 0) {
	//	if (teach) {
	//		std::cout << "\nAnother" << t1;
	//		ShowEq(1, g_B, g_A*g_C, 0, 0, 0, "m", "n");
	//		printf(" = 1 is: \n");
	//	}
	//	else {
	//		printf("\nas well as\n");
	//	}
	//	MultAddLargeNumbers(1, Bi_H2, g_B, Bi_K2, Bi_H2); /* r <- r + Bs */
	//	ChangeSign(Bi_K2);            /* s <- -s */
	//	ShowRecursionRoot(type);
	//}

	mpz_clears(Dp_A, Dp_B, Dp_C, NULL);
}

/* called from SolveSimpleHyperbolic */
void SolByFact(long long R, long long T, long long B, long long D, long long E) {
	if (teach) {
		std::cout << "\n" << T << " is a factor of " << R << ", so we can set:\n";
		Show1(E, Show(B, " x", 0));
		//w(" = " + temp);
		printf(" =  %lld", T);
		if (E != 0) {
			printf("\nThis means that ");
			Show(B, " x = ", 0);
			//w("" + (temp - E));
			printf("%lld", T - E);
		}
		putchar('\n');
		Show1(D, Show(B, " y", 0));
		std::cout << " = " << R << "/" << par(T) << " = " << (R / T);
		if (D != 0) {
			printf("\nThis means that ");
			Show(B, " y = ", 0);
			//w("" + numToStr(R / temp - D));
			printf("%lld", R / T - D);
		}
		putchar('\n');
	}
	if ((T - E) % B == 0 && (R / T - D) % B == 0) {
		if (teach) {
			//printf("<TABLE BORDER=1><TR><TH>");
		}
		ShowXY((T - E) / B, (R / T - D) / B);
		if (teach) {
			//printf("</TABLE>");
		}
	}
	else {
		if (teach) {
			printf("These equations are not valid in integers.\n");
		}
	}
}

/* The discriminant is a perfect square;
the equation can be represented as the product of 2 linear expressions.
uses global variables g_A, g_B, g_D, g_CY0, g_CY1 */
void SolveDiscIsSq(long long N0, std::string x, std::string y, long long SqrtD, long long g, long long h,
	long long N1, std::string y1, std::string x1) {
	long long Yc = llSqrt(g / h);
	long long Xc = llSqrt(abs(g_CY1 / h));
	long long temp, Fact1, Fact2, X1, Y1;

	if (teach) {
		printf("(");
		ShowLin(Yc, Xc, 0, y1, x1);
		printf(") (");
		ShowLin(Yc, -Xc, 0, y1, x1);
		//w(") = " + numToStr(N0 / h) + "\n");
		printf(") = %lld \n", (N0 / h));
	}

	if (N0 == 0) {
		if (teach) {
			printf("\nOne of the parentheses must be zero, so: \n");
		}

		/* If SqrtD is not zero perform loop twice, 2nd time with sign of Xc reversed
		if SqrtD is zero,  for loop is only execuded once. */
		for (temp = (SqrtD == 0 ? 1 : 0); temp<2; temp++) {
			if (teach) {
				//w("<LI>");
				ShowLin(Yc, Xc, 0, y1, x1);
				std::cout << " = 0 \n" << par(Yc) << " (";
				ShowLin(0, g_CY1, g_CY0, x, y);
				std::cout << ") + " << par1(Xc) << " (";
				ShowLin(2 * g_A, g_B, g_D, x, y);
				printf(") = 0 \n");
				ShowLin(2 * g_A*Xc, Xc*g_B + g_CY1*Yc, g_D*Xc + g_CY0*Yc, x, y);
				printf(" = 0 \n");
			}
			Linear(2 * g_A*Xc, Xc*g_B + g_CY1*Yc, g_D*Xc + g_CY0*Yc);
			Xc = -Xc;          // reverse sign of Xc
		}
		SqrtD = abs(SqrtD);  // pointless? SqrtD is always positive and also the changed value is not returned to the caller
		if (teach) {
			//w("</UL>");
		}
		return;
	}

	/* N0 is not zero, discriminant is a perfect square */
	if (teach) {
		//w("Now we have to find all factors of " + numToStr(N1) + ". \n");
		printf("Now we have to find all factors of %lld. \n", N1);
	}
	for (unsigned long long t1 = 1; t1 <= llSqrt(N1); t1++) {
		if (N1%t1 == 0) {
			Fact1 = t1;       // t1 is a factor of N1, so copy it to Fact1
			Fact2 = N0 / h / t1;
			if (teach) {
				std::cout << "Since " << (Fact1*Fact2) << " is equal to " << Fact1 <<
					" times " << Fact2 << ", we can set:\n";
				ShowLin(Yc, Xc, 0, y1, x1);
				//w(" = " + numToStr(Fact1) + "\n");
				printf(" = %lld \n", Fact1);
				ShowLin(Yc, -Xc, 0, y1, x1);
				//w(" = " + numToStr(Fact2) + "\n");
				printf(" = %lld \n", Fact2);
				if ((Fact1 - Fact2) % (2 * Xc) == 0 && (Fact1 + Fact2) % (2 * Yc) == 0) {
					X1 = (Fact1 - Fact2) / (2 * Xc);
					Y1 = (Fact1 + Fact2) / (2 * Yc);
					std::cout << x1 << " = " << X1 << "\n" << y1
						<< " = " << Y1 << "\n";
					ShowX1Y1(X1, Y1, g_A, g_B, g_D, g_CY1, g_CY0);  // display solution

				}
				else {
					//w("Solving this system we do not obtain integer values for " + x1 + " and " + y1 + ". \n");
					std::cout << "Solving this system we do not obtain integer values for "
						<< x1 << " and " << y1 << ". \n";
				}
			}
			else {
				if ((Fact1 - Fact2) % (2 * Xc) == 0 && (Fact1 + Fact2) % (2 * Yc) == 0) {
					X1 = (Fact1 - Fact2) / (2 * Xc);
					Y1 = (Fact1 + Fact2) / (2 * Yc);
					ShowX1Y1(X1, Y1, g_A, g_B, g_D, g_CY1, g_CY0);  // display solution
				}
			}
		}
	}
	if (teach) {
		putchar('\n');
	}
	return;
}

/* called from solveEquation 
uses global variables g_B, g_D, g_E and g_F*/
void SolveSimpleHyperbolic(void) {
	/* simple hyperbolic; A = C = 0; B ≠ 0 */
	long long R, S, T;

	R = g_D*g_E - g_B*g_F;
	if (teach) {
		//w("Multiplying by " + numToStr(B) + " we obtain:\n");
		printf("Multiplying by %lld we obtain; \n", g_B);
		ShowEq(0, g_B*g_B, 0, g_D*g_B, g_E*g_B, 0, "x", "y");
		//w(" = " + numToStr(-F*B) + "\n");
		printf(" = %lld \n", (-g_F*g_B));
		//w("Adding " + numToStr(D*E) + " to both sides of the equal sign:\n");
		printf("Adding %lld to both sides of the equal sign: \n", (g_D*g_E));
		ShowEq(0, g_B*g_B, 0, g_D*g_B, g_E*g_B, g_D*g_E, "x", "y");
		//w(" = " + numToStr(R) + "\n");
		printf(" = %lld \n", R);
		printf("Now the left side can be factored as follows:\n");
		printf("(");
		Show1(g_E, Show(g_B, " x", 0));
		printf(") (");
		Show1(g_D, Show(g_B, " y", 0));
		//w(") = " + numToStr(R) + "\n");
		printf(") = %lld \n", R);
		if (R != 0) {
			printf("Then ");
			Show1(g_E, Show(g_B, " x", 0));
			/*w(" must be a factor of " + numToStr(R) + ", so we must find the factors of "
			+ numToStr(R) + ":\n");*/
			printf(" must be a factor of %lld, so we must find the factors of %lld: \n", R, R);
		}
		else {
			printf("One of the parentheses must be zero, so: \n");
			Show1(g_E, Show(g_B, " x", 0));
			printf(" = 0");
			if (g_E%g_B == 0) {
				std::cout << " means that x = " << (-g_E / g_B) << " and y could be any integer.\n";
				also = true;
			}
			else {
				printf("This equation cannot be solved in integers.\n");
			}
			Show1(g_D, Show(g_B, " y", 0));
			printf(" = 0");
			if (g_D%g_B == 0) {
				//w(" means that y = " + numToStr(-D / B) + " and x could be any integer. \n");
				printf(" means that y = %lld and x could be any integer. \n", (-g_D / g_B));
				also = true;
			}
			else {
				printf("This equation cannot be solved in integers. \n");
			}
			return;
		}
	}
	if (R != 0) {
		S = llSqrt(abs(R));
		for (T = 1; T <= S; T++) {
			if (R%T == 0) {
				SolByFact(R, T, g_B, g_D, g_E);
				SolByFact(R, -T, g_B, g_D, g_E);
				if (T*T != abs(R)) {
					SolByFact(R, R / T, g_B, g_D, g_E);
					SolByFact(R, -R / T, g_B, g_D, g_E);
				}
			}
		} /* end for */
		if (teach) {
			//w("</UL>");
			putchar('\n');
		}
		return;
	}
	if (g_E%g_B == 0) {
		g_Xi = -g_E / g_B; g_Xl = 0; g_Yi = 0; g_Yl = 1;
		PrintLinear("t");
	}
	if (g_D%g_B == 0) {
		g_Xi = 0; g_Xl = 1; g_Yi = -g_D / g_B; g_Yl = 0;
		PrintLinear("t");
	}
	return;
}

/* uses global variables g_A, g_B, g_C, g_D, g_E, g_F*/
void SolveParabolic(std::string x, std::string y, std::string x1) {
	int t;
	std::string t1;
	long long r, r1 =1, r2, s, P, P1, P2, Q, Q1, Q2, R, R1, S, S1, S2, T, u;
	const long long E1 = 4 * g_A*g_E - 2 * g_B*g_D;
	const long long F1 = 4 * g_A*g_F - g_D*g_D;

	r = gcd(2 * g_A, g_B);
	s = 2 * g_A / r;
	P = r / 2;
	Q = g_D;
	R = (2 * g_A*g_E - g_B*g_D) / r;
	S = 2 * g_A*g_F / r;

	if (teach) {
		if (s != 1) {
			std::cout << "Multiplying the equation by " << par(s) << ":\n";
			ShowEq(g_A*s, g_B*s, g_C*s, g_D*s, g_E*s, g_F*s, x, y);
			printf(" = 0 \n");
		}
		if (r != 2) {
			std::cout << "Extracting the factor " << (r / 2) << " in the quadratic terms:\n";
			std::cout << par1(r / 2) << " (";
			ShowEq(s*s, 2 * s*g_B / r, 2 * s*g_C / r, 0, 0, 0, x, y);
			printf(")");
			t = Show(g_D*s, " " + x, 1);
			t = Show(g_E*s, " " + y, t);
			Show1(g_F*s, t);
			printf(" = 0 \n");
		}
		if (g_B != 0) {
			std::cout << par1(r / 2) << "(";
			ShowLin(s, g_B / r, 0, x, y);
			std::cout << ")" << sq;
			if (g_D != 0 || g_E != 0 || g_F != 0) {
				printf(" + (");
				ShowLin(g_D*s, g_E*s, g_F*s, x, y);
				printf(")");
			}
			printf(" = 0 \n");
			if (g_D != 0) {
				std::cout << "Adding and subtracting " << par1(g_B*g_D / r) << y << ":\n";
				std::cout << par1(r / 2) << "(";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << ")" << sq;
				std::cout << " + " << par1(g_D) << " (";
				ShowLin(s, g_B / r, 0, x, y);
				printf(")");
			}
			if (E1 != 0 || g_F != 0) {
				printf(" + (");
				ShowLin(0, R, g_F*s, x, y);
				printf(")");
			}
		}
		else {
			std::cout << par1(r / 2) << "(";
			ShowLin(s, g_B / r, 0, x, y);
			std::cout << ")" << sq;
			if (g_D != 0) {
				std::cout << " + " << par1(g_D) << " (";
				ShowLin(s, g_B / r, 0, x, y);
				printf(")");
			}
			if (E1 != 0 || g_F != 0) {
				printf(" + (");
				ShowLin(0, R, g_F*s, x, y);
				printf(")");
			}
		}
		printf(" = 0 \nNow we perform the substitution:\n");
		std::cout << x1 + " = ";
		ShowLin(s, g_B / r, 0, x, y);
		printf("\nThis gives:\n");
		ShowEq(r / 2, 0, 0, g_D, R, g_F*s, x1, y);
		printf(" = 0 \n");
	}

	if (E1 == 0) {
		if (teach) {
			printf("This can be solved by the standard quadratic equation formula: \n");
			std::cout << "The roots are: " << x1 << " = ";
			if (g_D != 0) {
				std::cout << "(" << par(-g_D) << " ± llSqrt(" << (-F1) + "))";
			}
			else {
				std::cout << " +/- Sqrt(" << (-F1) << ")";
			}
			std::cout << " / " << par(r) << "\n";
		}

		if (F1>0) {
			printf("This quadratic equation has no solution in reals, so it has no solution in integers. \n");
			also = true;
			return;
		}

		T = (long long)floor((-g_D + sqrt(-F1)) / r + 0.5);

		if (r*T*T / 2 + g_D*T + 2 * g_A*g_F / r == 0) {
			if (teach) {
				std::cout << "The first root is: " << x1 << " = " << T << "\n";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << " = " << T << "\n";
			}
			r1 = Linear(s, g_B / r, -T);
		}
		else {
			printf("The first root is not an integer. \n");
		}

		if (F1 == 0) {
			return;
		}

		T = (long long)floor((-g_D - sqrt(-F1)) / r + 0.5);

		if (r*T*T / 2 + g_D*T + 2 * g_A*g_F / r == 0) {
			if (teach) {
				std::cout << "The second root is: " << x1 << " = " << T << "\n";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << " = " << T << "\n";
			}
			r2 = Linear(s, g_B / r, -T);
		}
		else {
			printf("The second root is not an integer. \n");
		}

		if (r1 == 1 && r2 == 1)
			std::cout << msg;     // there are no solutions

		return;
	}

	long long N = DivideGcd(P, Q, R, S, x1, y);
	if (N == 0) {
		return;   // cannot divide S by gcd(P,Q,R)
	}

	P /= N; Q /= N; R /= N; S /= N;
	int numEq = 1;
	for (u = 0; u<abs(R); u++) {
		T = P*u*u + Q*u + S;
		if (T%R == 0) {
			int numEq2 = 3;
			P1 = g_B*P*R / r; Q1 = abs(R) + g_B*(2 * P*u + Q)*abs(R) / R / r;
			R1 = -s; S1 = g_B*T / R / r + u;
			P2 = -P*R; Q2 = (2 * P*u + Q)*abs(R) / -R; S2 = -T / R;
			if (teach) {
				t1 = ((u == 0) ? "" : numToStr(u) + " + ") + numToStr(abs(R)) + "t";
				std::cout << x1 << " = " << t1 << "   (" << numEq << ")";
				printf("\nReplacing this in the equation shown above:\n");
				std::cout << par1(-R) << y << " = ";
				ShowEq(P, 0, 0, Q, 0, S, "(" + t1 + ")", y);
				std::cout << "\n" << par1(-R) << y << " = ";
				ShowEq(P*R*R, 0, 0, (2 * P*u + Q)*abs(R), 0, T, "t", y);
				std::cout << "\n" << y << " = ";
				ShowEq(P2, 0, 0, Q2, 0, S2, "t", y);
				std::cout << "   (" << (numEq + 1) << ")\nFrom (" << numEq << "): ";
				ShowLin(s, g_B / r, 0, x, y);
				std::cout << " = " << t1 << "\nReplacing (" << (numEq + 1) << ") here:\n";
				std::cout << par1(s) << x << " = ";
				ShowEq(P1, 0, 0, Q1, 0, S1, "t", y);
				std::cout << "   (" << (numEq + 2) << ") \n";
			}
			long long N1 = DivideGcd(P1, Q1, R1, S1, "t", x);
			if (N1 == 0) {
				continue;   // cannot divide S by gcd(P,Q,R)
			}
			P1 /= N1; Q1 /= N1; R1 /= N1; S1 /= N1;
			for (long long u1 = 0; u1<abs(R1); u1++) {
				long long T1 = P1*u1*u1 + Q1*u1 + S1;
				if (T1%R1 == 0) {
					if (teach && abs(R1) != 1) {
						t1 = ((u1 == 0) ? "" : u1 + " + ") + numToStr(abs(R1)) + "u";
						std::cout << "t = " << t1 << "   (" << (numEq + numEq2) << ")";
						std::cout << "\nReplacing this in (" << (numEq + 2) << "):\n";
						std::cout << par1(-R1) << x << " = ";
						ShowEq(P1, 0, 0, Q1, 0, S1, "(" + t1 + ")", y);
						std::cout << "\n" << par1(-R1) << x << " = ";
						ShowEq(P1*R1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1), 0, T1, "u", "");
						std::cout << "\n" << x << " = ";
						ShowEq(-P1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1) / -R1, 0, -T1 / R1, "u", y);
						std::cout << "\nFrom (" << (numEq + 1) << ") and (" <<
							(numEq + numEq2) << "):\n" << y << " = ";
						ShowEq(P2, 0, 0, Q2, 0, S2, "(" + t1 + ")", y);
						std::cout << "\n" << y << " = ";
						ShowEq(P2*R1*R1, 0, 0, (2 * P2*u1 + Q2)*abs(R1), 0, P2*u1*u1 + Q2*u1 + S2, "u", "");
						putchar('\n');
						numEq2++;
					}
					if (teach) {
						//w("<TABLE BORDER=1><TR><TH>");
					}
					showAlso();
					printf("x = ");
					if (!ExchXY) {
						ShowEq(-P1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1) / -R1, 0, -T1 / R1, "u", "");
						printf("\ny = ");
						ShowEq(P2*R1*R1, 0, 0, (2 * P2*u1 + Q2)*abs(R1), 0, P2*u1*u1 + Q2*u1 + S2, "u", "");
					}
					else {
						ShowEq(P2*R1*R1, 0, 0, (2 * P2*u1 + Q2)*abs(R1), 0, P2*u1*u1 + Q2*u1 + S2, "u", "");
						printf("\ny = ");
						ShowEq(-P1*R1, 0, 0, (2 * P1*u1 + Q1)*abs(R1) / -R1, 0, -T1 / R1, "u", "");
					}
					if (teach) {
						//w("</TABLE>");
					}
					putchar('\n');
				}
			}
			if (teach && abs(R1) != 1) {
				//w("</UL>");
			}
			numEq += numEq2;
		}
	}

	if (teach) {
		//w("</UL>");
	}
	return;   /* end of specific processing for parabolic case*/
}

/* uses global variables g_A, g_B, g_C, g_D, g_E, g_F */
void SolveElliptical(std::string x, std::string x1, std::string y) {
	int b;
	long long u, w1, w2, g;

	long long NegDisc = 4 * g_A*g_C - g_B*g_B;   // get -ve of discriminant
	long long E1 = 4 * g_A*g_E - 2 * g_B*g_D;
	long long F1 = 4 * g_A*g_F - g_D*g_D;
	g = gcd(NegDisc, E1 / 2);
	g_CY1 = NegDisc / g;
	g_CY0 = E1 / 2 / g;
	long long N0 = g_CY0*g_CY0*g - g_CY1*F1;
	
	double sqrtgN0 = sqrt((double)g*(double)N0);

	double R3 = (-E1 / 2 - sqrtgN0) / NegDisc;    // get minimum for x as a real number
	double R4 = (-E1 / 2 + sqrtgN0) / NegDisc;    // get maximum for x as a real number
	long long R2 = (long long)floor(R4);          // get maximum for x as an integer
	long long R1 = (long long)ceil(R3);           // get minimum for x as a real number

	if (teach) {
		std::cout << "Since " << x1 << sq << " is always greater than, or equal to zero,\n";
		ShowEq(0, 0, NegDisc, 0, E1, F1, x1, y);
		printf(" must be less than, or equal to zero. \n"
			"This is verified in the segment limited by the roots. \n");
		if (N0 < 0) {
			std::cout << "The polynomial in " << y << " is always positive,";
			NoSol();
			return;
		}
		std::cout << "The roots are: (-" << par(E1) << " - llSqrt(" << numToStr(E1) << sq <<
			" - 4 * " << par(NegDisc) << " * " << par(F1) << ")) / (2 * " << par(NegDisc) << ") = " <<
			R3 << "\n";
		std::cout << "and: (-" << par(E1) << " + llSqrt(" << E1 << sq << " - 4 * " << par(NegDisc)
			<< " * " << par(F1) << ")) / (2 * " << par(NegDisc) + ") = " << R4 << " \n";
		if (R2<R1) {
			printf("There are no integers in this range,");
			NoSol();
			return;
		}

		std::cout << "All values of " << y << " from " << R1 << " to " << R2 << " should be replaced in \n";
		ShowEq(0, 0, NegDisc, 0, E1, F1, x1, y);
		printf(". The result should be the negative of a perfect square. \n");
		b = 0;  // flag used to control output formatting 
		for (u = R1; u <= R2; u++) {
			w1 = -NegDisc*u*u - E1*u - F1;
			w2 = llSqrt(w1);
			if (w2*w2 == w1) {
				if (b != 0) {
					/* not 1st solution */
					//w(", " + u);
					printf(", %lld", u);
				}
				else {
					/* first solution found*/
					std::cout << "The values of " << y << " are: " << u;
					b = 1;  /* set flag to show a solution was found */
				}
			}
		}
		if (b == 0) {
			std::cout << "This is not satisfied by any value of " << y;
			NoSol();
			return;
		}
		putchar('\n');
	}

	if (N0 < 0) {
		std::cout << "equation has no real solutions \n";
		NoSol();
		return;
	}

	/* try all integer values between R1 and R2*/
	bool found = false;
	for (u = R1; u <= R2; u++) {
		w1 = -NegDisc*u*u - E1*u - F1;
		w2 = llSqrt(w1);
		if (w2*w2 == w1) {
			/* we have found a solution */
			found = true;
			if (teach) {
				std::cout << y << " = " << u << "\n";
				std::cout << x1 << " = ";
				ShowLin(2 * g_A, g_B, g_D, x, y);
				std::cout << " = +/-Sqrt(" << (w2*w2) << ") = +/-" << w2 << "\n";
			}
			ShowElipSol(g_A, g_B, g_D, u, x, x1, y, w2);
			if (w2 != 0) {
				ShowElipSol(g_A, g_B, g_D, u, x, x1, y, -w2);
			}
		}
	}
	if (!found) {
		NoSol();
		return;
	}
}

/* determine type of equation. Also some basic checks for cases where there is no solution,
divide all coefficients by their gcd if gcd > 1, and calculate the discriminant */
equation_class	classify() {
	long long gcdA_E;

	/* get gcd of A, B, C, D, E. gcd = zero only if they are all zero*/
	gcdA_E = gcd(g_A, gcd(g_B, gcd(g_C, gcd(g_D, g_E))));
	if (teach) {
		std::cout << "First of all we must determine the gcd of all coefficients but the constant term." << "\n";
		std::cout << "that is : gcd (" << g_A << ", " << g_B << ", " << g_C << ", " << g_D << ", " << g_E << ") = "
			<< gcdA_E << "\n";
	}
	if (gcdA_E != 0) {
		/* protect against divide-by-zero error */
		if (g_F%gcdA_E != 0) {
			NoGcd(g_F);    // output message - no solution 
			return no_soln;
		}
		else {
			/* divide all coefficients by gcd */
			g_A /= gcdA_E;
			g_B /= gcdA_E;
			g_C /= gcdA_E;
			g_D /= gcdA_E;
			g_E /= gcdA_E;
			g_F /= gcdA_E;
			if (teach && (gcdA_E != 1)) {
				std::cout << divgcd;    // show new values after division by gcd
				ShowEq(g_A, g_B, g_C, g_D, g_E, g_F, "x", "y");
				printf(" = 0\n");
			}
		}
	}

	if (g_D == 0 && g_A != 0 && g_C != 0) 
		if (CheckMod(g_A, g_B, g_C, g_E, g_F))
			return no_soln;

	if (g_E == 0 && g_A != 0 && g_C != 0) 
		if (CheckMod(g_B, g_C, g_A, g_D, g_F))
			return no_soln;

	g_Disc = g_B*g_B - 4 * g_A*g_C; // get discriminant
	if (g_Disc > 0 &&
		llSqrt(g_Disc)*llSqrt(g_Disc) != g_Disc &&
		g_D == 0 && g_E == 0 && g_F != 0)
		return hyperbolic_homog;  /* Disc is not a perfect square  and D, E are 0, and F is non-zero.
		this is a type of homogeneous equation*/

	if (g_A == 0 && g_C == 0) {
		if (g_B == 0) {
			/* linear equation: A=B=C=0 */
			return linear;
		}
		else {
			/* simple hyperbolic; A = C = 0; B ≠ 0 */
			return simple_hyperbolic;
		}
	}

	/* not homogeneous, linear or simple hyperbolic equation */
	if (teach) {
		printf("We try now to solve this equation modulo 9, 16 and 25.\n");
	}
	if (Mod(9) || Mod(16) || Mod(25)) {
		return no_soln;
	}
	if (teach) {
		printf("There are solutions, so we must continue.\n");
	}

	if (g_Disc == 0)
		return parabolic;

	if (g_Disc < 0)
		return elliptical;

	return hyperbolic_gen;  // hyperbolic, not in other hyperbolic classes above
}

/* Solve Diophantine equations of the form: Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
This is the heart of the whole program.
there are six basic types which are solved using different methods:
•	Linear case:            A = B = C = 0.
•	Simple hyperbolic case: A = C = 0; B ≠ 0.
•	Elliptical case:        B^2 - 4AC < 0.
•	Parabolic case:			B^2 - 4AC = 0, A, B, C non-zero.
•	Hyperbolic case:        B^2 - 4AC > 0.
		(1) homogeneous equation Ax^2 + Bxy + Cy^2 + F = 0   i.e. D = 0, E = 0
		(2) general quadratic equation:  Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
Note: B^2 -4AC is known as the discriminant of the equation
*/
void solveEquation(const long long ax2, const long long bxy, const long long cy2,
	const long long dx, const long long ey, const long long f) {
	mpz_t Dp_A, Dp_B, Dp_C;
	bool teachaux;
	long long NegDisc, E1, F1, G, H, K, N0;
	long long T, g, h;
	std::string t1, x, y, x1, y1;

	/* code to check for heap corruption and memory leakage
	No code is generated unless compiled in debug mode */
	assert(_CrtCheckMemory());  //see https://msdn.microsoft.com/en-us/library/e73x0s4b.aspx
	_CrtMemCheckpoint(&memState1);  // see https://docs.microsoft.com/en-gb/cpp/c-runtime-library/reference/crtmemcheckpoint

	/* initialise some global variables */
	allSolsFound = false;
	NbrSols = 0;
	NbrCo = -1;
	also = false;
	g_A = ax2;           // x² coefficient
	g_B = bxy;           // xy coefficient 
	g_C = cy2;           // y² coefficient
	g_D = dx;            // x coefficient
	g_E = ey;            // y coefficient
	g_F = f;             // constant

	std::cout << "\nSolve Diophantine equation:   ";
	ShowEq(g_A, g_B, g_C, g_D, g_E, g_F, "x", "y");
	printf("  =  0\n by Dario Alejandro Alpern\n");

	equation_class eqnType = classify();   /* get equation type */

	if (g_A == 0 && g_C != 0) {
		T = g_A;   /* swap A and C */
		g_A = g_C;
		g_C = T;
		T = g_D;   /* swap D and E */
		g_D = g_E;
		g_E = T;
		ExchXY = true;
	}
	x = (ExchXY ? "y" : "x");
	y = (ExchXY ? "x" : "y");
	x1 = x + "'";
	y1 = y + "'";

	switch (eqnType) {
	case no_soln: 
		std::cout << msg;   // no solutions
		break; 

	case linear: { /* linear equation: A=B=C=0 */
		int rv = Linear(g_D, g_E, g_F);
		if (rv == 1) {
			std::cout << msg;   // no solutions
		}
		break;
	}

	case parabolic:   /* Parabolic case B^2 - 4AC = 0  A, B, C non-zero*/
		SolveParabolic(x, y, x1);
		break;

	case elliptical: {   /* B^2 - 4AC < 0 */
		SolveElliptical(x, x1, y);
		break;
	}

	case simple_hyperbolic:   // A = C = 0; B ≠ 0; therefore B^2 - 4AC > 0
		SolveSimpleHyperbolic();
		break;

	/* B^2 - 4AC > 0, D = 0, E = 0 */
	case hyperbolic_homog: {  

		/* solve using continued fractions */
		mpz_inits(Dp_A, Dp_B, Dp_C, NULL);  // local vars
		mpz_inits(Dp_NUM, Dp_DEN, Bi_H1, Bi_H2, Bi_K1, Bi_K2, Bi_L1, Bi_L2, NULL);   // global vars
		g_Disc = g_B*g_B - 4 * g_A*g_C; // get discriminant again

		LongToDoublePrecLong(g_A, Dp_A);   // Dp_A = A
		LongToDoublePrecLong(g_B, Dp_B);	// Dp_B = B
		LongToDoublePrecLong(g_C, Dp_C);	// Dp_C = C

		teachaux = teach;  // save value of teach
		if (abs(g_F) != 1) {
			teach = false;
		}
		GetRoot(Dp_A, Dp_B, Dp_C);       //  sneaky!! values returned in 
										 // global vars are used by SolContFrac
		teach = teachaux;    // restore saved value of teach

		G = H = g_F;
		K = 1;

		/* remove any double factors from G, add to K*/
		while (G % 4 == 0) {
			G /= 4;          // remove double 2
			K *= 2;
		}

		T = 3;             // remove odd factors
		while (abs(G) >= T*T) {
			while (G % (T*T) == 0) {
				G /= T*T;
				K *= T;
			}
			T += 2;   // advance to next odd number
		}

		/* K is the product of all the prime factors removed from G*/
		for (T = 1; T*T <= K; T++) {
			if (K%T == 0) {
				/* call SolContFrac for each factor of K  <= sqrt(K)*/
				SolContFrac(H, T, g_A, g_B, g_C, "");
			}
		}

		for (T = T - 1; T > 0; T--) {
			if (K%T == 0 && T*T < K) {
				/* call SolContFrac for each factor of K  > sqrt(K)*/
				SolContFrac(H, K / T, g_A, g_B, g_C, "");
			}
		}
		putchar('\n');

		if (also) {
			/* solutions found, so print them  */
			if (teach) {
				putchar('\n');
			}
			else {
				ShowAllLargeSolutions();
				printf("If (x,y) is a solution, (-x,-y) is also a solution.\n");
			}
			ShowRecursion(eqnType);
		}
		else
			std::cout << msg;   // no solutions

		mpz_clears(Dp_A, Dp_B, Dp_C, NULL);
		mpz_clears(Dp_NUM, Dp_DEN, Bi_H1, Bi_H2, Bi_K1, Bi_K2, Bi_L1, Bi_L2, NULL);
		break;             // finished processing homogeneous equation
	}

	/* none of the above */
	case hyperbolic_gen: {  
		/* B^2 - 4AC > 0. There are a couple of special cases:
		   (1)   B^2 - 4AC is a perfect square
		   (2) N0 (calculated below) is zero 
		   otherwise the equation is solved using continued fractions */

			NegDisc = 4 * g_A*g_C - g_B*g_B;   // get -ve of discriminant
			E1 = 4 * g_A*g_E - 2 * g_B*g_D;
			F1 = 4 * g_A*g_F - g_D*g_D;
			g = gcd(NegDisc, E1 / 2);
			g_CY1 = NegDisc / g;
			g_CY0 = E1 / 2 / g;
			N0 = g_CY0*g_CY0*g - g_CY1*F1;
			h = gcd(g_CY1, gcd(g, N0));
			if (teach) {
				printf("We want to convert this equation to one of the form:\n");
				std::cout << x1 << sq << " + B " << y << sq << " << C " << y << " + D = 0 \n";
				std::cout << "Multiplying the equation by " << par(4 * g_A) << ":\n";
				ShowEq(4 * g_A*g_A, 4 * g_A*g_B, 4 * g_A*g_C, 4 * g_A*g_D, 4 * g_A*g_E, 4 * g_A*g_F, x, y);
				printf(" = 0\n");
				ShowLin(4 * g_A*g_A, 0, 0, x + sq, y);
				if (g_B != 0 || g_D != 0) {
					printf(" + (");
					ShowLin(0, 4 * g_A*g_B, 4 * g_A*g_D, x, y);
					std::cout << ")" << x;
				}
				if (g_C != 0 || g_E != 0 || g_F != 0) {
					printf(" + (");
					ShowEq(0, 0, 4 * g_A*g_C, 0, 4 * g_A*g_E, 4 * g_A*g_F, x, y);
					printf(") = 0\n");
				}
				if (g_B != 0 || g_D != 0) {
					printf("To complete the square we should add and subtract:\n(");
					ShowLin(0, g_B, g_D, x, y);
					std::cout << ")" << sq << "\nThen the equation converts to:\n(";
					ShowLin(2 * g_A, g_B, g_D, x, y);
					std::cout << ")" << sq << " + (";
					ShowEq(0, 0, 4 * g_A*g_C, 0, 4 * g_A*g_E, 4 * g_A*g_F, x, y);
					printf(") - (");
					ShowEq(0, 0, g_B*g_B, 0, 2 * g_B*g_D, g_D*g_D, x, y);
					printf(") = 0 \n");
				}

				printf("(");
				ShowLin(2 * g_A, g_B, g_D, x, y);
				std::cout << ")" << sq << " + (";
				ShowEq(0, 0, NegDisc, 0, E1, F1, x, y);
				printf(") = 0\nNow we perform the substitution:\n");
				std::cout << x1 << " = ";
				ShowLin(2 * g_A, g_B, g_D, x, y);
				printf("\nThis gives:\n");
				ShowEq(1, 0, NegDisc, 0, E1, F1, x1, y);
				printf(" = 0 \n");
			}
			if (teach) {
				if (NegDisc != g*h) {
					//w("Multiplying the equation by " + numToStr(CY1 / h) + ":\n");
					printf("Multiplying the equation by %lld:\n", (g_CY1 / h));
					ShowEq(g_CY1 / h, 0, NegDisc*g_CY1 / h, 0, NegDisc*E1 / g / h, NegDisc*F1 / g / h, x1, y);
				}
				printf(" = 0 \n");
				if (E1 != 0) {
					Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
					ShowEq(g_CY1*g_CY1, 0, 0, 2 * g_CY0*g_CY1, 0, 0, y, "");
					printf(")");
					Show1(NegDisc*F1 / g / h, 1);
					printf(" = 0 \n");
					Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
					if (g_CY1 != 1) {
						std::cout << par(g_CY1) << sq << " ";
					}
					std::cout << y << sq << " + 2*";
					if (g_CY1 != 1) {
						std::cout << par(g_CY1) << "*";
					}
					std::cout << par(g_CY0) << " " << y << ")";
					Show1(NegDisc*F1 / g / h, 1);
					printf(" = 0 \n");
					std::cout << "Adding and subtracting " << (g == h ? "" : numToStr(g / h) + " * ")
						<< par(E1 / 2 / g) + sq + ":\n";
				}
				Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
				if (g_CY1 != 1) {
					std::cout << par(g_CY1) << sq << " ";
				}
				std::cout << y << sq << " + 2*";
				if (g_CY1 != 1) {
					std::cout << par(g_CY1) << "*";
				}
				std::cout << par(g_CY0) << " " << y << " + " << par(g_CY0) << sq << ")";
				Show1(NegDisc*F1 / g / h, 1);
				std::cout << " - " << (g == h ? "" : numToStr(g / h) + " * ")
					<< par(E1 / 2 / g) << sq << " = 0 \n";
				Show(g / h, "(", Show(g_CY1 / h, x1 + sq, 0));
				ShowLin(0, g_CY1, g_CY0, x, y);
				std::cout << ")" << sq;
				Show1(-N0 / h, 1);
				std::cout << " = 0 \nMaking the substitution " << y1 << " = ";
				ShowLin(0, g_CY1, g_CY0, x, y);
				printf(":\n");
				ShowLin(g_CY1 / h, g / h, -N0 / h, x1 + sq, y1 + sq);
				printf(" = 0 \n");
			}

			long long SqrtD = llSqrt(-NegDisc);
			long long N1 = abs(N0 / h);

			if (SqrtD*SqrtD == -NegDisc) {    // is discriminant a perfect square?
				SolveDiscIsSq(N0, x, y, SqrtD, g, h, N1, y1, x1);
				break;
			}

			/* discriminant is not a perfect square */
			if (N0 == 0) {
				/* solution is x=0 (y any int) or y=0 (x any int)*/
				ShowX1Y1(0, 0, g_A, g_B, g_D, NegDisc, E1 / 2);   // display solution
				break;
			}

			mpz_inits(Dp_A, Dp_B, Dp_C, NULL);  // local vars
			mpz_inits(Dp_NUM, Dp_DEN, Bi_H1, Bi_H2, Bi_K1, Bi_K2, Bi_L1, Bi_L2, NULL);   // global vars

			/* Solve by continued fractions.
			Firstly, test if we need two cycles or four cycles */
			LongToDoublePrecLong(g_A, Dp_A);     // Dp_A = A
			LongToDoublePrecLong(g_C, Dp_B);     // Dp_B = C
			Mult2LargeNumbers(Dp_A, Dp_B, Dp_C); // DP_C = A*C
			LongToDoublePrecLong(1, Dp_A);	     // Dp_A = 1
			LongToDoublePrecLong(g_B, Dp_B);	 // Dp_B = B
			GetRoot(Dp_A, Dp_B, Dp_C);           //  sneaky!! values returned in 
												 // global vars are used by ContFrac
			LongToDoublePrecLong(g_A, Dp_A);     // Dp_A = A

			ContFrac(Dp_A, 5, 1, 0, g_B*g_B - 4 * g_A*g_C, 1, g_A); /* A2, B2 solutions */

			g_Disc = g_B*g_B - 4 * g_A*g_C;        // calculate discrininant again
			G = (2 * g_A2 + g_B*g_B2) % g_Disc;
			H = (g_B*g_A2 + 2 * g_A*g_C*g_B2) % g_Disc;

			if (((g_C*g_D*(G - 2) + g_E*(g_B - H)) % g_Disc != 0 ||
				(g_D*(g_B - H) + g_A*g_E*(G - 2)) % g_Disc != 0) && ((g_C*g_D*(-G - 2) + g_E*(g_B + H)) % g_Disc != 0 ||
				(g_D*(g_B + H) + g_A*g_E*(-G - 2)) % g_Disc != 0)) {
				NbrCo *= 2;
				//std::cout << "**temp SolveEquation  NbrCo=" << NbrCo << "\n";
			}

			LongToDoublePrecLong(NegDisc / g / h, Dp_A);     // Dp_A = NegDisc/g/h
			LongToDoublePrecLong(0, Dp_B);                   // DP_B = 0:
			LongToDoublePrecLong(g / h, Dp_C);               // DP_C = g/h
			GetRoot(Dp_A, Dp_B, Dp_C);
			G = H = -N0 / h;
			K = 1;
			T = 3;

			/* remove any double factors from G, add to K*/
			while (G % 4 == 0) {
				G /= 4;
				K *= 2;
			}
			while (abs(G) >= T*T) {     // remove odd factors
				while (G % (T*T) == 0) {
					G /= T*T; K *= T;
				}
				T += 2;    // advance to next odd number
			}

			for (T = 1; T*T <= K; T++) {
				if (K%T == 0) {
					/* call SolContFrac for each divisor of K  <= sqrt(k) */
					SolContFrac(H, T, NegDisc / g / h, 0, g / h, "'");
				}
			}
			for (T = T - 1; T > 0; T--) {
				if (K%T == 0 && T*T < K) {
					/* call SolContFrac for each divisor of K  > sqrt(k) */
					SolContFrac(H, K / T, NegDisc / g / h, 0, g / h, "'");
				}
			}
			putchar('\n');
			if (also) {
				/* solutions found */
				if (!teach) {
					ShowAllLargeSolutions();
				}
				ShowRecursion(eqnType);
			}
			else
				std::cout << msg;		// no solutions

			mpz_clears(Dp_A, Dp_B, Dp_C, NULL);
			mpz_clears(Dp_NUM, Dp_DEN, Bi_H1, Bi_H2, Bi_K1, Bi_K2, Bi_L1, Bi_L2, NULL);
			break;
		}

	default:
		assert(1 == 0);			// cause breakpoint if compiled in debug mode
		abort();				// must not happen
	}
	
	assert(_CrtCheckMemory());        /* check for heap corruption &  */
	_CrtMemCheckpoint(&memState2); /* check for memory leakage*/
	if (_CrtMemDifference(&memStatDiff, &memState1, &memState2)) {
		std::cerr << "** there seems to have been a memory leak";
		_CrtMemDumpStatistics(&memStatDiff);
		_CrtDumpMemoryLeaks();
	}
}


//int Eval_Exception(int);

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
			std::cout << "Invalid Input! Please input a numerical value.  " ;
			std::cin.clear();               // clear error flags
			std::cin.ignore(100000, '\n');  // discard invalid input up to newline/enter
		}
	}
	return result;
}

/* Entry point for program. Get parameters and call solveEquation */
int main(int argc, char* argv[]) {
	try {
		long long int a, b, c, d, e, f;
		char yn = '\0';
		bool test = false;

		while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
			std::cout << "Enter own data (else run standard tests)? (Y/N): ";
			std::cin >> yn;
			test = (toupper(yn) == 'N');
		}

		if (!test) {
			std::cout << "Solve Diophantine equations of the form: Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0" << "\n";
			a = getnumber("Enter value for A ");
			b = getnumber("Enter value for B ");
			c = getnumber("Enter value for C ");
			d = getnumber("Enter value for D ");
			e = getnumber("Enter value for E ");
			f = getnumber("Enter value for F ");
			//std::cout << "solve " << a << "x² + " << t << "xy + " << cy2 << "y² + " << dx << "x + " << ey << "y + " << f << " = 0" << endl;
			printf("solve %lldx^2 + %lldxy + %lldy^2 + %lldx + %lldy + %lld = 0\n", a, b, c, d, e, f);

			yn = '\0';
			while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
				std::cout << "Detailed explanation required? (Y/N): ";
				std::cin >> yn;
				teach = (toupper(yn) == 'Y');
			}
			solveEquation(a, b, c, d, e, f);
			system("PAUSE");   // press any key to continue
			return EXIT_SUCCESS;
		}

		else {
			yn = '\0';
			while (toupper(yn) != 'Y' && toupper(yn) != 'N') {
				std::cout << "Detailed explanation required? (Y/N): ";
				std::cin >> yn;
				teach = (toupper(yn) == 'Y');
			}

			printf("Run standard tests\n");

			printf("\n1st test. No solution exists \n");
			a = 0; b = 0; c = 0; d = 10; e = 84; f = 15;
			solveEquation(a, b, c, d, e, f);

			/* example numbers below refer to Dario Alpert's document. */
			printf("\nExample 1 (linear)\n");
			/* x = -136 + 42t,  y = 16 - 5t,  where t is any integer number.
			equivalent to: x = -10 +42t     y = 1-5t      (substitute t+3 for t) */
			a = 0; b = 0; c = 0; d = 10; e = 84; f = 16;
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 2 (simple hyperbolic)\n");
			a = 0; b = 2; c = 0; d = 5; e = 56; f = 7;
			/*
			x = (266-56)/2  =  105, y = (266/266-5)/2    =  -2
			x = (-266-56)/2 = -161, y = [266/(-266)-5]/2 =  -3
			x = (2-56)/2    =  -27, y = (266/2-5)/2      =  64
			x = (-2-56)/2   =  -29, y = [266/(-2)-5]/2   = -69
			x = (38-56)/2   =   -9, y = (266/38-5)/2     =   1
			x = (-38-56)/2  =  -47, y = [266/(-38)-5]/2  =  -6
			x = (14-56)/2   =  -21, y = (266/14-5)/2     =   7
			x = (-14-56)/2  =  -35, y = [266/(-14)-5]/2  = -12
			The only 8 solutions to the equation are above
			*/
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 3 (elliptical)\n");
			a = 42; b = 8; c = 15; d = 23; e = 17; f = -4915;
			/* x = -11, y=-1 is the only solution */
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 4 (parabolic)\n");
			a = 8; b = -24; c = 18; d = 5; e = 7; f = 16;
			/*
			set 1: x = -174t^2 - 17t - 2, y = -116t^2 - 21t - 2
			set 2: x = -174t^2 - 41t - 4, y = -116t^2 - 37t - 4
			*/
			solveEquation(a, b, c, d, e, f);

			printf("\nProject Euler problem 140 (general hyperbolic)\n");
			a = -1; b = 0; c = 5; d = 0; e = 14; f = 1;
			solveEquation(a, b, c, d, e, f);
			/*
			basic solutions are:
			x=±1, 	y=0
			x=±2, 	y=-3
			x=±5, 	y=-4
			x=±7, 	y=2

			also, for each of the solutions above we can generate an infinite number of other
			solutions using the formulas:
			Xn+1 = PXn + QYn + K
			Yn+1 = RXn + SYn + L
			where:
			P = -9		    
			Q = 20		or -20
			K = 28		or -28
			R = 4		or -4
			S = -9
			L = -14
			*/

			printf("\nExample 5 (hyperbolic homogeneous)\n");
			//teach = true;
			a = 18; b = 41; c = 19; d = 0; e = 0; f = -24;
			/*
			basic solutions are:
			x = -10,							y =  6
			X = -7								Y =  11
			X = -202							Y =  312
			X = -14267							Y =  8751
			x = -10130							Y =  15646
			X = -284123							Y =  438834
			X =  -4 680127 (7 digits)			Y =  2 870666 (7 digits
			X =  -14 247838 (8 digits)			Y =  22 006088 (8 digits)
			If (x,y) is a solution, (-x,-y) is also a solution.
			*/
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 6\n");
			a = 3; b = 13; c = 5; d = -11; e = -7; f = -92;
			/* x=-4, y=0, or  x=2, y=3 etc
			P=8351, Q=32625, R=-19575, S=-76474, K=-28775, L=67450*/
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 6 - modified to produce 2 sets of values for P, Q, R, ...\n");
			/* x=-4, y=0, x=12, y=-2 etc
			   P= 8351, Q=32625, R=-19575, S=-76474, K=-775, L= 1825
			or P=-8351, Q=32625, R= 19575, S= 76474, K= 783, L=-1827*/
			a = 3; b = 13; c = 5; d = -11; e = -42; f = -92;
			
			solveEquation(a, b, c, d, e, f);

			printf("\nExample 7\n");
			/* x=4, y=7 etc
			P = -1188641 Q = -4979520 K = 5146869 R = 2489760 S = 10430239 L = -10780770 */
			a = 3; b = 14; c = 6; d = -17; e = -23; f = -505;
			solveEquation(a, b, c, d, e, f);

			a = 2; b = 5; c = 2; d = 6; e = 6; f = 4;
			std::cout << "\nproduct of two linear expressions\n";
			std::cout << "(general hyperbolic with discriminant a perfect square)";
			solveEquation(a, b, c, d, e, f);

			/* D == 0 && A != 0 && C != 0*/
			std::cout << "Hyperbolic homogeneous\n";
			a = 18; b = 41; c = 19; d = 0; e = 0; f = 13;
			solveEquation(a, b, c, d, e, f);

			std::cout << "\nElliptical, no solutions\n";
			a = -1; b = -1; c = -2; d = -3; e = -4; f = -5;
			solveEquation(a, b, c, d, e, f);

			std::cout << "\nElliptical, no solutions\n";
			a = 1; b = 0; c = 1; d = -0; e = 0; f = -6;
			solveEquation(a, b, c, d, e, f);
		}
		
		system("PAUSE");   // press any key to continue
		return EXIT_SUCCESS;
	}

	/* code below is executed if an exception occurs*/
	catch (const std::exception &exc)
	{
		// catch anything thrown within try block that derives from std::exception
		std::cerr << "\n** Exception:" << exc.what() << std::endl;
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");               // press any key to continue
		exit(EXIT_FAILURE);
	}

	catch (const char * str)
	{   
		std::cerr << "Caught exception: " << str << std::endl;
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");               // press any key to continue
		exit(EXIT_FAILURE);
	}

	catch (...) {   // catch block should only be executed under /EHa 
		/* most likely to be a SEH-type exception */
		//Eval_Exception(GetExceptionCode());
		std::cerr << "Caught unknown exception in catch(...)." << std::endl;
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");   // press any key to continue
		exit(EXIT_FAILURE);
	}
}

/* SEH exception handler */
//int Eval_Exception(int n_except) {
//	if (n_except != STATUS_INTEGER_OVERFLOW &&
//		n_except != STATUS_INTEGER_DIVIDE_BY_ZERO )   // Pass on most exceptions  
//		return EXCEPTION_CONTINUE_SEARCH;
//
//	// Execute some code to handle problems  
//	if (n_except == STATUS_INTEGER_OVERFLOW) {
//		std::cout << "\n** INTEGER_OVERFLOW - program continues ** \n";
//		Beep(1000, 1000);
//		return EXCEPTION_CONTINUE_EXECUTION;
//	}
//
//	if (n_except == STATUS_INTEGER_DIVIDE_BY_ZERO) {
//		std::cout << "\n** INTEGER DIVIDE BY ZERO - program terminates ** \n";
//		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
//		exit(EXIT_FAILURE);
//	}
//
//
//	return EXCEPTION_CONTINUE_EXECUTION;
//}
