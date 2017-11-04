#include "stdafx.h"

const std::string sq = "^2";
extern mpz_t Bi_L1, Bi_L2, Bi_H1, Bi_H2, Bi_K1, Bi_K2;
extern long long g_Disc;
extern long long g_A, g_B, g_D, g_F;
extern mpz_t Dp_NUM, Dp_DEN;
extern unsigned long long SqrtDisc;
extern long long g_CY1, g_CY0;
extern long long g_A2, g_B2;
extern int NbrSols, NbrCo;
extern bool teach;

long long g_A1, g_B1;
int NbrEqs, EqNbr;
std::string UU = "";
std::string VU = "";
std::string UL = "";
std::string VL = "";
std::string UL1 = "";
std::string VL1 = "";
std::string FP = "";

/*******************************************/
/* NextConv:                               */
/*  BigInteger tmp = Prev * A1 + Act * B1; */
/*  Act = Prev * A2 + Act * B2;            */
/*  Prev = Tmp;                            */
/*******************************************/
void NextConv(mpz_t Bi_Prev, mpz_t Bi_Act, const long long A1,
	const long long A2, const long long B1, const long long B2) {
	mpz_t t1, t2, tmp;

	/*std::cout << "**temp NextConv: Prev = ";
	ShowLargeNumber(Bi_Prev);
	std::cout << " Act =";
	ShowLargeNumber(Bi_Act);
	std::cout << " A1 =" << A1 << " A2=" << A2 << " B1=" << B1 << " B2=" << B2 << "\n";	*/

	mpz_inits(t1, t2, tmp, NULL);
	mpz_mul_si(t1, Bi_Prev, A1);   // t1 = Prev * A1
	mpz_mul_si(t2, Bi_Act, B1);    // t2 = Act * B1
	mpz_add(tmp, t1, t2);          // tmp = Prev * A1 + Act * B1

	mpz_mul_si(t1, Bi_Prev, A2);   // t1 = Prev * A2
	mpz_mul_si(t2, Bi_Act, B2);    // t2 = Act * B2
	mpz_add(Bi_Act, t1, t2);       // act = Prev * A2 + Act * B2

	mpz_set(Bi_Prev, tmp);        // Prev = Tmp

	mpz_clears(t1, t2, tmp, NULL);   // avoid memory leakage

	/* temporary */
	/*std::cout << "**exit NextConv: Prev = ";
	ShowLargeNumber(Bi_Prev);
	std::cout << " Act =";
	ShowLargeNumber(Bi_Act);
	std::cout << "\n\n";*/
	/* end temporary */
}

/* uses global variables Bi_L1, Bi_L2, g_A, g_B, g_D, g_F, g_CY0, g_CY1*/
bool ShowHomoSols(int type, mpz_t Bi_SHH, mpz_t Bi_SHK, long long s, long long T,
	const long long MagnifY, const std::string eqX, const std::string eqY) {

	assert(_CrtCheckMemory());

	/*std::cout << "**temp ShowHomoSols: s=" << s << "  T=" << T << "  MagnifY=" << MagnifY << "\n";
	std::cout << "Bi_SH1="; ShowLargeNumber(Bi_SH1);
	std::cout << "  Bi_SHK1="; ShowLargeNumber(Bi_SHK1); std::cout << "\n";*/

	int i;
	std::string U = (type == 4 ? "'" : "");
	std::string X1 = "X";
	std::string Y1 = (MagnifY == 1 ? "Y" : "Y'");

	if (teach) {
		ShowLargeXY(Y1, "Z", Bi_SHH, Bi_SHK, true, eqX, eqY);
		std::cout << "Since " << X1 << U << "0 = ";
		if (T>1) {
			std::cout << T << " " << UU << U << "0 = " << T << " (";
		}
		ShowLin(s, -g_F, 0, VU + "0", "Z0");
		if (T>1) {
			std::cout << ")\nand " << Y1 << "0 = " << T << " V0";
		}
		printf(":\n");
	}
	MultAddLargeNumbers(s*T, Bi_SHH, -g_F*T, Bi_SHK, Bi_L1);  // result stored in Bi_L1
	MultLargeNumber(T, Bi_SHH, Bi_L2);

	/*std::cout << "**temp ShowHomoSols(1) Bi_L1="; 	ShowLargeNumber(Bi_L1);
	std::cout << "  Bi_L2="; ShowLargeNumber(Bi_L2); std::cout << "\n";*/

	if (type == 4) {
		if (teach) {
			ShowLargeXY(X1 + U, Y1 + U, Bi_L1, Bi_L2, false, "", "");
		}
		for (i = (IsZero(Bi_L1) && IsZero(Bi_L2) ? 1 : 0); i<2; i++) {
			AddLarge(Bi_L2, -g_CY0, Bi_L2);
			//std::cout << "**temp ShowHomoSols(2) Bi_L2="; 	ShowLargeNumber(Bi_L2); std::cout << "\n";

			if (teach) {
				std::cout << Y1 << "0 = (";
				ShowLin(0, i == 0 ? 1 : -1, -g_CY0, "", Y1 + "0");
				std::cout << ")/" << par(g_CY1) << "\n";
			}
			if (tDivLargeNumber(Bi_L2, g_CY1, Bi_L2) != 0) {
				//std::cout << "**temp ShowHomoSols(3) Bi_L2="; 	ShowLargeNumber(Bi_L2); std::cout << "\n";
				if (teach) {
					printf("It is not an integer number. \n");
				}
			}
			else {
				//std::cout << "**temp ShowHomoSols(4) Bi_L2="; 	ShowLargeNumber(Bi_L2); std::cout << "\n";
				if (teach) {
					std::cout << X1 << "0 = (";
					ShowLin(i == 0 ? 1 : -1, -g_B, -g_D, X1 + "0", Y1 + "0");
					std::cout << ") / " << par(2 * g_A) << "\n";
				}
				AddLarge(Bi_L1, -g_D, Bi_L1);  // L1 -= D
				MultAddLargeNumbers(1, Bi_L1, -g_B, Bi_L2, Bi_L1);  // store result in L1
				//std::cout << "**temp ShowHomoSols(5) Bi_L1="; 	ShowLargeNumber(Bi_L1); std::cout << "\n";
				if (tDivLargeNumber(Bi_L1, 2 * g_A, Bi_L1) != 0) {
					//std::cout << "**temp ShowHomoSols(6) Bi_L1="; 	ShowLargeNumber(Bi_L1); std::cout << "\n";
					if (teach) {
						printf("It is not an integer number. \n");
					}
				}
				else {
					//std::cout << "**temp ShowHomoSols(7) Bi_L1="; 	ShowLargeNumber(Bi_L1); std::cout << "\n";
					if (teach) {
						if (MagnifY != 1) {
							ShowLargeXY(X1, Y1, Bi_L1, Bi_L2, false, "", "");
							std::cout << "Since Y = " << MagnifY << Y1;
						}
						putchar('\n');
					}
					MultLargeNumber(MagnifY, Bi_L2, Bi_L2);
					ShowLargeXY("X", "Y", Bi_L1, Bi_L2, false, "", "");

					/*std::cout << "**temp ShowHomoSols(8) returns true: Bi_L1="; ShowLargeNumber(Bi_L1);
					std::cout << "  Bi_L2="; ShowLargeNumber(Bi_L2); std::cout << "\n";*/

					return true;
				}
			}
			//gmp_printf("**temp ShowHomoSols(9) Bi_L1= %lld*%Zd + %lld*%Zd =", -s*T, Bi_SH1, g_F*T, Bi_SHK1);

			MultAddLargeNumbers(-s*T, Bi_SHH, g_F*T, Bi_SHK, Bi_L1);
			MultLargeNumber(-T, Bi_SHH, Bi_L2);

			/*ShowLargeNumber(Bi_L1);
			std::cout << "  Bi_L2="; ShowLargeNumber(Bi_L2); std::cout << "\n";*/
		}
	}
	else {
		if (teach) {
			if (MagnifY != 1) {
				ShowLargeXY(X1, Y1, Bi_L1, Bi_L2, false, "", "");
				std::cout << "Since " << Y1 << " = " << MagnifY << "Y";
			}
			putchar('\n');
		}
		MultLargeNumber(MagnifY, Bi_L2, Bi_L2);
		
		if (teach) {
			ShowLargeXY("X", "Y", Bi_L1, Bi_L2, false, "", "");
			ChangeSign(Bi_L1);
			ChangeSign(Bi_L2);
			ShowLargeXY("X", "Y", Bi_L1, Bi_L2, false, "", "");
			ChangeSign(Bi_L1);
			ChangeSign(Bi_L2);
		}

		/*std::cout << "**temp ShowHomoSols(10) returns true: Bi_L1=";
		ShowLargeNumber(Bi_L1);
		std::cout << "  Bi_L2=";
		ShowLargeNumber(Bi_L2);
		std::cout << "\n";*/

		return true;
	}
	//std::cout << "**temp ShowHomoSols(11) returns false \n";
	return false;   // solution not found
}

/*************************************************************************** 
* type = 1: Find convergents                                               *
* type = 2: Find convergents for x^2 + Bxy + ACy^2 = 1 (recursion)         *
* type = 3: Find convergents for modified equation in homogeneous equation *
* type = 4: Find convergents for modified equation in complete solution    *
* type = 5: Find convergents for x^2 + Bxy + ACy^2 = 1 (mod B^2-4AC)       *
* returns true if there are solutions, otherwise false                     *                      
* uses global variables DP_NUM, DP_DEN, Bi_H1, Bi_H2, Bi_K1, Bi_K2, NbrCo  *
*     NbrEqs, Eqnbr, NbrSols, g_F, g_A1, g_A2, g_B1, g_B2                  *    
****************************************************************************/
bool ContFrac(const mpz_t Dp_A, int type, const int SqrtSign, long long s, long long T,
	long long MagnifY, long long A) {

	/*std::cout << "**temp ContFrac: type=" << type << "  SqrtSign=" << SqrtSign << "  s=" << s;
	std::cout << "  T=" << T << "  MagnifY=" << MagnifY << "  A=" << A << "\n";
	std::cout << "Dp_NUM=" << numToString(Dp_NUM);
	std::cout << "  Dp_DEN=" << numToString(Dp_DEN) << "\n";*/

	long long P, Z, M, P1, M1, Tmp, K, L, Mu;
	mpz_t Dp_P,  Dp_M, Dp_Z, Dp_G, Dp_Mu, Dp_K, Dp_L, Dp_M1, Dp_P1, Dp_zz;
	long long H1ModCY1 = 1, H2ModCY1 = 0, K1ModCY1 = 0, K2ModCY1 = 1;
	bool Sols = true, secondDo = true;
	int Conv;
	int Co = -1;
	std::string U = (type == 4 ? "'" : "");
	std::string X1 = "X";
	std::string Y1 = (MagnifY == 1 ? "Y" : "Y'");

	assert(_CrtCheckMemory());

	mpz_inits(Dp_P,Dp_M, Dp_Z, Dp_G, Dp_Mu, Dp_K, Dp_L, Dp_M1, Dp_P1, Dp_zz, NULL);

	if (IsOne(Dp_A)) {
		/* Dp_A = 1 */
		mpz_set_si(Bi_H1, SqrtSign);
		mpz_set_si(Bi_K1, 0);
		/*std::cout << "**temp ContFrac Bi_H1="; ShowLargeNumber(Bi_H1);
		std::cout << "  Bi_K1="; ShowLargeNumber(Bi_K1); std::cout << "\n";*/

		if (type == 1) {
			ShowLargeXY(X1, Y1, Bi_H1, Bi_K1, true, "", "");
			mpz_clears(Dp_P, Dp_M, Dp_Z, Dp_G, Dp_Mu, Dp_K, Dp_L, Dp_M1, Dp_P1, Dp_zz, NULL);
			//std::cout << "** temp ContFrac returns true (1) - solution(s) found\n";
			assert(_CrtCheckMemory());
			return true;            /* Indicate there are solutions */
		}

		if ((type == 3 || type == 4) && (g_Disc != 5 || A*g_F<0)) {
			if (ShowHomoSols(type, Bi_H1, Bi_K1, s, T, MagnifY, "", "")) {
				mpz_clears(Dp_P, Dp_M, Dp_Z, Dp_G, Dp_Mu, Dp_K, Dp_L, Dp_M1, Dp_P1, Dp_zz, NULL);
				//std::cout << "**temp ContFrac(2) - solution found\n";
				assert(_CrtCheckMemory());
				return true;          /* Indicate there are solutions */
			}
		}
	}

	/* Paso = 1: Quick scan for solutions */
	/* Paso = 2: Show actual solutions */
	for (int Paso = (type == 2 || g_Disc == 5 && A*g_F>0 && (type == 3 || type == 4) ? 2 : 1);
		Sols && Paso <= 2; Paso++) {
		Conv = 0;
		Sols = false;
		
		mpz_set(Dp_P, Dp_DEN);   // P = DEN
		if (SqrtSign < 0) {
			ChSign(Dp_P);
		}
		LongToDoublePrecLong(SqrtDisc + (IsNeg(Dp_P) ? 1 : 0), Dp_K);
		//std::cout << "  Dp_K=" << numToString(Dp_K) << "(1)\n";
		if (SqrtSign < 0) {
			SubtDoublePrecLong(Dp_K, Dp_NUM, Dp_K);  // K -= NUM
			//std::cout << "  Dp_K=" << numToString(Dp_K) << "(2)\n";
		}
		else {
			AddLarge(Dp_K, Dp_NUM, Dp_K);      // K += NUM
			//std::cout << "  Dp_K=" << numToString(Dp_K) << "(3)\n";
		}

		//std::cout << "  Dp_DEN=" << numToString(Dp_DEN);
		Z = DivDoublePrec(Dp_K, Dp_P);           // Z = K/P
		//std::cout << "  Dp_K=" << numToString(Dp_K) << "  Dp_P=" << numToString(Dp_P) << "\n";
		LongToDoublePrecLong(Z, Dp_M);           // M = Z (=K/P)
		Mult2LargeNumbers(Dp_M, Dp_DEN, Dp_K);  // K = M*DEN
		//std::cout << "  Dp_K=" << numToString(Dp_K) << " Z=" << Z;

		SubtDoublePrecLong(Dp_K, Dp_NUM, Dp_M);      // M = K-NUM
		//std::cout << "  Dp_M=" << numToString(Dp_M) << "(4)\n";

		if (SqrtSign < 0) {
			ChSign(Dp_M);
		}

		/* type = 4: Find convergents for modified equation in complete solution */
		if (type == 4) {
			H2ModCY1 = Z%g_CY1;
		}

		/* type = 5: Find convergents for x^2 + Bxy + ACy^2 = 1 (mod B^2-4AC) */
		if (type == 5) {
			g_A1 = g_B2 = 1;
			g_A2 = Z%T;
			g_B1 = 0;
		}
		else {
			mpz_set_si(Bi_H1, SqrtSign);
			mpz_set_si(Bi_H2, Z*SqrtSign);
			mpz_set_si(Bi_K1, 0);
			mpz_set(Bi_K2, Bi_H1);
			/*std::cout << "**temp ContFrac Bi_H1="; ShowLargeNumber(Bi_H1);
			std::cout << "  Bi_H2="; ShowLargeNumber(Bi_H2);
			std::cout << "  Bi_K1="; ShowLargeNumber(Bi_K1);
			std::cout << "  Bi_K2="; ShowLargeNumber(Bi_K2); std::cout << "\n";*/

			g_A1 = g_B2 = 1;
			g_A2 = g_B1 = 0;
		}

		Co = -1;

		mpz_set_si(Dp_K, -1);  // K = -1
		mpz_set_si(Dp_L, -1);  // L = -1
		//std::cout << "  Dp_K=" << numToString(Dp_K) << "(6)\n";
		/* set Mu */
		switch (type) {
		case 1:   // find convergents 
			LongToDoublePrecLong(-2 * g_F*SqrtSign, Dp_Mu);
			break;

		case 3:     // find convergents for modified equation in homogeneous equation
		case 4:     // find convergents for modified equation in complete solution
			LongToDoublePrecLong(-2 * SqrtSign, Dp_Mu);
			break;

		default:  // type = 2 or 5
			mpz_set(Dp_Mu, Dp_DEN); // Mu = DEN
			if (SqrtSign > 0) {
				ChSign(Dp_Mu);
				//std::cout << "**temp ContFrac(3A) Dp_Mu= " << numToString(Dp_Mu) << "  (revsign)\n";
			}
		}

		/*std::cout << "**temp ContFrac(3B) Dp_Mu= " << numToString(Dp_Mu);
		std::cout << "  Dp_K=" << numToString(Dp_K);
		std::cout << "  Dp_P1=" << numToString(Dp_P1) << "\n";*/

		do {
			LongToDoublePrecLong(g_Disc, Dp_Z);     // Z = Disc  
			//std::cout << "**temp ContFrac  Dp_Z=" << numToString(Dp_Z) << "(3D)\n";

			Mult2LargeNumbers(Dp_M, Dp_M, Dp_G);  // G = M*M
			//std::cout << "**temp ContFrac  Dp_G=" << numToString(Dp_G) << "\n";

			SubtDoublePrecLong(Dp_Z, Dp_G, Dp_G);    // G = Z-G  = Disc -M*M
			//std::cout << "**temp ContFrac  Dp_G=" << numToString(Dp_G) << "\n";

			DivideDoublePrecLong(Dp_G, Dp_P, Dp_P1);  // P1 = (Disc-M*M)/P
			//std::cout << "**temp ContFrac Dp_P1=" << numToString(Dp_P1) << "\n";
			/* Z = SqrtDisc +(1 or 0, depending on sign of P1) */
			LongToDoublePrecLong(SqrtDisc + (IsNeg(Dp_P1) ? 1 : 0), Dp_Z);
			AddLarge(Dp_Z, Dp_M, Dp_K);      // K = M+Z
			//std::cout << "**temp ContFrac  Dp_K=" << numToString(Dp_K) << "\n";
			/* round Z to a multiple of P1 */

			Z = DivDoublePrec(Dp_K, Dp_P1);        // Z = K/P1
			LongToDoublePrecLong(Z, Dp_G);         // G = Z = K/P1
			Mult2LargeNumbers(Dp_G, Dp_P1, Dp_Z); // Z=G*P1  i.e Z rounded to multiple of P1
			//std::cout << "**temp ContFrac  Dp_Z=" << numToString(Dp_Z) << "\n";

			SubtDoublePrecLong(Dp_Z, Dp_M, Dp_M1); // M1 = Z-M
			//std::cout << " **temp ContFrac Dp_M1=" << numToString(Dp_M1) << "\n";

			mpz_add_si(Dp_zz, Dp_M, SqrtDisc);  // Dp_zz = SqrtDisc + DP_M
			if (Co<0 &&
				(mpz_cmp_si(Dp_P, 0)        >  0) &&
				(mpz_cmp   (Dp_P, Dp_zz)    <= 0) &&
				(mpz_cmp_si(Dp_M, 0)        >  0) &&
				(mpz_cmp_si(Dp_M, SqrtDisc) <= 0) ) 	{
				/* if (P > 0 & P <= SqrtDisc+M & M > 0 & M < SqrtDisc) */

				Co = 0;
				mpz_set(Dp_K, Dp_P);   // K = P
				mpz_set(Dp_L, Dp_M);   // L = M
				//std::cout << "**temp ContFrac(4) Dp_K=" << numToString(Dp_K);
				//std::cout << "  Dp_L=" << numToString(Dp_L) << "\n";
			}

			/*std::cout << "**temp ContFrac(5)  type=" << type << "  Co=" << Co;
			std::cout << "  Dp_P=" << numToString(Dp_P);
			std::cout << "  Dp_Mu=" << numToString(Dp_Mu);
			std::cout << "  Dp_K=" << numToString(Dp_K);
			std::cout << "  Dp_P1=" << numToString(Dp_P1) << "\n";*/

			if (type == 1 && AreEqual(Dp_P, Dp_Mu)	) { 
				// Solution found
				if (Co % 2 == 0 || !AreEqual(Dp_K, Dp_P1)) {
					if (Paso == 2) {
						if (g_A2 != 0) {
							NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
							NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
							g_A1 = g_B2 = 1;
							g_A2 = g_B1 = 0;
						}
						ShowLargeXY(X1, Y1, Bi_H2, Bi_K2, true, "NUM(" + numToStr(Conv) + ") = ",
							"DEN(" + numToStr(Conv) + ") = ");
					}
					Sols = true;
				}
				secondDo = false;
				//std::cout << "**temp ContFrac(5A)\n";
				break;
			}

			if (type == 3 || type == 4) {
				if (Co == 0 && A*g_F>0 && g_Disc == 5) {  /* Solution found */
					if (Paso == 1) {
						secondDo = false;
						Sols = true;
						//std::cout << "**temp ContFrac(6)\n";
						break;
					}
					else {
						NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
						NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
						g_A1 = g_B2 = 1;
						g_A2 = g_B1 = 0;
						ChangeSign(Bi_H2);
						AddLarge(Bi_H1, Bi_H2, Bi_H2);
						ChangeSign(Bi_H2);
						ChangeSign(Bi_K2);
						AddLarge(Bi_K1, Bi_K2, Bi_K2);
						ChangeSign(Bi_K2);
						if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY,
							"NUM(" + numToStr(Conv) + ") - NUM(" + numToStr(Conv - 1) + ") = ",
							"DEN(" + numToStr(Conv) + ") - DEN(" + numToStr(Conv - 1) + ") = ")) {
							secondDo = false;
							Sols = true;
							//std::cout << "**temp ContFrac(7) - solution found\n";
							break;
						}
						AddLarge(Bi_H1, Bi_H2, Bi_H2);
						AddLarge(Bi_K1, Bi_K2, Bi_K2);
					}
				}
				if (AreEqual(Dp_P1, Dp_Mu)) {  
					// Solution found
					if (Co % 2 == 0 || !AreEqual(Dp_K, Dp_P1) || !AreEqual(Dp_L, Dp_M1)	) {
						if (Paso == 2) {
							if (g_A2 != 0) {
								NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
								NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
								g_A1 = g_B2 = 1;
								g_A2 = g_B1 = 0;
							}
							if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY, "NUM(" + numToStr(Conv) + ") = ",
								"DEN(" + numToStr(Conv) + ") = ")) {
								secondDo = false;
								Sols = true;
								//std::cout << "**temp ContFrac(8) - solution found\n";
								break;
							}
						}
						else {
							if (type == 4) {
								Tmp = H2ModCY1*T;
								if ((Tmp - g_CY0) % g_CY1 == 0 || (Tmp + g_CY0) % g_CY1 == 0) {
									secondDo = false;
									Sols = true;
									//std::cout << "**temp ContFrac(9)\n";
									break;
								}
							}
							else {
								secondDo = false;
								Sols = true;
								//std::cout << "**temp ContFrac(10)\n";
								break;
							}
						}
					}
				}

				if (Paso == 1 && type == 4) {
					Tmp = (H1ModCY1 + Z*H2ModCY1) % g_CY1;
					H1ModCY1 = H2ModCY1;
					H2ModCY1 = Tmp;
					Tmp = (K1ModCY1 + Z*K2ModCY1) % g_CY1;
					K1ModCY1 = K2ModCY1;
					K2ModCY1 = Tmp;
				}
			}
			mpz_set(Dp_M, Dp_M1);   // M = M1
			mpz_set(Dp_P, Dp_P1);   // P = P1
			if (Co == 0) {
				Co = 1;
			}
			if (type == 5) {
				Tmp = (g_A1 + Z*g_A2) % T;
				g_A1 = g_A2; g_A2 = Tmp;
				Tmp = (g_B1 + Z*g_B2) % T;
				g_B1 = g_B2; g_B2 = Tmp;
			}
			ChSign(Dp_Mu);
			if (Paso == 2) {
				if (g_A2 != 0 && Z>(quintillion / 10 - g_A1) / g_A2 || 
					g_B2 != 0 && Z>(quintillion / 10 - g_B1) / g_B2) {
					NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
					NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
					g_A1 = g_B2 = 1;
					g_A2 = g_B1 = 0;
				}
				g_A1 += Z*g_A2;
				g_B1 += Z*g_B2;
				Tmp = g_A1; g_A1 = g_A2; g_A2 = Tmp;   // swap A1 and A2
				Tmp = g_B1; g_B1 = g_B2; g_B2 = Tmp;   // swap B1 and B2
			}
			Conv++;
		} while (Co<0);

		if (!secondDo) {
			continue;    // go to next step (paso = 2)
		}

		Mu = DoublePrecToLong(Dp_Mu);
		L = DoublePrecToLong(Dp_L);
		K = DoublePrecToLong(Dp_K);
		M = DoublePrecToLong(Dp_M);
		P = DoublePrecToLong(Dp_P);

		do {
			P1 = (g_Disc - M*M) / P;    /* P & Q should be > 0 (See Knuth Ex 4.5.3-12) */
			Z = (SqrtDisc + M) / P1;
			M1 = Z*P1 - M;
			if (type == 1 && P == Mu) {    /* Solution found */
				if (Co % 2 == 0 || K != P1 || L != M1) {
					if (Paso == 2) {
						if (g_A2 != 0) {
							NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
							NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
							g_A1 = g_B2 = 1;
							g_A2 = g_B1 = 0;
						}
						ShowLargeXY(X1, Y1, Bi_H1, Bi_K1, true, "NUM(" + numToStr(Conv) + ") = ",
							"DEN(" + numToStr(Conv) + ") = ");
					}
					Sols = true;
				}
				//std::cout << "**temp ContFrac(11)\n";
				break;
			}
			if (type == 3 || type == 4) {
				if ((Co & 1) == 0 && A*g_F>0 && g_Disc == 5) {   /* Solution found */
					if (Paso == 1) {
						Sols = true;
						//std::cout << "**temp ContFrac(12)\n";
						break;
					}
					else {
						NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
						NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
						g_A1 = g_B2 = 1;
						g_A2 = g_B1 = 0;
						ChangeSign(Bi_H2);
						AddLarge(Bi_H1, Bi_H2, Bi_H2);
						ChangeSign(Bi_H2);
						ChangeSign(Bi_K2);
						AddLarge(Bi_K1, Bi_K2, Bi_K2);
						ChangeSign(Bi_K2);
						if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY,
							"NUM(" + numToStr(Conv) + ") - NUM(" + numToStr(Conv - 1) + ") = ",
							"DEN(" + numToStr(Conv) + ") - DEN(" + numToStr(Conv - 1) + ") = ")) {
							Sols = true;
							//std::cout << "**temp ContFrac(13) - solution found\n";
							break;
						}
						AddLarge(Bi_H1, Bi_H2, Bi_H2);  // H2 = H1 + H2
						AddLarge(Bi_K1, Bi_K2, Bi_K2);  // K2 = K1 + K2
					}
				}
				if (P1 == Mu) {   /* Solution found */
					if (Co % 2 == 0 || K != P1 || L != M1) {
						if (Paso == 2) {
							if (g_A2 != 0) {
								NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
								NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
								g_A1 = g_B2 = 1;
								g_A2 = g_B1 = 0;
							}
							if (ShowHomoSols(type, Bi_H2, Bi_K2, s, T, MagnifY,
								"NUM(" + numToStr(Conv) + ") = ",
								"DEN(" + numToStr(Conv) + ") = ")) {
								Sols = true;
								//std::cout << "**temp ContFrac(14) - solution found\n";
								break;
							}
						}
						else {
							if (type == 4) {
								Tmp = H2ModCY1*T;
								if ((Tmp - g_CY0) % g_CY1 == 0 || (Tmp + g_CY0) % g_CY1 == 0) {
									Sols = true;
									break;
								}
							}
							else {
								Sols = true;
								//std::cout << "**temp ContFrac(15)\n";
								break;
							}
						}
					}
				}

				if (Paso == 1 && type == 4) {
					Tmp = (H1ModCY1 + Z*H2ModCY1) % g_CY1;
					H1ModCY1 = H2ModCY1;
					H2ModCY1 = Tmp;
					Tmp = (K1ModCY1 + Z*K2ModCY1) % g_CY1;
					K1ModCY1 = K2ModCY1;
					K2ModCY1 = Tmp;
				}
			}

			Co++;
			if (Co % 5000 == 0) {
				std::cout << "Conv: " << Co << " (Eq " << EqNbr << " of " << NbrEqs << ") \n";
				std::cout << NbrSols << " solution" << (NbrSols == 1 ? "\n" : "s\n");
			}
			if (Co % 5000 == 2500) {
				std::cout << NbrSols << " solution" << (NbrSols == 1 ? "\n" : "s\n");
			}

			M = M1;
			P = P1;
			if (type == 2 && P1 == Mu) {
				NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
				NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
				Sols = true;
				//std::cout << "**temp ContFrac(17)\n";
				break;
			}
			if (type == 5) {
				if (P1 == Mu) {
					NbrCo = Co;
					Sols = true;
					//std::cout << "**temp ContFrac(18) NbrCo=" << NbrCo << "\n";
					break;
				}
				else {
					Tmp = (g_A1 + Z*g_A2) % T;
					g_A1 = g_A2; g_A2 = Tmp;
					Tmp = (g_B1 + Z*g_B2) % T;
					g_B1 = g_B2; g_B2 = Tmp;
				}
			}
			Mu = -Mu;
			if (Paso == 2) {
				if (g_A2 != 0 && Z>(quintillion / 10 - g_A1) / g_A2 || g_B2 != 0 && Z>(quintillion / 10 - g_B1) / g_B2) {
					NextConv(Bi_H1, Bi_H2, g_A1, g_A2, g_B1, g_B2);
					NextConv(Bi_K1, Bi_K2, g_A1, g_A2, g_B1, g_B2);
					g_A1 = g_B2 = 1;
					g_A2 = g_B1 = 0;
				}
				g_A1 += Z*g_A2; g_B1 += Z*g_B2;
				Tmp = g_A1; g_A1 = g_A2; g_A2 = Tmp;   // swap A1 and A2
				Tmp = g_B1; g_B1 = g_B2; g_B2 = Tmp;   // swap B1 and B2
			}
			Conv++;
			/*std::cout << "**temp ContFrac(18A) NbrCo=" << NbrCo << "  Co=" << Co;
			std::cout << "  K=" << K << "  P=" << P << "  L=" << L << "  M=" << M << "\n";*/
		} while (NbrCo>0 ? Co != NbrCo : Co % 2 != 0 || K != P || L != M);

		//std::cout << "**temp ContFrac(18B)\n";
		/* type = 5: Find convergents for x^2 + Bxy + ACy^2 = 1 (mod B^2-4AC) */
		if (type == 5) {
			//std::cout << "**temp ContFrac(19)\n";
			break;   // break out of for loop
		}
	}                        /* end for */

	//std::cout << "**temp ContFrac returns " << Sols << "\n";
	mpz_clears(Dp_P, Dp_M, Dp_Z, Dp_G, Dp_Mu, Dp_K, Dp_L, Dp_M1, Dp_P1, Dp_zz, NULL);
	assert(_CrtCheckMemory());
	return Sols;
}

/******************************************************************************
* H = Constant term                                                           *
* T = Divisor of the square part of the constant term                         *
* A = X^2 coefficient                                                         *
* B = XY coefficient                                                          *
* C = Y^2 coefficient                                                         *
* SCFstr = Nothing or apostrophe (complete quad equation)                     *
* called from solveEquation                                                   *
* overwrites global variable g_F, Bi_L1, Bi_L2, DP_NUM, DP_DEN, disc,         *
*            sqrtdisc, NbrSols                                                *
*******************************************************************************/
void SolContFrac(long long H, long long T, long long A, long long B, long long C, 
	std::string SCFstr) {
	long long factor[64] = { 0 };
	long long P[64] = { 0 };
	long long Q[64] = { 0 };
	long long Dif[64] = { 0 };   /* Holds difference */
	long long mod[64] = { 0 };
	long long pos[64] = { 0 };
	long long Tmp, q, s, t, v, Pp, dif, Sol1, Sol2, Modulo;
	long long Tmp1 = SqrtDisc;
	mpz_t Dp_Tmp2, Dp_Tmp3;
	long long Tmp4 = g_Disc;
	long long ValA, ValB, ValC, ValF, ValAM, ValBM, ValCM;
	long long VarD, VarK, VarQ, VarR, VarV, VarW, VarX, VarY, VarY1;
	mpz_t Dp_A, Dp_B, Dp_C, Dp_R, Dp_S, Dp_T;
	mpz_t Bi_Xcopy, Bi_Ycopy;
	int index, index2, cont;
	int NbrFactors;
	long long gcdAF, MagnifY;
	int cuenta = 0;
	long long OrigA, OrigC;
	bool ShowHR = false;

	mpz_inits(Dp_A, Dp_B, Dp_C, Dp_R, Dp_S, Dp_T, Dp_Tmp2, Dp_Tmp3, Bi_Xcopy, Bi_Ycopy, NULL);

	//std::cout << "**temp SolContFrac:  H=" << H << "  T=" << T << "  A=" << A << "  B=" << B << "  C=" << C << "\n";;

	mpz_set(Dp_Tmp2, Dp_NUM);          // copy global NUM to local Tmp2
	mpz_set(Dp_Tmp3, Dp_DEN);          // copy global DEN to l.ocal Tmp3

	g_F = H / T / T;
	if (teach && T>1) {
		std::cout << "Since " << T << " * " << T << " is a divisor of the constant term ("
			<< H << "), the solutions should be " << T << " times the solutions of ";
		ShowEq(A, B, C, 0, 0, g_F, "u", "v");
		printf(" = 0.\n");
		if (abs(g_F) != 1) {
			printf(" Let F be the constant term.");
		}
		putchar('\n');
		UU = "U";
		VU = "V";
		UL = "u";
		VL = "v";
		FP = "F";
	}
	if (teach && T == 1) {
		UU = "X" + SCFstr;
		VU = "Y" + SCFstr;
		UL = "x" + SCFstr;
		VL = "y" + SCFstr;
		FP = "f" + SCFstr;
	}
	gcdAF = gcd(A, g_F);
	OrigA = A;
	OrigC = C;
	if (teach && gcdAF > 1) {
		std::cout << "Since gcd(A,F) = gcd(" << A << "," << g_F << ") = " << gcdAF
			<< " > 1, we have to replace y = ny' where n is a divisor of gcd(A,F).\n";
	}
	for (MagnifY = 1; MagnifY*MagnifY <= gcdAF; MagnifY++) {
		do {
			if (gcdAF / MagnifY*MagnifY != gcdAF) {
				continue;   // if MagnifY^2 is a not factor of gcdAF skip to next value 
			}
			if (teach) {
				if (ShowHR) {
					//w("<HR>"); 
				}
				else {
					ShowHR = true;
				}
			}
			MagnifY = gcdAF / MagnifY;
			g_F = H / T / T / MagnifY;
			ValF = abs(g_F);
			A = OrigA / MagnifY;
			C = OrigC*MagnifY;
			ValA = (A + ValF) % ValF;
			ValB = (B + ValF) % ValF;
			ValC = (C + ValF) % ValF;
			if (teach) {
				if (MagnifY != 1) {
					std::cout << "Let y = " << MagnifY << "y' We obtain: ";
				}
				putchar('\n');
				ShowEq(A, B, C, 0, 0, g_F, "x", "y'");
				printf(" = 0\n");
			}
			/* Find factors of F, store in array factors */
			NbrFactors = 0;
			Tmp = ValF;
			if (Tmp == 1) {
				factor[NbrFactors++] = 1;
			}
			else {
				while ((Tmp % 2) == 0) {
					factor[NbrFactors++] = 2;
					Tmp /= 2;
				}
				while ((Tmp % 3) == 0) {
					factor[NbrFactors++] = 3;
					Tmp /= 3;
				}
				s = 5;        /* Sequence of divisors 5, 7, 11, 13, 17, 19,... */
				do {
					while ((Tmp%s) == 0) {
						factor[NbrFactors++] = s;
						Tmp /= s;
					}
					s += 2;
					while ((Tmp%s) == 0) {
						factor[NbrFactors++] = s;
						Tmp /= s;
					}
					s += 4;
				} while (s*s <= Tmp);
				if (Tmp != 1) {
					factor[NbrFactors++] = Tmp;
				}
			}
			/* complete list of prime factors of F now in array F */

			mod[NbrFactors] = Tmp = 1;
			Pp = (2 * ValA) % ValF;
			for (index = NbrFactors - 1; index >= 0; index--) {
				P[index] = Pp;
				Tmp *= factor[index];
				mod[index] = Tmp;
				Pp = MultMod(MultMod(Pp, factor[index], ValF), factor[index], ValF);
			}
			Modulo = factor[NbrFactors - 1];  // get largest prime factor
			ValAM = (ValA + Modulo) % Modulo;
			ValBM = (ValB + Modulo) % Modulo;
			ValCM = (ValC + Modulo) % Modulo;
			if (ValAM == 0) {  /* Linear equation: sol=-C/B */
				Sol1 = Sol2 = MultMod(Modulo - ValCM, ModInv(ValBM, Modulo), Modulo);
			}
			else {    /* Quadratic equation Ax^2+Bx+C=0 (mod F) */
				if (Modulo>2) {
					Sol1 = MultMod(ValBM, ValBM, Modulo) - MultMod(4 * ValAM, ValCM, Modulo);
					if (Sol1<0) { Sol1 += Modulo; }
					/* Find square root of Sol1 mod Modulo */
					if (Sol1 == 0) {                 /* if double root: sol = -t/2a */
						Sol1 = Sol2 = MultMod(ModInv((2 * ValAM + Modulo) % Modulo, Modulo), ((-ValBM) + Modulo) % Modulo, Modulo);
					}
					else {
						if (ModPow(Sol1, (Modulo - 1) / 2, Modulo) == 1) { /* if sols exist */
							if (Modulo % 8 == 5) {
								s = ModPow(2 * Sol1, (Modulo - 5) / 8, Modulo);
								Sol1 = MultMod(MultMod(MultMod(MultMod(2 * Sol1, s, Modulo), s, Modulo) - 1, Sol1, Modulo), s, Modulo);
							}
							else {
								if (Modulo % 8 != 1) {
									Sol1 = ModPow(Sol1, (Modulo + 1) / 4, Modulo);
								}
								else {
									VarR = 1;
									VarQ = Modulo - 1;
									while (VarQ % 2 == 0) {
										VarQ /= 2;
										VarR *= 2;
									}
									VarX = 2;
									while (true) {
										VarY = ModPow(VarX, VarQ, Modulo);
										if (ModPow(VarY, VarR / 2, Modulo) != 1) { break; }
										VarX++;
									}
									VarX = ModPow(Sol1, (VarQ - 1) / 2, Modulo);
									VarV = MultMod(Sol1, VarX, Modulo);
									VarW = MultMod(VarV, VarX, Modulo);
									while (VarW != 1) {
										VarK = 1; VarD = VarW;
										while (VarD != 1) {
											VarD = MultMod(VarD, VarD, Modulo);
											VarK *= 2;
										}
										VarD = ModPow(VarY, VarR / VarK / 2, Modulo);
										VarY1 = MultMod(VarD, VarD, Modulo);
										VarR = VarK;
										VarV = MultMod(VarV, VarD, Modulo);
										VarW = MultMod(VarW, VarY1, Modulo);
										VarY = VarY1;
									}   /* end while */
									Sol1 = VarV;
								}     /* end modulo 8 = 1 */
							}
							s = ModInv((2 * ValAM) % Modulo, Modulo);
							Sol2 = MultMod((Modulo + Sol1 - ValBM) % Modulo, s, Modulo);
							Sol1 = MultMod((2 * Modulo - Sol1 - ValBM) % Modulo, s, Modulo);
						}
						else {   /* No solution exists */
							Sol1 = Sol2 = -1;
						}
					}
				}
				else {         /* Modulo <= 2 */
					if (Modulo == 2) {
						switch ((int)ValBM * 2 + (int)ValCM) {
						case 0:           /* A = 1, B = 0, C = 0 */
							Sol1 = Sol2 = 0;    /* Solution only for s=0 */
							break;
						case 1:           /* A = 1, B = 0, C = 1 */
							Sol1 = Sol2 = 1;    /* Solution only for s=1 */
							break;
						case 2:           /* A = 1, B = 1, C = 0 */
							Sol1 = 0;         /* Solution for s=0 and s=1 */
							Sol2 = 1;
							break;
						default:          /* A = 1, B = 1, C = 1 */
							Sol1 = Sol2 = -1;   /* No solutions */
							break;
						}
					}                   /* End Modulo = 2 */
					else {                /* Modulo = 1 */
						Sol1 = Sol2 = 0;
					}
				}
			}               /* End Quadratic Equation */

			//std::cout << "**temp SolContFrac: Sol1 =" << Sol1 << "  Sol2=" << Sol2 << "\n";

			ValAM = (ValA + ValF) % ValF;
			ValBM = (ValB + ValF) % ValF;
			ValCM = (ValC + ValF) % ValF;
			if (teach) {
				UL1 = UL;
				VL1 = VL + (MagnifY == 1 ? "" : "'");
				std::cout << "Let " << UL1 << " = s" << VL1 << " - " << FP << "z, so [-(as" << sq << " + bs + c)/"
					<< FP << "]" << VL1 << sq << " + (2as + b)" << VL1 << "z - a" << FP << "z" << sq << " = 1.\nSo \n";
				ShowEq(A, 0, 0, B, 0, C, "s", "");
				std::cout << " should be multiple of " << ValF << "\n";
			}

			NbrEqs = EqNbr = 0;
			int sol = 0;
			t = Sol1;
			/* if Sol1 >= 0 execute loop twice, otherwise not at all */
			for (cont = (Sol1<0 ? 2 : 0); cont<2; cont++) {
				index = NbrFactors - 1;
				v = mod[index];
				dif = 0;
				q = (MultMod((MultMod(ValAM, t, ValF) + ValBM) % ValF, t, ValF) + ValCM) % ValF;  /* q%v = 0 */
				while (true) {
					if (q%v == 0) {
						if (index == 0) {          /* Solution found */
							NbrEqs++;
							if (teach) {
								s = t*mod[1];
								for (index2 = 1; index2<NbrFactors; index2++) {
									s += pos[index2] * mod[index2 + 1];
								}
								s = s%ValF;
								if (sol == 0) {
									sol = 1;
									std::cout << "This holds for s = " << s;
								}
								else {
									std::cout << ", " << s;
								}
							}
							else {
								if (sol == 0)
									sol = 1;
							}
						}
						else {  /* solution not found */
							pos[index] = t;
							t = 0;
							for (index2 = index; index2<NbrFactors; index2++) {
								t += pos[index2] * mod[index2 + 1];
							}
							t = t%ValF;
							Dif[index] = dif;
							Q[index] = q;
							dif = MultMod((MultMod((2 * t + v) % ValF, ValAM, ValF) + ValBM) % ValF, v, ValF);
							Pp = P[--index];
							t = 0; v = mod[index];
							continue;
						}
					}

					if (index != NbrFactors - 1 && ++t < factor[index]) {
						q = (q + dif) % ValF;
						dif = (dif + Pp) % ValF;
						continue;
					}

					else {
						while (++index < NbrFactors) {
							t = pos[index];
							v = mod[index];   /* Restore previous values */
							if (index < NbrFactors - 1 && ++t < factor[index]) {
								Pp = P[index];
								dif = Dif[index];
								q = (Q[index] + dif) % ValF;
								dif = (dif + Pp) % ValF;
								break;  // exit from 'while' loop
							}
						}
						if (index >= NbrFactors) {
							break;
						}
					}
				}

				if (Sol1 == Sol2) {
					break;   /* Do not process double root */
				}
				t = Sol2;
			}
			if (teach) {
				if (sol == 0) {
					printf("No values of s makes the previous assertion true. \n");
					SqrtDisc = Tmp1;
					mpz_set(Dp_NUM, Dp_Tmp2);         // NUM = Tmp2
					mpz_set(Dp_DEN, Dp_Tmp3);         // DEN = Tmp3
					g_Disc = Tmp4;
					//std::cout << "**temp SolContFrac exit(1)\n";
					mpz_clears(Dp_A, Dp_B, Dp_C, Dp_R, Dp_S, Dp_T, Dp_Tmp2, Dp_Tmp3, Bi_Xcopy, Bi_Ycopy, NULL);
					assert(_CrtCheckMemory());
					return;
				}
				printf(". \n");
			}     /* end if (teach) */

			t = Sol1;
			/* if Sol1 >= 0 execute loop twice, otherwise not at all */
			for (cont = (Sol1<0 ? 2 : 0); cont<2; cont++) {
				index = NbrFactors - 1; v = mod[index];
				dif = 0;
				q = (MultMod((MultMod(ValAM, t, ValF) + ValBM) % ValF, t, ValF) + ValCM) % ValF;  /* q%v = 0 */
				while (true) {
					//std::cout << "**temp SolContFrac(2)\n";
					if (q%v == 0) {
						if (index == 0) {          /* Solution found */
							EqNbr++;
							s = t*mod[1];
							for (index2 = 1; index2<NbrFactors; index2++) {
								s += pos[index2] * mod[index2 + 1];
							}
							s = s%ValF;
							//  Calculate Dp_A = -((As+B)s+C)/F
							//            Dp_B = 2As+B
							//            Dp_C = -AF
							//            Dp_R = gcd(Dp_A, Dp_B, Dp_C)
							LongToDoublePrecLong(A, Dp_R);         // Dp_R = A
							LongToDoublePrecLong(s, Dp_S);         // Dp_S = s
							LongToDoublePrecLong(B, Dp_B);         // Dp_B = B
							Mult2LargeNumbers(Dp_R, Dp_S, Dp_C);     // Dp_C = As
							AddLarge(Dp_C, Dp_B, Dp_T);      // Dp_T = As+B
							Mult2LargeNumbers(Dp_T, Dp_S, Dp_A);     // Dp_A = (As+B)s
							LongToDoublePrecLong(C, Dp_R);         // Dp_R = C
							AddLarge(Dp_A, Dp_R, Dp_T);      // Dp_T = (As+B)s+C
							LongToDoublePrecLong(-g_F, Dp_R);        // Dp_R = -F
							DivideDoublePrecLong(Dp_T, Dp_R, Dp_A);   // Dp_A = -((As+B)s+C)/F
							AddLarge(Dp_C, Dp_C, Dp_C);      // Dp_C = 2As
							AddLarge(Dp_C, Dp_B, Dp_B);      // Dp_B = 2As+B
							LongToDoublePrecLong(A, Dp_T);         // Dp_T = A
							Mult2LargeNumbers(Dp_R, Dp_T, Dp_C);     // Dp_C = -AF
							gcd(Dp_A, Dp_B, Dp_T);      // Dp_T = gcd(Dp_A, Dp_B)
							gcd(Dp_T, Dp_C, Dp_R);      // Dp_R = gcd(Dp_A, Dp_B, Dp_C)
							if (teach) {
								std::cout << "Let s = " << s << ". Replacing in the above equation:\n";
								ShowBigEq(Dp_A, Dp_B, Dp_C, VL1, "z");
								printf(" = 1\n");
							}
							//if (Dp_R[0]>1 || Dp_R[1]>0 || Dp_R[2]>0 || Dp_R[3]>0 || Dp_R[4]>0 ||
							//	Dp_R[5]>0 && Dp_R[5]<DosALa31) {  // if Dp_R > 1 ...
							if (mpz_cmp_si(Dp_R, 1) >0) {  // if Dp_R > 1 ...
								if (teach) {
									std::cout << "Since the gcd of the three coefficients is "
										<< numToStr(Dp_R) <<
										" there are no integer solutions.\n";
								}
							}
							else {
								//LongArrayX = LongArrayY = NULL;
								GetRoot(Dp_A, Dp_B, Dp_C);
								if (SCFstr == "") {
									if (ContFrac(Dp_A, 3, 1, s, T, MagnifY, A)) {
										if (IsNeg(Bi_L2)) {
											ChangeSign(Bi_L1);
											ChangeSign(Bi_L2);
										}
										//std::cout << "**temp SolContFrac (2C) copy sol to ArrayX and ArrayY";
										mpz_set(Bi_Xcopy, Bi_L1);
										mpz_set(Bi_Ycopy, Bi_L2);
									}

									if (ContFrac(Dp_A, 3, (-1), s, T, MagnifY, A)) {
										if (IsNeg(Bi_L2)) {
											ChangeSign(Bi_L1);
											ChangeSign(Bi_L2);
										}
										if (teach) {
											NbrSols++;
											printf("Choosing the solution with minimum y:  \n%d solutions\n", NbrSols);

										}
										if (Compare(Bi_L2, Bi_Ycopy) > 0) {
											//std::cout << "**temp SolContFrac (2A)\n";
											ShowLargeXY("X", "Y", Bi_Xcopy, Bi_Ycopy, true, "", "");
										}
										else {
											//std::cout << "**temp SolContFrac (2B)\n";
											ShowLargeXY("X", "Y", Bi_L1, Bi_L2, true, "", "");
										}
										//w("</TABLE>");
									}
									else {
										if (!IsZero(Bi_Ycopy)) {
											//std::cout << "**temp SolContFrac(3)\n";
											ShowLargeXY("X", "Y", Bi_Xcopy, Bi_Ycopy, true, "", "");
										}
									}
								}
								else {
									ContFrac(Dp_A, 4, 1, s, T, MagnifY, A);
									ContFrac(Dp_A, 4, -1, s, T, MagnifY, A);
								}
							}
						}
						else {
							pos[index] = t;
							t = 0;
							for (index2 = index; index2<NbrFactors; index2++) {
								t += pos[index2] * mod[index2 + 1];
							}
							t = t%ValF;
							Dif[index] = dif;
							Q[index] = q;
							dif = MultMod((MultMod((2 * t + v) % ValF, ValAM, ValF) + ValBM) % ValF, v, ValF);
							Pp = P[--index];
							t = 0;
							v = mod[index];
							continue;
						}
					}
					if (index < NbrFactors - 1 && ++t < factor[index]) {
						q = (q + dif) % ValF;
						dif = (dif + Pp) % ValF;
						//std::cout << "**temp SolContFrac(4)\n";
						continue;
					}
					else {
						while (++index < NbrFactors) {
							//std::cout << "**temp SolContFrac(5)\n";
							t = pos[index]; v = mod[index];   /* Restore previous values */
							if (index < NbrFactors - 1 && ++t < factor[index]) {
								Pp = P[index];
								dif = Dif[index];
								q = (Q[index] + dif) % ValF;
								dif = (dif + Pp) % ValF;
								break;
							}
						}
						if (index >= NbrFactors) {
							//std::cout << "**temp SolContFrac(6)\n";
							break;   // break out of while loop
						}
					}
				}
				//std::cout << "**temp SolContFrac(7)\n";
				if (Sol1 == Sol2) {
					break;  /* Do not process double root */
				}
				//std::cout << "**temp SolContFrac(8)\n";
				t = Sol2;
			}
			if (teach) {
				putchar('\n');
			}
			SqrtDisc = Tmp1;
			mpz_set(Dp_NUM, Dp_Tmp2);        // NUM = Tmp2
			mpz_set(Dp_DEN, Dp_Tmp3);         // DEN = Tmp3
			g_Disc = Tmp4;
			//std::cout << "**temp SolContFrac(9)\n";
		} while (MagnifY*MagnifY > gcdAF);
		//std::cout << "**temp SolContFrac(10)\n";
	}   // end of for loop

		//std::cout << "**temp SolContFrac exit(11) \n";
	mpz_clears(Dp_A, Dp_B, Dp_C, Dp_R, Dp_S, Dp_T, Dp_Tmp2, Dp_Tmp3, Bi_Xcopy, Bi_Ycopy, NULL);
	assert(_CrtCheckMemory());
}

/* calculate base^exp mod (Mod)*/
long long ModPow(long long Base, long long Exp, long long Mod) {
	long long Pot, Pwr, mask, value;
	if (Exp == 0) { return 1LL; }

	assert(Exp >= 0);

	mask = 1LL;
	Pot = 1LL;
	Pwr = Base;
	value = 0;
	while (true) {
		if ((Exp & mask) != 0) {
			Pot = MultMod(Pot, Pwr, Mod);
			value += mask;
			if (value == Exp) { return Pot; }
		}
		mask *= 2LL;
		Pwr = MultMod(Pwr, Pwr, Mod);
	}
}

long long ModInv(long long Val, long long Mod) {
	long long U1, U3, V1, V3, Aux, Q;
	U1 = 1; U3 = Val; V1 = 0; V3 = Mod;
	while (V3 != 0) {
		Q = U3 / V3;
		Aux = U1 - V1*Q;
		U1 = V1;
		V1 = Aux;
		Aux = U3 - V3*Q;
		U3 = V3;
		V3 = Aux;
	}
	return (U1 + Mod) % Mod;
}

/* Calculate Num1*Num2 mod Mod.
NB will overflow and give wrong answer if Num1*Num2 > 64 bits!! */
long long MultModOld(long long Num1, long long Num2, long long Mod) {
	long long aux;
	aux = Num1*Num2 - Mod*(long long)(((double)Num1*(double)Num2) / (double)Mod);
	if (aux >= Mod) { 
		return aux - Mod; 
	}
	if (aux<0) { 
		return aux + Mod; 
	}
	return aux;
}

/* removes overflow risk. sign of result is same as sign of Num1*Num2 */
long long MultMod(long long Num1, long long Num2, long long Mod) {
	mpz_t N1, Prod;
	long long x;

	Mod = abs(Mod);   // ensure Mod is +ve (there is no mpz__tdiv_si function)
	mpz_init_set_ui(N1, Num1);
	mpz_init(Prod); 
	mpz_mul_si(Prod, N1, Num2);  // result = Num1*Num2
	x = mpz_tdiv_ui(Prod, Mod);  //  Prod%Mod is returned
	mpz_clears(N1, Prod, NULL);  // avoid memory leakage
	return x;
}
