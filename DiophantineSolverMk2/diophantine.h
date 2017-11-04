#pragma once

#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>  

/* constants */
const long long Bill = 1000000000;  // one billion = 10^9
const long long quintillion = Bill*Bill;  // one quintillion = 10^18
//const long long DosALa31 = (long long)1 << 31;
//const long long DosALa32 = (long long)1 << 32;
//const double dDosALa32 = (double)DosALa32;
//const long long DosALa32_1 = DosALa32 - 1;

/* functions that use bigints */
void mpz_add_si(mpz_t rop, const mpz_t op1, mpir_si op2);
void mpz_sub_si(mpz_t rop, const mpz_t op1, mpir_si op2);
bool IsZero(const mpz_t Bi_Nbr);
void ShowLargeXY(std::string x, std::string y, mpz_t Bi_H1, mpz_t Bi_K1,
	bool sol, const std::string eqX, const std::string eqY);
bool ShowHomoSols(int type, mpz_t Bi_H1, mpz_t Bi_K1, long long s, long long T,
	long long MagnifY, std::string eqX, std::string eqY);
void ShowLargeNumber(const mpz_t Bi_Nbr);
void ChangeSign(mpz_t Bi_Nbr);
void AddLarge(const mpz_t Bi_Nbr1, const mpz_t Bi_Nbr2, mpz_t Bi_Sum);
void AddLarge(const mpz_t Bi_Src, long long Nbr, mpz_t Bi_Dest);
//void AdjustSign(mpz_t Bi_Nbr);
void MultAddLargeNumbers(long long CPrev, const mpz_t Bi_Prev,
	long long CAct, const mpz_t Bi_Act, mpz_t Bi_Dest);
void MultLargeNumber(long long Coef, const mpz_t Bi_Nbr, mpz_t Bi_Dest);

long long DivLargeNumber(const mpz_t Bi_Nbr, long long Coef, mpz_t Bi_Dest);
long long tDivLargeNumber(const mpz_t n, long long d, mpz_t q);

/* other functions that need a separate declaration */
bool IsOne(const mpz_t Dp_A);
bool IsNeg(const mpz_t Dp_N1);
bool AreEqual(const mpz_t Dp_N1, const mpz_t  Dp_N2);
void ChSign(mpz_t Dp_Nbr);
void LongToDoublePrecLong(long long Nbr, mpz_t Dp_Out);
long long DoublePrecToLong(const mpz_t x);
void SubtDoublePrecLong(const mpz_t Dp_Nbr1, const mpz_t Dp_Nbr2, mpz_t Dp_Diff);

long long DivDoublePrec(const mpz_t Dp_Dividend, const mpz_t Dp_Divisor);
void Mult2LargeNumbers(const mpz_t Dp_Nbr1, const mpz_t Dp_Nbr2, mpz_t Dp_Prod);
void DivideDoublePrecLong(const mpz_t Dp_Dividend, const mpz_t Dp_Divisor, mpz_t Dp_Quotient);
std::string numToStr(long long num);
std::string par(long long num);
std::string numToStr(const mpz_t Dp_Nbr);
void ShowLin(long long D, long long E, long long F, std::string x, std::string y);
bool ContFrac(const mpz_t Dp_A, int type, int SqrtSign, long long s, long long T,
	long long MagnifY, long long A);
void ShowAllLargeSolutions();
void ShowEq(long long A, long long B, long long C, long long D, long long E, long long F,
	std::string x, std::string y);
void SolContFrac(long long H, long long T, long long A, long long B, long long C, std::string SCFstr);
void gcd(const mpz_t Dp_Nbr1, const mpz_t Dp_Nbr2, mpz_t Dp_Gcd);
long long gcd(long long M, long long N);
long long MultMod(long long factor1, long long factor2, long long Mod);
long long ModInv(long long Val, long long Mod);
long long ModPow(long long Base, long long Exp, long long Mod);
void ShowBigEq(mpz_t Dp_A, mpz_t Dp_B, mpz_t Dp_C, std::string x, std::string y);
void GetRoot(mpz_t Dp_A, mpz_t Dp_B, mpz_t Dp_C);
int Compare(const mpz_t Bi_array, const mpz_t Bi_K1);

