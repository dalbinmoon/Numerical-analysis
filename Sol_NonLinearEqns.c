#ifndef _SOL_NONLINEAREQUATION
#define _SOL_NONLINEAREQUATION

#include "basicHeader.h"
#include "UserDefine.h"
#include "Sol_NonLinearEqns.h"

void print_MethodOpt(OUT unsigned int *_METHOD) {

	printf("+++++++++++++++++++++++++++++++++\n");
	printf("+ Enter the number of Method    +\n");
	printf("+++++++++++++++++++++++++++++++++\n");
	printf("+ 1. BISECT - Bisection Method  +\n");
	printf("+ 2. NEWTON - Newton's Method   +\n");
	printf("+ 3. SECANT - Secant method     +\n");
	printf("+++++++++++++++++++++++++++++++++\n");

	scanf("%d", &Method);
}

float Diff0_Func(IN float x) {
	float y = pow(x, 3.0) - 3.0 * x + 1.0;
	return y;
}

float Diff1_Func(IN float x) {
	float y = 3.0 * pow(x, 2.0) - 3.0;
	return y;
}

void RootFinding(IN int Method) {
	switch (Method) {

	case BISECT:
		RootFinding_BISECT(lowerB, upperB);
		break;

	case NEWTON:
		RootFinding_NEWTON(initN);
		break;

	case SECANT:
		RootFinding_SECANT(initS0, initS1);
		break;

	default: 
		break;
	}
}

float RootFinding_BISECT(IN float lowerB, IN float upperB) {
	
	float ans1 = Diff0_Func(lowerB);
	float ans2 = Diff0_Func(upperB);

	int a = SIGN(ans1);
	int b = SIGN(ans2);

	if (a*b > 0) { State = ERR_DIVERGE; Error_Message(State); }
	
	for (; seq < MAX_ITER ; seq ++) {
		
		halfB = (lowerB + upperB) / 2;
		residual = Diff0_Func(halfB);
		bound = (upperB - lowerB) / 2;

		int a = SIGN(residual);
		int b = SIGN(ans1);

		Disp(BISECT, seq, lowerB, upperB, halfB, fabs(residual), bound);

		if (fabs(residual) < TOL_RESIDUAL) {
			State = CONVERGE; Error_Message(State); break;
		}

		if (a*b < 0 ) {
			upperB = halfB;
			ans2 = residual;
		}
		else {
			lowerB = halfB;
			ans1 = residual;
		}

		if (seq >= MAX_ITER) { State = ERR_ITER; Error_Message(State); }
	}
}

float RootFinding_NEWTON(IN float initN) {

	preN = initN - Diff0_Func(initN) / Diff1_Func(initN);

	for ( ; seq < MAX_ITER; seq++) {
		
		if(seq > 1) curN = preN - Diff0_Func(preN) / Diff1_Func(preN);

		Disp(NEWTON, seq, curN, Diff1_Func(curN), Diff0_Func(curN), 0.0, 0.0);

		if (fabs(curN - preN) < TOL_RESIDUAL) {
			State = CONVERGE; Error_Message(State); break;
		}

		preN = curN;
	}

	if (seq >= MAX_ITER) { State = ERR_ITER; Error_Message(State); }
}

float RootFinding_SECANT(IN float initS0, IN float initS1) {

	preS = initS1 - Diff0_Func(initS1) * (initS1 - initS0) / (Diff0_Func(initS1) - Diff0_Func(initS0));

	for (; seq < MAX_ITER; seq++) {

		nextS = curS - Diff0_Func(curS) * (curS - preS) / (Diff0_Func(curS) - Diff0_Func(preS));

		Disp(SECANT, seq, preS, curS, Diff0_Func(nextS), 0.0, 0.0);

		if (fabs(Diff0_Func(nextS)) < TOL_RESIDUAL) {
			State = CONVERGE; Error_Message(State); break;
		}

		preS = curS;
		curS = nextS;
		
	}
	
	if (seq >= MAX_ITER) { State = ERR_ITER; Error_Message(State); }
}

void Disp(IN int _Method, IN int seq, IN float disp1, IN float disp2, IN float disp3, IN float disp4, IN float disp5) {
	
	if (_Method == 1) {
		if (printflag == 0) printf("Step    LowerBound    UpperBound    Half Point    F(Half Point)    Bound Interval\n");
		printflag = 1;
		printf("%d        %f      %f      %f      %f      %f\n", seq, disp1, disp2, disp3, disp4, disp5);
	}

	if(_Method == 2) {
		if (printflag == 0) printf("Step        X            diff(X)        f(X)\n");
		printflag = 1;
		printf("%d        %f      %f      %f\n", seq, disp1, disp2, disp3);
	}

	if (_Method == 2) {
		if (printflag == 0) printf("Step        X            diff(X)        f(X)\n");
		printflag = 1;
		printf("%d        %f      %f      %f\n", seq, disp1, disp2, disp3);
	}

	if (_Method == 3) {
		if (printflag == 0) printf("Step       preX          curX           f(X)\n");
		printflag = 1;
		printf("%d        %f      %f      %f\n", seq, disp1, disp2, disp3);
	}
	
}

void Error_Message(IN unsigned int _Errortype) {
	if (_Errortype == 1) printf("Converge\n");
	if (_Errortype == 2) printf("Over Iterate\n");
	if (_Errortype == 3) printf("Diverge\n");
	if (_Errortype == 4) printf("Order is too high\n");
}
#endif