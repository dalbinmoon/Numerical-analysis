#ifndef _SOL_NONLINEAREQUATION
#define _SOL_NONLINEAREQUATION

#include "basicHeader.h"
#include "UserDefine.h"
#include "Sol_NonLinearEqns.h"
#include "Linear_Algebra.h"

void readFile(OUT float result[READ_ROW][READ_COLUMN]) {
	//printf("%f \n", RR[0][0]);
	FILE* myfile;
	myfile = fopen("hw3_3b.dat", "r");
	for (int i = 0; i < READ_ROW; i++) {
		fscanf(myfile, "%f %f %f", &result[i][0], &result[i][1], &result[i][2]);
	}

	fclose(myfile);
}


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

	case NEWTONRAPHON:
		RootFinding_NEWTONRAPHON(initNR);
		break;

	case RECURSIVE_LMS:
		Recursive_LMS(prex);
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
		
		curN = preN - Diff0_Func(preN) / Diff1_Func(preN);

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

float Recursive_LMS(OUT float out[READ_ROW]) {
	readFile(RR);

	//initalize and calulate usually used parameter
	for (int j = 0; j < READ_COLUMN; j++) {
		for (int k = 0; k < READ_COLUMN; k++) {
			preP[j][k] = initP[j][k];
			prex[j] = initx[j];
		}
	}
	
	float trans_H[READ_COLUMN][READ_COLUMN] = { 0, };
	transpose_matrix(H_LMS, trans_H);

	float inv_R[READ_COLUMN][READ_COLUMN] = { 0, };
	matrix_inverse_3(R_LMS, inv_R);

	float invR_transH[READ_COLUMN][READ_COLUMN] = { 0, };
	matrix_mul(trans_H, inv_R, invR_transH);

	float kf[READ_COLUMN][READ_COLUMN] = { 0, };

	for (int i = 0; i < 100; i++) {
		//calulate y
		for (int w = 0; w < READ_COLUMN; w++) {
			y[w] = sqrt(RR[i][w]) + minus_y[w];
		}

		float temp[READ_COLUMN][READ_COLUMN] = { 0, };

		//Gramian matrix update
		float invP[READ_COLUMN][READ_COLUMN] = { 0 };
		matrix_mul(invR_transH, H_LMS, temp);
		matrix_inverse_3(preP, invP);
		matrix_add(temp, invP, invcurP);

	
		//calculation of the kalman gain
		matrix_inverse_3(invcurP, curP);
		matrix_mul(curP, invR_transH, kf);

		//calculation of the residual sequence using sensor measurement vector y
		float temp1[READ_COLUMN] = { 0 };
		for (int i = 0; i < READ_COLUMN; i++) {
			for (int k = 0; k < READ_COLUMN; k++) {
				temp1[i] = temp1[i] + H_LMS[i][k] * prex[k];
			}
		}

		for (int j = 0; j < READ_COLUMN; j++) {
			y_tilda[j] = y[j] - temp1[j];
		}


		//update the LS solution
		float temp2[READ_COLUMN] = { 0 };
		for (int i = 0; i < READ_COLUMN; i++) {
			for (int k = 0; k < READ_COLUMN; k++) {
				temp2[i] = temp2[i] + kf[i][k] * y_tilda[k];
			}
		}

		for (int l = 0; l < READ_COLUMN; l++) {
			curx[l] = prex[l] + temp2[l];
		}

		//update value
		for (int j = 0; j < READ_COLUMN; j++) {
			prex[j] = curx[j];
			for (int k = 0; k < READ_COLUMN; k++) {
				preP[j][k] = curP[j][k];
				
			}
		}

		

	}

	for (int k = 0; k < READ_COLUMN; k++) {
		printf("%f ", curx[k]);
	}

	printf("\n");
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

	if (_Method == 4){
		if (printflag == 0) printf("Iteration # |\n");
		printflag = 1;
		printf("%2d          |   %2f\n", seq, disp1);
	}
}

void Error_Message(IN unsigned int _Errortype) {
	if (_Errortype == 1) printf("Converge\n");
	if (_Errortype == 2) printf("Over Iterate\n");
	if (_Errortype == 3) printf("Diverge\n");
	if (_Errortype == 4) printf("Order is too high\n");
	if (_Errortype == 5) printf("There is No Solution\n");
}
#endif